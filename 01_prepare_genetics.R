library(tidyverse)

# liftOver can produce multiple overlapping output.
# This function merges ranges to their simplest form possible
merge_overlapping_genetic_ranges <- function(ranges)
{
  ranges %>%
    gather(end, pos, hg38.Start, hg38.End) %>%
    arrange(pos) %>%
    mutate(depth = cumsum(end == "hg38.Start") -
           cumsum(end == "hg38.End") %>%
           lag %>%
           replace_na(0)) %>%
    filter(depth == 1) %>%
    mutate(num = cumsum(end == "hg38.Start")) %>%
    spread(end, pos) %>%
    select(rowname, num, hg38.Start, hg38.End)
}

# Write the liftOver.bed file on disk, for the provided data
write_bedfile <- function(data)
{
  data %>%
    mutate(Chr = "chr22",
           id = str_c(Patient.ID, "_", Array.num)) %>%
  select(Chr, Array.Start, Array.End, id) %>%
  write.table(file = "liftOver.bed" , sep = "\t", row.names = F, col.names = F, quote = F)
}

# Stand-alone version of liftOver
liftOver_alone <- function(build, data)
{
  chainfile = str_c(build, "ToHg38.over.chain")

# First try with unique mapping
  data %>%
    write_bedfile

  system2("./liftOver", c("liftOver.bed", chainfile, "liftOver.out", "liftOver.unmap"))

  read_tsv("liftOver.out", col_names = c("Chr", "hg38.Start", "hg38.End", "rowname")) -> out

  if (nrow(out) > 0)
  {
    out %>%
      separate(rowname, into = c("Patient.ID", "Array.num"), sep = "_") %>%
      mutate_at(vars(Patient.ID, Array.num), as.numeric) %>%
      select(-Chr) %>%
      set_names(~ str_c("Pass1.", .x)) -> out
  }

  data %>%
    left_join(out, by = c("Patient.ID" = "Pass1.Patient.ID",
                          "Array.num" = "Pass1.Array.num")) -> result

# Second liftOver pass with multiple alignment and auto merging
  read_tsv("liftOver.unmap", col_names = c("Chr", "Array.Start", "Array.End", "rowname"), comment = "#") -> unmap

  if (nrow(unmap) > 0)
  {
    data %>%
      anti_join(out, by = c("Patient.ID" = "Pass1.Patient.ID",
                            "Array.num" = "Pass1.Array.num")) %>%
      write_bedfile

    system2("./liftOver", c("liftOver.bed", chainfile, "liftOver.out", "liftOver.unmap", "-multiple"))

    read_tsv("liftOver.out", col_names = c("Chr", "hg38.Start", "hg38.End", "rowname", "num")) -> outunmap
    if (nrow(outunmap) > 0)
    {
      outunmap %>%
        group_by(rowname) %>%
        do(merge_overlapping_genetic_ranges(.)) %>%
        separate(rowname, into = c("Patient.ID", "Array.num"), sep = "_") %>%
        mutate_at(vars(Patient.ID, Array.num), as.numeric) %>%
        set_names(~ str_c("Pass2.", .x)) -> out

      result %>%
        left_join(out, by = c("Patient.ID" = "Pass2.Patient.ID",
                              "Array.num" = "Pass2.Array.num")) -> result
    }
  }

  result
}

# Execute liftOver in two passes.
# - with unique mappings
# - with multiple outputs on the firt failed ones, and merging
liftOver <- function(build, data)
{
  chainfile = str_c(build, "ToHg38.over.chain")

# First try with unique mapping
  data %>%
    write_bedfile

  system2("./liftOver", c("liftOver.bed", chainfile, "liftOver.out", "liftOver.unmap"))

  read_tsv("liftOver.out", col_names = c("Chr", "hg38.Start", "hg38.End", "rowname")) -> out

  if (length(out) > 0)
  {
    out %>%
      separate(rowname, into = c("Patient.ID", "Array.num"), sep = "_") %>%
      mutate_at(vars(Patient.ID, Array.num), as.numeric) %>%
      select(-Chr) -> out
  }

# Second liftOver pass with multiple alignment and auto merging
  read_tsv("liftOver.unmap", col_names = c("Chr", "Array.Start", "Array.End", "rowname"), comment = "#") -> unmap

  if (length(unmap) > 0)
  {
    data %>%
      anti_join(out) %>%
      write_bedfile

    system2("./liftOver", c("liftOver.bed", chainfile, "liftOver.out", "liftOver.unmap", "-multiple"))

    read_tsv("liftOver.out", col_names = c("Chr", "hg38.Start", "hg38.End", "rowname", "num")) -> outunmap
    if (length(outunmap) > 0)
    {
      outunmap %>%
        group_by(rowname) %>%
        do(merge_overlapping_genetic_ranges(.)) %>%
        separate(rowname, into = c("Patient.ID", "Array.num"), sep = "_") %>%
        mutate_at(vars(Patient.ID, Array.num), as.numeric) %>%
        bind_rows(out) -> out
    }
  }

  out
}

# Process each row according to the result type
process <- function(type, data)
{
  if (type == "Array")
  {
    data %>%
      group_by(Array.BrowserBuild) %>%
      nest %>%
      mutate(liftover = map2(Array.BrowserBuild, data, liftOver),
             data2 = map2(data, liftover, left_join)) %>%
      unnest(data2) %>%
      # Compare with pre-filled mappings and replace
      mutate(startok = Array.Start.latest == hg38.Start,
             endok = Array.End.latest == hg38.End) %>%
      arrange(Patient.ID, Array.num) %>%
      mutate(Array.Start.latest = case_when(Array.Start.latest %>% is.na & hg38.Start %>% is.na ~ NA_real_,
                                            Array.Start.latest %>% is.na ~ hg38.Start,
                                            hg38.Start %>% is.na ~ Array.Start.latest,
                                            T ~ hg38.Start),
             Array.End.latest = case_when(Array.End.latest %>% is.na & hg38.End %>% is.na ~ NA_real_,
                                          Array.End.latest %>% is.na ~ hg38.End,
                                          hg38.End %>% is.na ~ Array.End.latest,
                                          T ~ hg38.End)) %>%
      select(-hg38.Start, -hg38.End, -startok, -endok) %>%
      mutate(num = (num - 1) %>% replace_na(0),
             Array.num = Array.num + num) %>%
      select(-num)
  } else
  {
    data
  }
}

# Verify processing of liftOver
read_csv("data/genetic.csv") %>%
  mutate(Array.BrowserBuild = Array.BrowserBuild %>% fct_recode(hg17 = "NCBI35/hg17",
                                                                hg18 = "NCBI36/hg18",
                                                                hg18 = "GRCh36/hg18",
                                                                hg19 = "GRCh37/hg19")) %>%
  filter(Result.type == "Array") %>%
  filter(Patient.ID == 10131) %>%
  select(Patient.ID, Array.num, Array.Type, Array.BrowserBuild, Array.Start, Array.End) %>%
  group_nest(Array.BrowserBuild) %>%
  mutate(data = map2(Array.BrowserBuild, data, liftOver_alone)) -> liftovered

liftovered %>% unnest %>% print(n = Inf)

# Process all the genetic results
read_csv("data/genetic.csv") %>%
  mutate(Array.BrowserBuild = Array.BrowserBuild %>% fct_recode(hg17 = "NCBI35/hg17",
                                                                hg18 = "NCBI36/hg18",
                                                                hg18 = "GRCh36/hg18",
                                                                hg19 = "GRCh37/hg19")) %>%
  group_by(Result.type) %>%
  nest %>%
  mutate(data = map2(Result.type, data, process)) %>%
  unnest -> genetic

# Write the processed genetic results
genetic %>%
  write_csv("data/genetic_hg38.csv")


# # Plot check
# genetic %>%
#   filter(Result.type == "Array") %>%
#   mutate(Patient.ID = Patient.ID %>% factor %>% fct_reorder(Array.Start, min)) %>%
#   ggplot() +
#   aes(x = Patient.ID, ymin = Array.Start, ymax = Array.End, color = Array.Type) +
#   geom_linerange() +
#   coord_flip()
