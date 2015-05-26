library(dplyr)

CNVplot <- function(genetics_ranges)
{
  genetics_ranges <- filter(genetics_ranges, Genome.Browser.Build == "GRCh38/hg38", Chr.Gene == "22")

  genetics_ranges$Start <- as.numeric(genetics_ranges$Start)
  genetics_ranges$End <- as.numeric(genetics_ranges$End)
  genetics_ranges <- mutate(genetics_ranges, size = End - Start + 1)

  patients <- group_by(genetics_ranges, Patient.ID) %>% summarize(size = sum(size)) %>% arrange(desc(size))
  patients <- add_rownames(patients, var = "y")

  min.pos <- min(genetics_ranges$Start)
  max.pos <- max(genetics_ranges$End)
  max.size <- max(patients$size)

  height <- nrow(patients) * 10

  plot(0, xlim = c(min.pos - max.size / 10, max.pos + max.size / 10), ylim = c(0, height), main = "Chromosomal abnormalities", yaxt = "n", xaxt = "n", xlab = "", ylab = "")

  for (pat in patients$Patient.ID)
  {
    range <- filter(genetics_ranges, Patient.ID == pat)
    for (i in 1:nrow(range))
    {
      if (range$Result.type[i] == "coordinates" | range$Result.type[i] == "gene")
      {
        if (range$Gain.Loss[i] == "Gain")
          col <- "blue"
        else if (range$Gain.Loss[i] == "Loss")
          col <- "red"
      }
      else if (range$Result.type[i] == "mutation")
        col <- "darkred"

      y <- as.numeric(patients$y[patients$Patient.ID == pat])

      if (nrow(range) > 1 & range$Result.type[i] != "gene" & i == 1)
        abline(h = height - y * 10 + 3)
      rect(range$Start[i], height - y * 10 + 7, range$End[i] + 1, height - y * 10, col = col, border = col)
    }
  }

  rect("50674641", -15, "50733212", height + 5, col = rgb(0,0,0,.1), border = F, lwd = 1)
}

delPlot <- function(genetics_ranges, depvar, noOutput = T)
{
  genetics_ranges <- filter(genetics_ranges, Genome.Browser.Build == "GRCh38/hg38", Chr.Gene == "22", Gain.Loss == "Loss")

  genetics_ranges$Start <- as.numeric(genetics_ranges$Start)
  genetics_ranges$End <- as.numeric(genetics_ranges$End)
  genetics_ranges <- mutate(genetics_ranges, size = End - Start + 1)

  patients <- group_by(genetics_ranges, Patient.ID) %>% summarize(size = sum(size)) %>% arrange(desc(size))
  patients <- add_rownames(patients, var = "y")

  depvar <- filter(depvar, Patient.ID %in% patients$Patient.ID)
  depvar[2] <- factor(depvar[[2]], exclude = "")

  min.pos <- min(genetics_ranges$Start)
  max.pos <- max(genetics_ranges$End)
  max.size <- max(patients$size)

  df <- merge(depvar,patients)

  height <- nrow(patients) * 10

  p <- 1
  tryCatch(p <- kruskal.test(df$size ~ df[[2]])$p.value,
           error = function(e) e)

  if (!noOutput)
  {
    filename <- gsub(".", " ", names(depvar[2]), fixed = T)
    filename <- gsub("_", " ", filename)
    filename <- gsub("?", " ", filename, fixed = T)
    filename <- gsub("/", " ", filename, fixed = T)
    png(filename = paste0("delplots/", substr(filename,1,128), ".png"), width = 1280, height = 1024)
  }

  main <- gsub(".", " ", names(depvar[2]), fixed = T)
  main <- gsub("_", "\n", main)
  plot(0, xlim = c(min.pos - max.size / 10, max.pos + max.size / 10), ylim = c(0, height), main = paste0(main, "\np = ", p), yaxt = "n", xaxt = "n", xlab = "", ylab = "")

  for (pat in df$Patient.ID)
  {
    range <- filter(genetics_ranges, Patient.ID == pat)
    for (i in 1:nrow(range))
    {

      col <- as.numeric(df[[2]][df$Patient.ID == pat]) + 1
      border <- F
      if (is.na(col))
      {
        col <- "white"
        border <- "black"
      }
      y <- as.numeric(patients$y[patients$Patient.ID == pat])

      if (nrow(range) > 1 & range$Result.type[i] != "gene" & i == 1)
        abline(h = height - y * 10 + 3)
      rect(range$Start[i], height - y * 10 + 7, range$End[i] + 1, height - y * 10, col = col, border = border)
    }
  }

  rect("50674641", -15, "50733212", height + 5, col = rgb(0,0,0,.1), border = F, lwd = 1)
  legend("topleft", c(levels(depvar[[2]]), "Missing"), fill = c(2:(nlevels(depvar[[2]]) + 1), 0), title = "Legend")

  if (!noOutput)
    dev.off()

  p
}