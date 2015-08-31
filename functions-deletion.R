library(epicalc)
library(dplyr)

CNVplot <- function(genetics_ranges)
{
  genetics_ranges <- filter(genetics_ranges, Chr_Gene == "22")

  genetics_ranges <- mutate(genetics_ranges, size = End - Start + 1)

  patients <- group_by(genetics_ranges, Patient.ID) %>% filter(Gain_Loss == "Loss" | Result.type == "mutation") %>% summarize(size = sum(size)) %>% arrange(desc(size))
  patients <- add_rownames(patients, var = "y")

  min.pos <- min(genetics_ranges$Start)
  max.pos <- max(genetics_ranges$End)
  max.size <- max(patients$size)

  height <- nrow(patients) * 10

  par(mar = c(6,0,0,0), oma = c(0,0,0,0), bty = "n")
  plot(0, xlim = c(min.pos - max.size / 10, max.pos + max.size / 10), ylim = c(-200, height), main = NULL, yaxt = "n", xaxt = "n", xlab = "", ylab = "")

  for (pat in patients$Patient.ID)
  {
    range <- filter(genetics_ranges, Patient.ID == pat)
    for (i in 1:nrow(range))
    {
      if (range$Result.type[i] == "coordinates" | range$Result.type[i] == "gene")
      {
        if (range$Gain_Loss[i] == "Gain")
          col <- "blue"
        else if (range$Gain_Loss[i] == "Loss")
          col <- "red"
      } else if (range$Result.type[i] == "mutation")
      {
        col <- "darkred"
        range$End[i] <- range$End[i]+50000
      }

      y <- as.numeric(patients$y[patients$Patient.ID == pat])

      if (nrow(range) > 1 & range$Result.type[i] != "gene" & i == 1)
        abline(h = height - y * 10 + 3)
      rect(range$Start[i], height - y * 10 + 7, range$End[i] + 1, height - y * 10, col = col, border = F)
    }
  }

  rect("50674641", -15, "50733212", height + 5, col = rgb(0,0,0,.1), border = F, lwd = 1)
}

delAnalysis <- function(genetics_ranges, data, depvar)
{
  genetics_ranges <- filter(genetics_ranges, Chr_Gene == "22", Gain_Loss == "Loss" | Result.type == "mutation")
  genetics_ranges <- mutate(genetics_ranges, size = End - Start + 1)

  patients <- group_by(genetics_ranges, Patient.ID) %>% summarize(size = sum(size)) %>% arrange(desc(size))
  patients <- add_rownames(patients, var = "y")

  df <- merge(data, patients)

  if (is.numeric(df[[depvar]]))
  {
    model <- lm(df[[depvar]] ~ df$size + df$Gender)
    table <- regress.display(model)$table[-1,]
  } else if (is.ordered(df[[depvar]]))
  {
    df2 <<- df
    depvar <<- depvar
    model <- polr(df2[[depvar]] ~ scale(df2$size) + df2$Gender + scale(df2$Age))
    table <- ordinal.or.display(model)
  } else
  {
    if (depvar == "Is.the.patient's.menstrual.cycle.regular?_ currently")
      model <- glm(df[[depvar]] ~ df$size + df$Age, family = binomial)
    else
      model <- glm(df[[depvar]] ~ df$size + df$Gender + df$Age, family = binomial)
    table <- logistic.display(model)$table
  }

  or  <- table[1, 1]
  icl <- table[1, 2]
  icu <- table[1, 3]
  p   <- table[1, 4]

  c(p,icl,or,icu)
}

delPlotGenes <- function(genetics_ranges, data, results_ranges, results_genes_p, results_lasso, genes, depvar, noOutput = T)
{
  #
  #genetics_ranges <- genetics_ranges
  #genes <- Genes
  #depvar <- vars$Variable[37]
  #
  genetics_ranges <- filter(genetics_ranges, Chr_Gene == "22", Gain_Loss == "Loss" | Result.type == "mutation")
  genetics_ranges <- mutate(genetics_ranges, size = End - Start + 1)

  results_genes_p_corr <- t(apply(results_genes_p[-(1:2)], 1, p.adjust, method = "bonferroni"))
  rownames(results_genes_p_corr) <- results_genes_p$Variable
  results_genes_p_corr <- -log(results_genes_p_corr, base = 10)

  patients <- group_by(genetics_ranges, Patient.ID) %>% summarize(size = sum(size)) %>% arrange(desc(size))
  patients <- add_rownames(patients, var = "y")

  min.pos <- min(genetics_ranges$Start)
  min.pos <- ifelse (min.pos < 40000000, 40000000, min.pos)
  max.pos <- max(genetics_ranges$End)
  max.size <- max(patients$size)

  df <- merge(data, patients)

  height <- nrow(patients) * 10

  if (!noOutput)
  {
    filename <- gsub(".", " ", depvar, fixed = T)
    filename <- gsub("_", " ", filename)
    filename <- gsub("?", " ", filename, fixed = T)
    filename <- gsub("/", " ", filename, fixed = T)
    png(filename = paste0("delplots/", substr(filename,1,128), ".png"), width = 1280, height = 1024)
  }

  main <- gsub(".", " ", depvar, fixed = T)
  main <- gsub("_", "\n", main)
  plot(0, xlim = c(min.pos, max.pos ), ylim = c(-(height/4 + 50), height), main = paste0(main, "\np = ", format(results_ranges$Range_p_corrected[results_ranges$Variable == depvar]), digits = 3), yaxt = "n", xaxt = "n", xlab = "", ylab = "")

  #Colors
  if (is.numeric(df[[depvar]]))
  {
    cols <- colorRampPalette(c("red", "blue"))(100)
  } else
    cols <- colorRampPalette(c("red", "blue"))(nlevels(df[[depvar]]))

  for (pat in df$Patient.ID)
  {
    range <- filter(genetics_ranges, Patient.ID == pat)
    for (i in 1:nrow(range))
    {
      if (is.numeric(df[[depvar]]))
      {
        col <- cols[1 + floor(100 * (df[[depvar]][df$Patient.ID == pat] - min(df[[depvar]], na.rm = T)) / max(df[[depvar]], na.rm = T))]
      } else
        col <- cols[unclass(df[[depvar]][df$Patient.ID == pat])]

      border <- F
      if (is.na(col))
      {
        col <- "white"
        #border <- "black"
      }
      y <- as.numeric(patients$y[patients$Patient.ID == pat])

      if (nrow(range) > 1 & range$Result.type[i] != "gene" & i == 1)
        abline(h = height - y * 10 + 3)
      rect(range$Start[i], height - y * 10 + 7, range$End[i] + 1, height - y * 10, col = col, border = border)
    }
  }

  maxp <- max(results_genes_p_corr[depvar,], na.rm = T)
  maxp <- ifelse(maxp < 2, 2, maxp)

  for (gene in genes$Gene)
  {
    if (results_lasso[results_lasso$Variable == depvar, gene] < 0)
    {
      col <- "red"
    } else if (results_lasso[results_lasso$Variable == depvar, gene] > 0)
    {
      col <- "green"
    } else
    {
      col <- "orange"
    }

    ygene <- results_genes_p_corr[depvar, gene] * (height/4)/maxp - (height/4 + 50)
    rect(genes$txStart[genes$Gene == gene], ygene, genes$txEnd[genes$Gene == gene], ygene + 10, col = col, border = F)
  }

  abline(h = -log(.05, base = 10) * (height/4)/maxp - (height/4 + 50), lty = 2)
  text(x = min.pos + 100, y = -log(.05, base = 10) * (height/4)/maxp - (height/4 + 50) + 20, label = "-log(p = .05)")

  rect("50674641", -15, "50733212", height + 5, col = rgb(0,0,0,.1), border = F, lwd = 1)
  legend("topleft", c(levels(df[[depvar]]), "Missing"), fill = c(cols, 0), title = "Legend")

  if (!noOutput)
    dev.off()
}

delPlotRange <- function(genetics_ranges, data, results_ranges, depvar, noOutput = T, bnw = F)
{
  #
  #genetics_ranges <- Genetics_ranges
  #depvar <- vars$Variable[37]
  #

  genetics_ranges <- filter(genetics_ranges, Chr_Gene == "22", Gain_Loss == "Loss" | Result.type == "mutation")
  genetics_ranges <- mutate(genetics_ranges, size = End - Start + 1)

  patients <- group_by(genetics_ranges, Patient.ID) %>% summarize(size = sum(size)) %>% arrange(desc(size))
  patients <- add_rownames(patients, var = "y")

  df <- merge(data, patients)

  height <- nrow(patients) * 10

  if (!noOutput)
  {
    filename <- gsub(".", " ", depvar, fixed = T)
    filename <- gsub("_", " ", filename)
    filename <- gsub("?", " ", filename, fixed = T)
    filename <- gsub("/", " ", filename, fixed = T)
    png(filename = paste0("delplotsrange/", substr(filename,1,128), ".png"), width = 400, height = 800)
  }

  main <- gsub(".", " ", depvar, fixed = T)
  main <- gsub("_", "\n", main)

  par(mar = c(6,0,0,0), oma = c(0,0,0,0), cex=.7)
  plot(NULL, xlim = c(0,10), ylim = c(-400, height), sub = paste0(main, "\np = ", format(results_ranges$Range_p_corrected[results_ranges$Variable == depvar], digits = 3)), yaxt = "n", xaxt = "n", xlab = "", ylab = "", asp = 1/50, bty = "n")

  #Colors
  if (bnw)
  {
    col1 <- rgb(.9,.9,.9)
    col2 <- rgb(0,0,0)
  } else {
    col1 <- "red"
    col2 <- "blue"
  }
  if (is.numeric(df[[depvar]]))
  {
    cols <- colorRampPalette(c(col1, col2))(100)
  } else
    cols <- colorRampPalette(c(col1, col2))(1 + nlevels(df[[depvar]]))

  for (pat in df$Patient.ID)
  {
    range <- filter(genetics_ranges, Patient.ID == pat)
    for (i in 1:nrow(range))
    {
      if (is.numeric(df[[depvar]]))
      {
        col <- cols[1 + floor(100 * (df[[depvar]][df$Patient.ID == pat] - min(df[[depvar]], na.rm = T)) / max(df[[depvar]], na.rm = T))]
      } else
        col <- cols[unclass(df[[depvar]][df$Patient.ID == pat])]

      border <- F
      if (is.na(col))
      {
        col <- "white"
        #border <- "black"
      }
      y <- as.numeric(patients$y[patients$Patient.ID == pat])
      rect(0, height - y * 10 + 7, 10, height - y * 10, col = col, border = border)
    }
  }

  legend("bottom", c(levels(df[[depvar]]), "Missing"), fill = c(cols[1:(nlevels(df[[depvar]]))], 0), title = NULL, bty = "n")

  if (!noOutput)
    dev.off()
}