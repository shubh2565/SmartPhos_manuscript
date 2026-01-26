###################### Function for processing the assay #######################
medianNorm <- function(x, method = "median") {
  if (method == "median") {
    # Calculate the median of each column, ignoring NA values
    mVal <- matrixStats::colMedians(x, na.rm = TRUE)
    # Adjust by the overall median of these medians
    mVal <- mVal - median(mVal, na.rm = TRUE)
  } else if (method == "mean") {
    # Calculate the mean of each column, ignoring NA values
    mVal <- colMeans(x, na.rm = TRUE)
    # Adjust by the overall mean of these means
    mVal <- mVal - mean(mVal, na.rm = TRUE)
  }
  mMat <- matrix(rep(mVal, each = nrow(x)), ncol =ncol(x))
  return(x-mMat)
}

getOneSymbol <- function(Gene) {
  # Apply a function to each element in the Gene vector
  outStr <- sapply(Gene, function(x) {
    # Split the string by semicolons
    sp <- str_split(x, ";")[[1]]
    # Return the last element in the split string
    sp[length(sp)]
  })
  # Remove names from the resulting vector
  names(outStr) <- NULL
  return(outStr)
}

preprocessPhos <- function(seData, filterList = NULL, missCut = 50,
                           transform = c("log2", "vst", "none"),
                           normalize = FALSE, getFP = FALSE,
                           removeOutlier = NULL, assayName = NULL, batch = NULL,
                           scaleFactorTab = NULL,
                           impute =  c("none", "QRILC", "MLE", "bpca",
                                       "missForest", "MinDet"),
                           verbose = FALSE) {
  
  transform <- match.arg(transform)
  impute <- match.arg(impute)
  
  # Retrieve the desired sample type or specified assay
  if (is.null(assayName)) {
    if (getFP) {
      # Retrieve FullProteome samples if specified
      ppe <- seData[,seData$sampleType %in% c("FullProteome", "FP")]
      colData(ppe) <- colData(seData)[colnames(ppe),]
    } else {
      # Otherwise, retrieve Phospho samples
      ppe <- seData[,seData$sampleType %in% c("Phospho", "PP")]
      colData(ppe) <- colData(seData)[colnames(ppe),]
    }
  } else {
    ppe <- seData[[assayName]]
    colData(ppe) <- colData(seData[,colnames(ppe)])
  }
  
  # Remove specified outliers
  if (length(removeOutlier) > 0) {
    if (length(removeOutlier) > 1) {
      for (i in removeOutlier) {
        ppe <- ppe[, !grepl(i, ppe$sample)]
      }
    }
    else {
      ppe <- ppe[, !grepl(removeOutlier, ppe$sample)]
    }
  }
  
  # Apply specified filters
  if (!is.null(filterList)) {
    for (n in names(filterList)) {
      ppe <- ppe[,ppe[[n]] %in% filterList[[n]]]
    }
  }
  
  # Rename columns to sample names
  colnames(ppe) <- ppe$sample
  # Get last gene name
  rowData(ppe)$Gene <- getOneSymbol(rowData(ppe)$Gene)
  # Get last phosphorylation site
  rowData(ppe)$Residue <- getOneSymbol(rowData(ppe)$Residue)
  rowData(ppe)$Position <- getOneSymbol(rowData(ppe)$Position)
  # Remove features without gene symbols
  ppe <- ppe[!rowData(ppe)$Gene %in% c(NA,""),]
  # Rename phosphorylation sites
  rowData(ppe)$site <- paste0(rowData(ppe)$Gene,"_",rowData(ppe)$Residue,
                              rowData(ppe)$Position)
  # Filter features based on missing values
  countMat <- assay(ppe)
  missPer <- rowSums(is.na(countMat))/ncol(countMat)*100
  ppeSub <- ppe[missPer < missCut,]
  
  # Apply transformation and normalization
  if (transform=="log2") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- medianNorm(log2(assay(ppeSub)))
      } else {
        assay(ppeSub) <- log2(t(t(
          assay(ppeSub))/scaleFactorTab[match(paste0(
            ppeSub$sample),scaleFactorTab$sample),]$scaleFactor))
      }
    } else {
      assay(ppeSub) <- log2(assay(ppeSub))
    }
  } else if (transform == "vst") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- vsn::justvsn(assay(ppeSub))
      } else {
        normMat <- t(t(
          assay(ppeSub))/scaleFactorTab[match(paste0(
            ppeSub$sample),scaleFactorTab$sample),]$scaleFactor)
        assay(ppeSub) <- vsn::justvsn(normMat, calib="none")
      }
    } else {
      assay(ppeSub) <- vsn::justvsn(assay(ppeSub), calib="none")
    }
  } else if (transform == "none") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- medianNorm(assay(ppeSub))
      } else {
        assay(ppeSub) <- t(t(
          assay(ppeSub))/scaleFactorTab[match(paste0(
            ppeSub$sample),scaleFactorTab$sample),]$scaleFactor)
      }
    } else {
      assay(ppeSub) <- assay(ppeSub)
    }
  }
  
  # Impute missing values
  if (impute != "none") {
    rowData(ppeSub)$name <- rowData(ppeSub)$site
    rowData(ppeSub)$ID <- rowData(ppeSub)$site
    if (impute == "missForest") {
      doParallel::registerDoParallel(cores = 6)  # Set number of CPU cores
      doRNG::registerDoRNG(seed = 123)
      mf <- missForest::missForest(t(assay(ppeSub)), parallelize = "forests",
                                   maxiter = 2, ntree = 50)
      imp <- t(mf$ximp)
    }
    else {
      imp <- DEP::impute(ppeSub, fun = impute)
    }
    assays(ppeSub)[["imputed"]] <- assay(imp)
    rowData(ppeSub)$name <- NULL
    rowData(ppeSub)$ID <- NULL
    # Show number of samples and features
    if (verbose) {
      message("Number of proteins and samples:")
      print(dim(ppeSub))
    }
  }
  
  # Remove batch effects if specified
  if(!is.null(batch)) {
    if(length(batch) == 1) {
      remBatchImp <- limma::removeBatchEffect(assays(ppeSub)[["imputed"]],
                                              batch = colData(ppeSub)[,batch])
      remBatch <- limma::removeBatchEffect(assay(ppeSub),
                                           batch = colData(ppeSub)[,batch])
    }
    else {
      remBatchImp <- limma::removeBatchEffect(assays(ppeSub)[["imputed"]],
                                              batch = colData(
                                                ppeSub)[,batch[1]],
                                              batch2 = colData(
                                                ppeSub)[,batch[2]])
      remBatch <- limma::removeBatchEffect(assay(ppeSub),
                                           batch = colData(ppeSub)[,batch[1]],
                                           batch2 = colData(ppeSub)[,batch[2]])
    }
    assays(ppeSub)[["imputed"]] <- assay(remBatchImp)
    assay(ppeSub) <- assay(remBatch)
  }
  
  return(ppeSub)
}

##################### Plot for the completeness of the assay ###################
plotCompleteness <- function(se, colorByCol = "none", 
                             title = "Percentage of sample completeness") {
  # Extract the assay data from the SummarizedExperiment object
  countMat <- assay(se)
  
  # Create a table with sample names and their corresponding percentage of
  # non-missing values
  plotTab <- tibble(
    sample = se$sample,
    perNA = colSums(is.na(countMat)) / nrow(countMat)
  )
  
  # Extract metadata from the SummarizedExperiment object
  meta <- as.data.frame(colData(se))
  # Join the count data with metadata
  plotTab <- left_join(plotTab, meta, by = "sample")
  
  # Generate the bar plot using ggplot2
  g <- ggplot(plotTab, aes(x = sample, y = 1 - perNA)) +
    ggtitle(title) +
    ylab("completeness") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Add color to the barplot if a valid metadata column is specified
  if (colorByCol == "none") {
    g <- g + geom_bar(stat = "identity")
  } else {
    g <- g + geom_bar(stat = "identity", aes(fill = !!sym(colorByCol)))
  }
  
  return(g)
}

################### Plot of total number of sites detected #####################
plotTotalSites <- function(se, colorByCol = "none",
                           title = "Total number of phosphosites") {
  # Extract the assay data from the SummarizedExperiment object
  countMat <- assay(se)
  
  # Create a table with sample names and their corresponding percentage of
  # non-missing values
  plotTab <- tibble(
    sample = se$sample,
    total = colSums(!is.na(countMat)) 
  )
  
  # Extract metadata from the SummarizedExperiment object
  meta <- as.data.frame(colData(se))
  # Join the count data with metadata
  plotTab <- left_join(plotTab, meta, by = "sample")
  
  # Generate the bar plot using ggplot2
  g <- ggplot(plotTab, aes(x = sample, y = total)) +
    ggtitle(title) +
    ylab("Total sites") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Add color to the barplot if a valid metadata column is specified
  if (colorByCol == "none") {
    g <- g + geom_bar(stat = "identity")
  } else {
    g <- g + geom_bar(stat = "identity", aes(fill = !!sym(colorByCol)))
  }
  
  return(g)
}


################################# PCA Plot #####################################

plotPCAnalysis <- function(pca, se, xaxis = "PC1", yaxis = "PC2", color = "none",
                    shape = "none",
                    title = "Principal Components") {
  # Calculate the proportion of variance explained by each principal component
  varExplained <- pca$sdev^2 / sum(pca$sdev^2)
  # Convert the PCA result to a data frame
  pcaDf <- as.data.frame(pca[["x"]])
  # Convert the metadata to a data frame
  meta <- as.data.frame(colData(se))
  # Join the PCA scores with the metadata
  pcaMeta <- left_join(rownames_to_column(pcaDf),
                       meta,
                       by = c("rowname" = "sample")
  )
  
  # Create the initial ggplot object with labels for variance explained
  g <- ggplot(pcaMeta, aes(
    x = !!sym(xaxis), y = !!sym(yaxis),
    text = paste("sample:", meta$sample)
  )) +
    theme_bw() +
    theme(legend.position = "top") +
    labs(
      x = paste0(
        xaxis, ": ",
        round(varExplained[as.numeric(strsplit(xaxis, "PC")[[1]][2])] * 100,
              1), "%"),
      y = paste0(
        yaxis, ": ",
        round(varExplained[as.numeric(strsplit(yaxis, "PC")[[1]][2])] * 100,
              1), "%")
    ) +
    scale_shape(solid = FALSE) +
    ggtitle(title)
  
  # Add points to the plot with optional color and shape aesthetics
  if (color == "none" & shape == "none") {
    g <- g + geom_point(size = 3)
  } else if (color == "none") {
    g <- g + geom_point(aes(shape = !!sym(shape)), size = 3, alpha = 0.75)
  } else if (shape == "none") {
    g <- g + geom_point(aes(color = !!sym(color)), size = 3, alpha = 0.75, 
                        stroke = 0.8)
  } else {
    g <- g + geom_point(
      aes(
        color = !!sym(color),
        shape = !!sym(shape)
      ),
      size = 3, alpha = 0.75, stroke = 0.8
    )
  }
  
  return(g)
}

############################ Intensity Boxplot DEA #############################

intensityBoxPlotDEA <- function(se, id, symbol) {
  
  group <- value <- subjectID <- NULL
  
  exprMat <- assays(se)[["Intensity"]]
  
  # Check if the SE object contains subject-specific data
  if (is.null(se$subjectID)) {
    # Prepare data frame for plotting without subject-specific information
    plotTab <- data.frame(
      group = se$comparison,
      value = exprMat[id, ]
    )
    p <- ggplot(plotTab, aes(x = group, y = value))
  } else {
    # Prepare data frame for plotting with subject-specific information
    plotTab <- data.frame(
      group = se$comparison,
      value = exprMat[id, ],
      subjectID = se$subjectID
    )
    p <- ggplot(plotTab, aes(x = group, y = value)) +
      geom_line(aes(group = subjectID), linetype = "dotted", color = "grey50")
  }
  
  # Create the boxplot with additional formatting
  p <- p + geom_boxplot(aes(fill = group),
                        width = 0.5, alpha = 0.5,
                        outlier.shape = NA
  ) +
    geom_point() +
    ylab("Normalized Intensities") + xlab("") +
    ggtitle(symbol) + theme_bw() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 15)
    )
  
  return(p)
}

################## Volcano plot for Differential Expression ####################

plotVolcanoDEA <- function(resDE, pFilter = 0.05, fcFilter = 0.5, 
                           title = "Volcano plot", usePadj = FALSE) {
  # Convert the input table to a data frame and ensure the 'ID' column is of
  # type character
  dataVolcano <- data.frame(resDE)
  dataVolcano$ID <- as.character(dataVolcano$ID)
  # Categorize each gene based on the provided p-value and log2 fold-change
  # thresholds
  if (usePadj) {
    dataVolcano <- mutate(dataVolcano, expression = case_when(
      dataVolcano$log2FC >= as.numeric(fcFilter) &
        dataVolcano$padj <= as.numeric(pFilter) ~ "Up",
      dataVolcano$log2FC <= -as.numeric(fcFilter) &
        dataVolcano$padj <= as.numeric(pFilter) ~ "Down",
      dataVolcano$padj > as.numeric(pFilter) |
        (dataVolcano$log2FC < as.numeric(fcFilter) &
           dataVolcano$log2FC > -as.numeric(fcFilter)) ~ "Not Sig"
    ))
  }
  else {
    dataVolcano <- mutate(dataVolcano, expression = case_when(
      dataVolcano$log2FC >= as.numeric(fcFilter) &
        dataVolcano$pvalue <= as.numeric(pFilter) ~ "Up",
      dataVolcano$log2FC <= -as.numeric(fcFilter) &
        dataVolcano$pvalue <= as.numeric(pFilter) ~ "Down",
      dataVolcano$pvalue > as.numeric(pFilter) |
        (dataVolcano$log2FC < as.numeric(fcFilter) &
           dataVolcano$log2FC > -as.numeric(fcFilter)) ~ "Not Sig"
    ))
  }
  
  
  counts <- table(dataVolcano$expression)
  
  # redefine legend labels with counts
  labels_with_counts <- c(
    paste0("Up (", counts["Up"], ")"),
    paste0("Down (", counts["Down"], ")"),
    paste0("Not Sig (", counts["Not Sig"], ")")
  )
  
  # find top regulated genes/sites
  # top_up <- dataVolcano %>%
  #   filter(expression == "Up") %>%
  #   arrange(pvalue) %>%
  #   head(5)
  # 
  # top_down <- dataVolcano %>%
  #   filter(expression == "Down") %>%
  #   arrange(pvalue) %>%
  #   head(5)
  # 
  # top_labels <- bind_rows(top_up, top_down)
  # print(top_labels)
  
  # Create the volcano plot
  if (usePadj) {
    p <- ggplot(dataVolcano, aes(x = log2FC, y = -log10(padj)))
  }
  else {
    p <- ggplot(dataVolcano, aes(x = log2FC, y = -log10(pvalue))) 
  }
  v <- p +
    # Plot the points and color them based on their expression status
    geom_point(aes(color = expression, fill = expression), shape = 21, size = 1.5) +
    scale_color_manual(values = c("Up" = "firebrick3",
                                  "Down" = "navy", "Not Sig" = "darkgrey"),
                       breaks = c("Up", "Down", "Not Sig"),  # order of legend items
                       labels = labels_with_counts,
                       name   = "Significance") +
    scale_fill_manual(values = c("Up" = "#FE8484", "Down" = "#5494DA", "Not Sig" = "lightgrey"),
                      breaks = c("Up", "Down", "Not Sig"),  # order of legend items
                      labels = labels_with_counts,
                      name   = "Significance")+
    # Add vertical lines for fold-change thresholds
    geom_vline(xintercept = 0, color = "brown", linetype = "solid",
               linewidth = 0.25) +
    geom_vline(xintercept = as.numeric(fcFilter), color = "brown",
               linetype = "dashed") +
    geom_vline(xintercept = -as.numeric(fcFilter), color = "brown",
               linetype = "dashed") +
    annotate(x = 3.0, y = -log10(as.numeric(pFilter)) - 0.15,
             label = paste("P-value = ", as.numeric(pFilter)),
             geom = "text", size = 4, color = "brown") +
    # Add horizontal lines for p-value thresholds
    geom_hline(yintercept = -log10(as.numeric(pFilter)), color = "brown",
               linetype = "dashed") +
    # geom_text_repel(
    #   data = top_labels,
    #   aes(label = site),
    #   size = 3,
    #   max.overlaps = 10
    # ) +
    
    # annotate("text", 
    #          x = max(dataVolcano$log2FC, na.rm = TRUE) * 0.9, 
    #          y = max(-log10(dataVolcano$pvalue), na.rm = TRUE) * 0.95, 
    #          label = paste0("Up: ", counts["Up"]), 
    #          color = "firebrick", hjust = 1) +
    # annotate("text", 
    #          x = min(dataVolcano$log2FC, na.rm = TRUE) * 0.9, 
    #          y = max(-log10(dataVolcano$pvalue), na.rm = TRUE) * 0.95, 
    #          label = paste0("Down: ", counts["Down"]), 
    #          color = "navyblue", hjust = 0) +
    xlab("absolute log2(Quantity) difference") +
    ggtitle(title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(v)
}






################### Combine different GMT files into one #######################

combineGMTfiles <- function(path) {
  # List all .gmt files in the folder
  gmt_files <- list.files(path, pattern = "\\.gmt$", full.names = TRUE)
  
  # Load all GMT files into a list
  gsc_list <- lapply(gmt_files, loadGSC)
  
  # Combine gsc and addInfo from all files
  combined_gsc <- list(
    gsc = do.call(c, lapply(gsc_list, function(x) x$gsc)),
    addInfo = do.call(c, lapply(gsc_list, function(x) x$addInfo))
  )
  
  # Ensure gene set names are unique to avoid overwriting
  names(combined_gsc$gsc) <- make.unique(names(combined_gsc$gsc))
  
  # Build GMT lines as character vector
  gmt_lines <- vapply(seq_along(combined_gsc$gsc), function(i) {
    paste(c(names(combined_gsc$gsc)[i],
            combined_gsc$addInfo[[i]],
            combined_gsc$gsc[[i]]),
          collapse = "\t")
  }, FUN.VALUE = character(1))
  
  writeLines(gmt_lines, "combinedGeneSet.gmt")
  cat("Combined", length(gmt_files), "GMT files into 'combinedGeneSet.gmt'\n")
}