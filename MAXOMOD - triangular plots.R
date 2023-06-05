################ load libraries ###############
suppressMessages( library(ggplot2))
suppressMessages( library(gplots))
suppressMessages( library(RColorBrewer))
suppressMessages( library(dplyr))
suppressMessages( library(readxl))
suppressMessages( library(pheatmap))
suppressMessages( library(circlize))
suppressMessages( library(tools))        # for file path, basename_without_extension(file_path_sans_ext) and extensions(file_ext)
suppressPackageStartupMessages(library("argparse"))
BiocManager::install("BHC")
suppressPackageStartupMessages(library("BHC"))
suppressPackageStartupMessages(library("dendsort"))
suppressPackageStartupMessages(library("coop"))
suppressPackageStartupMessages(library("effects"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("gtable"))
suppressPackageStartupMessages(library("stats"))
suppressPackageStartupMessages(library("grDevices"))
suppressPackageStartupMessages(library("graphics"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ComplexHeatmap"))

########### Get data ###########

# -> proteomics
proteomics_data <- read_excel("~/Documents/HMGU/MAXOMOD/data to be used/Proteomics_final_table.xls")
colnames(proteomics_data) <- c("Protein ID",colnames(proteomics_data[2:ncol(proteomics_data)]))
protein_names <- proteomics_data[,c(1,105)]
count_proteomics <- proteomics_data[,2:104]
input_proteomics <- as.data.frame(t(count_proteomics))
input_proteomics$Sample <- rownames(input_proteomics)
clinical_info <- read_excel("~/Documents/HMGU/MAXOMOD/data to be used/Updated list for MAXOMOD CSF Proteomics and Metabolomics_2023.04.28 (in use).xlsx")
clinical_info <- clinical_info[1:51,]

## -> some pre-processing of NA data not coded as NA 
clinical_info$`pNFh (pg/ml)` <- ifelse(clinical_info$`pNFh (pg/ml)` == "na",NA,clinical_info$`pNFh (pg/ml)`)
clinical_info$`ALS progresion rate per month (delta ALS-FRSr / days *30)` <- ifelse(clinical_info$`ALS progresion rate per month (delta ALS-FRSr / days *30)` == "na",
                                                                                    NA,
                                                                                    clinical_info$`ALS progresion rate per month (delta ALS-FRSr / days *30)`)
clinical_info$`slow vital capacity in % (month of CSF sampling)` <- ifelse(clinical_info$`slow vital capacity in % (month of CSF sampling)` == "na",
                                                                           NA,
                                                                           clinical_info$`slow vital capacity in % (month of CSF sampling)`) 
clinical_info_interest <- clinical_info %>%
  filter(Group == "als") %>%
  dplyr::select(`CSF analyses ID (proteomic and metabolomic)`, # get only the relevant columns
                Gender,
                `Age at collection`,
                `pNFh (pg/ml)`,
                `Disease onset (location of first symptoms)`,
                `Age at onset`,
                `ALS progresion rate per month (delta ALS-FRSr / days *30)`,
                `slow vital capacity in % (month of CSF sampling)`) %>%
  dplyr::rename(Sample = `CSF analyses ID (proteomic and metabolomic)`, # alter the names as you wish to appear in the plot
                `age collection (y)` = `Age at collection`,
                `age onset (y)` = `Age at onset`,
                sex = Gender,
                `disease onset` = `Disease onset (location of first symptoms)`,
                `ALS progression rate` = `ALS progresion rate per month (delta ALS-FRSr / days *30)`,
                `slow vital capacity (%)` = `slow vital capacity in % (month of CSF sampling)`)

# get only the IDs available in the clinical info
input_proteomics <- input_proteomics[input_proteomics$Sample %in% clinical_info_interest$Sample,]
input_proteomics <- input_proteomics[match(clinical_info_interest$Sample,input_proteomics$Sample),] # put IDs in the same order
 
# -> phospho-metabolites

# -> smallRNA or miRNA (first DE analysis)


#### final tables #### 
input_proteomics_final <- merge(input_proteomics,clinical_info_interest)

##### Parse command line arguments ############
warnings    <- warnings();
args        <- commandArgs(TRUE);
inputfile   <- args[1];
outputfile  <- args[2];
rnaStartCol <- as.numeric(args[3]);
rnaEndCol   <- as.numeric(args[4]);  # from the end
colName     <- args[5]               # header name of the column which will be used as names of rows and columns
cmethod     <- args[6]
annCols     <- as.vector(unlist(strsplit(args[7], ','))) # "initals,age,gender,detailed_diagnosis,library_size,diagnosis"
annCols <- names(clinical_info_interest)[-1]
cohort      <- as.vector(unlist(strsplit(args[8], ','))) # "1 or 2 or 1,2 or 1,3"
diagnosis   <- as.vector(unlist(strsplit(args[9], ','))) # "0 or 1 or 0,1"
cmethod <- "p"
allData <- input_proteomics_final
rnaStartCol <- 1
rnaEndCol <- length(names(input_proteomics))

## Main logic  
  ## ############## DOWNLOAD AND LOAD THE DATA ###############
  
  ## Get the input data
  # cat("- Reading input file ...\n")
  # allDataOrig <- data.frame(read.table(inputfile,header=TRUE,sep="\t", na.strings="NA"))
  # 
  # # Subset the data for cohort and/or diagnosis
  # # allData <- allDataOrig[allDataOrig$cohort%in%cohort & allDataOrig$diagnosis%in%diagnosis,]
  # allData <- allDataOrig[allDataOrig$diagnosis%in%diagnosis,]
  
  # Get the number of colums
  ncols  <- length(names(allData))
  nrows  <- length(rownames(allData))
  # cat("- Total samples : ", nrows,"\n- Total features: ", ncols, "\n")
  
  # Get the smallncRNA data
  cat("- rnaStartCol: ", rnaStartCol, "\n- rnaEndCol: ", rnaEndCol, "\n\n")
  Xrna <- allData[,(rnaStartCol+1):rnaEndCol]
  
  # Get the rows and column headers
  column_header = names(Xrna) # ['hsa-let-7b-5p', 'hsa-let-7c-5p', 'hsa-let-7d-3p' ,...]
  # row_header    = allData[[colName]]      # [1, ..., 72, 73, 74, 78, 80, 83, 84, 85, 86, 90, 91,...]
  row_header = allData$Sample
  # Set index as colName (patient information)
  rownames(Xrna) <- allData$Sample
  
  # Get the annotation df needed to be ploted on the heatmap
  annDF <- subset(allData, select = names(allData) %in% annCols)
  annDF$`pNFh (pg/ml)` <- as.numeric(annDF$`pNFh (pg/ml)`)
  annDF$`ALS progression rate` <- as.numeric(annDF$`ALS progression rate`)
  annDF$`slow vital capacity (%)` <- as.numeric(annDF$`slow vital capacity (%)`)
  
  # Get the correlation matrix
  cat("- Computing pairwise correlation of columns, excluding NA/null values\n")
  if (cmethod == 'p'){
    corr_method <- 'pearson'
    cat ("\t- using pearson  : standard correlation coefficient\n")
  } else if (cmethod == 'k'){
    corr_method = 'kendall'
    cat("\t- using kendall  : Kendall Tau correlation coefficient\n")
  } else if (cmethod == 's'){
    corr_method = 'spearman'
    cat("\t- using spearman : Spearman rank correlation\n")
  }
  
  Xrna <- as.data.frame(sapply(Xrna, as.numeric)) # to convert count data to numeric
  
  # Calculate the correlation matrix
  D            <- as.data.frame(coop::cosine(t(Xrna)))
  rownames(D)  <- row_header
  colnames(D)  <- row_header
  itemLabels   <- row_header
  percentiles  <- FindOptimalBinning(D, itemLabels, transposeData=TRUE, verbose=TRUE)
  discreteData <- DiscretiseData(D, percentiles=percentiles)
  hc1          <- bhc(t(discreteData), itemLabels, verbose=TRUE)
  bhc_order    <- order.dendrogram(hc1)
  bhc_ordered_discreteData <- D[bhc_order,bhc_order]
  
  # Plot the heatmap with the new cluster and order and save it to the output file
  
  # Flatten matrix to scale and then reshape
  matbhc              <- as.matrix(bhc_ordered_discreteData)
  scaled_ordered_data <- matrix(scale(c(matbhc)), nrow=nrows, byrow=TRUE)
  scaled_df           <- as.data.frame(scaled_ordered_data)
  rownames(scaled_df) <- rownames(bhc_ordered_discreteData)
  colnames(scaled_df) <- rownames(bhc_ordered_discreteData)
 # colorP              <- colorRampPalette(c('dodgerblue3','deepskyblue2','lightblue','white','orange','red','darkred'))(256)
  colorP              <- circlize::colorRamp2(breaks = c(-4,-3,-2,0,2,3,4),
                                               colors = c('dodgerblue3',
                                                          'deepskyblue2',
                                                          'lightblue',
                                                          'white',
                                                          'orange',
                                                          'red',
                                                          'darkred'))
  colorAgeOnset       <- circlize::colorRamp2(c(20,100),
                                              c("#C6DBEF","#4292C6")) # blues
  colorAgeCollection  <- circlize::colorRamp2(c(40,100),
                                              c("#DADAEB","#807DBA")) # purples
  colorVitalCapacity  <- circlize::colorRamp2(c(0,150),
                                              c("#FDD0A2","#F16913")) # oranges
  colorALSRate        <- circlize::colorRamp2(c(0,50),
                                              c("#C7E9C0","#41AB5D")) # greens
  colorSex            <- c("Female" = "#57a8a2","Male" = "#A8575D") 
  colorpNFh           <- circlize::colorRamp2(c(0,15000),
                                              c("#FEE0D2","#FB6A4A")) # reds
  colorDiseaseOnset   <- c("spinal" = "#9EA9C4","bulbar" = "#7AEDCC","Both" = "#EDC06E")
  names(scaled_df) <- factor(names(scaled_df),levels = rownames(scaled_df))
  rownames(scaled_df) <- factor(rownames(scaled_df),levels = rownames(scaled_df))
  ha                  <- HeatmapAnnotation(df = as.data.frame(annDF),
                                           col = list(sex = colorSex,
                                                      `age onset (y)` = colorAgeOnset,
                                                      `age collection (y)` = colorAgeCollection,
                                                      `slow vital capacity (%)` = colorVitalCapacity,
                                                      `ALS progression rate` = colorALSRate,
                                                      `pNFh (pg/ml)` = colorpNFh,
                                                      `disease onset` = colorDiseaseOnset))
  res                 <- pheatmap(as.matrix(scaled_df) , cluster_rows=T, cluster_cols=T, 
                                  color = colorP, annotation_col=annDF,
                                  border_color = "grey50",treeheight_row = 0,
                                  main = "proteomics")
  res2                 <- ComplexHeatmap::Heatmap(as.matrix(scaled_df),
                                                  name = "correlation \nmagnitude",
                                                  col = colorP,show_column_dend = T,
                                                  show_row_dend = F,
                                                  top_annotation = ha)
  res2
  

############ USER DEFINED FUNCTIONS ##########
##Function to write out the items labels for each cluster.
WriteClusterLabels <- function(dendro, outputFile="", verbose=FALSE){
  ##----------------------------------------------------------------------
  ## DEFINE SOME FUNCTIONS TO USE RECURSIVELY ON THE DENDROGRAM NODES ----
  ##----------------------------------------------------------------------
  ##for ease, we use discrete height labels here
  ##this hardwires the logEvidence threshold at zero, which is the point
  ##where merge and not-merge hypotheses are equally likely  
  WhereToCut <- function(n){
    attr(n,"height") <- 1##default
    if (!is.leaf(n)){
      attr(n,"height") <- 2   
      if (attr(n, "logEvidence")<0)
        attr(n,"height") <- 3
    }
    n
  }
  ##----------------------------------------------------------------------
  ## PROCESS THE DENDROGRAM NODES RECURSIVELY ----------------------------
  ##----------------------------------------------------------------------
  dendro <- dendrapply(dendro, WhereToCut);
  str(dendro)
  
  ##----------------------------------------------------------------------
  ## CUT THE DENDROGRAM AND PRINT THE LABELS IN EACH CLUSTER -------------
  ##----------------------------------------------------------------------
  cutDendro     <- cut(dendro, 2)
  nClusters     <- length(cutDendro$lower)
  nTotalLabels  <- length(labels(dendro))
  outputStrings <- rep("", nTotalLabels+nClusters)
  counter       <- 1
  print(outputStrings)
  quit()  
  for (i in 1:nClusters) {
    print(cutDendro$lower[[i]])
    ##extract the current dendrogram
    currentCluster <- cutDendro$lower[[i]]
    currentLabels  <- labels(currentCluster) 
    nLabels        <- length(currentLabels) 
    ##for each cluster, construct and store the labels
    outputStrings[counter] <- paste("---CLUSTER", i, "---")
    counter                <- counter + 1
    for (j in 1:nLabels){
      outputStrings[counter] <- currentLabels[j]
      counter                <- counter + 1
    }
  }
  ##----------------------------------------------------------------------
  ## IF REQUIRED, WRITE OUT THE CLUSTER LABELS TO A FILE -----------------
  ##----------------------------------------------------------------------
  if (outputFile!="") write.table(outputStrings, file=outputFile, quote=FALSE, row.names=FALSE)
  ##----------------------------------------------------------------------
  ## IF REQUIRED, PRINT THE CLUSTER LABELS OUT TO SCREEN -----------------
  ##----------------------------------------------------------------------
  if (verbose) for (i in 1:length(outputStrings)) print(outputStrings[i], quote=FALSE)
}

## Call the main function in the end
main()