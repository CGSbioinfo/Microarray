# The script to call the data engine and the graph engine functions.
# Created by: scw65
# Created on: 26/04/2018

# Loading packages
cat("Loading required packages...\n")
suppressMessages(require(oligo))
suppressMessages(require(pd.hugene.2.1.st))
suppressMessages(require(pd.mogene.2.1.st))
suppressMessages(require(genefilter))
suppressMessages(require(limma))
suppressMessages(require(ggplot2))
suppressMessages(require(reshape))
suppressMessages(require(dendextend))
suppressMessages(require(colorspace))
suppressMessages(require(scatterplot3d))
suppressMessages(require(dplyr))
suppressMessages(require(gridExtra))
suppressMessages(require(gplots))
suppressMessages(require(plotly))

# Source the data and graph engine scripts to enable the functions to be run
cat("Sourcing the data processing and graph plotting functions...\n")
source("bin/data_engine.R")
source("bin/graph_engine.R")


# Get the analysis info file from the command line argument
analysis_info_file = commandArgs(TRUE)[1]
#analysis_info_file <- "/mnt/cgs-fs3/SERVICES/Bioinformatic/Steph_Wenlock/MicroarrayPipeline/NGS_Example_1000/Analysis/analysis_info.txt"
analysis_info <- read.table(analysis_info_file)
#print(analysis_info)

### Retrieve analysis directory and sample_info and comparision files 
analysis_directory <- analysis_info$V3[1]
array_info_folder <- as.character(analysis_info$V3[3])
sample_info <- paste0(analysis_directory,"/sample_info.csv")
comparisons <- read.csv(file=paste0(analysis_directory,"/comparisons.csv"), header=T)
target <- as.character(analysis_info$V3[6])


### Create the required output directory for the apt-cel-transformation cel files
if (analysis_info$V3[7] == TRUE){
  
  cels_dir <- paste0(analysis_directory,"/GCCN_cels")
  
  # Create a cels text file
  cels_list <- list.files(path=as.character(analysis_info$V3[2]), full.names=TRUE)
  cel_data_file <- paste0(analysis_directory,"/cels_list.txt")
  write.table(cels_list, file=cel_data_file, sep="\n", quote=F, row.names=F, col.names="cel_files")
  
  # Call transformation function to run GCCN-SST transformation 
  transform_data(cel_data_file, array_info_folder, cels_dir)
  
}else{
  cels_dir <- as.character(analysis_info$V3[2])
}

  
# Call function to load the cel file data using sample_info file
raw.data <- load_data(cels_dir, sample_info)


### Check whether any samples to remove
if (analysis_info$V3[13]== "None"){
   cat("No samples being excluded\n")
}else{
  samples_to_remove <- as.character(analysis_info$V3[13])
  raw.data <- exclude_samples(raw.data, samples_to_remove) # Run function to remove samples
}


### Run graph functions to produce raw data plots - if QC_raw == TRUE 
run_rawQC <- as.character(analysis_info$V3[7])

if (run_rawQC == TRUE){
  
  raw_data_folder <- paste0(analysis_directory,"/raw_data_plots")
  if (dir.exists(raw_data_folder)==FALSE){ # If the raw data plots directory does not exist, create it.
    dir.create(raw_data_folder, recursive=T)
  }

  PCA_plotting(raw.data, raw_data_folder, "raw") 
  MA_plotting(raw.data, raw_data_folder, "raw")
  density_plotting(raw.data, raw_data_folder, "raw")
  GC_content_plots(raw.data, raw_data_folder)
  clust_and_heatmap_plots(raw.data, raw_data_folder, "raw") 
  hierarch_cluster_plots(raw.data, raw_data_folder, "raw")
  expression_boxplot(raw.data, raw_data_folder, "raw")
}



### Run the control QC - if control_QC == TRUE
run_controlQC <- as.character(analysis_info$V3[14])

if (run_controlQC == TRUE) {
  # Get the QCC file
  QCC_file <- list.files(path=array_info_folder, pattern=".qcc", full.names=TRUE)
  QCC_info <- read.table(QCC_file, sep="\t", header=T)
  controlLevels <- get_controlLevels(raw.data, QCC_info, target) # Get the control levels by calling the function
  
  controls_output_dir <- paste0(analysis_directory,"/control_QC")
  controlQC_plots(raw.data, QCC_file, controls_output_dir)
}

### Run RMA normalization
rma.data <- RMAnorm(raw.data, target)


# Run graph functions to produce rma normalized data plots - if QC_norm == TRUE
run_normQC <- as.character(analysis_info$V3[8])

if (run_normQC == TRUE){
  
  rma_data_folder <- paste0(analysis_directory,"/rma_data_plots")
  if (dir.exists(rma_data_folder)==FALSE){ # If the raw data plots directory does not exist, create it.
    dir.create(rma_data_folder, recursive=T)
  }
  
  PCA_plotting(rma.data, rma_data_folder, "norm")
  MA_plotting(rma.data, rma_data_folder, "norm")
  density_plotting(rma.data, rma_data_folder, "norm")
  clust_and_heatmap_plots(rma.data, rma_data_folder, "norm")
  hierarch_cluster_plots(rma.data, rma_data_folder, "norm")
  expression_boxplot(rma.data, rma_data_folder, "norm")
}


# Run comparisons if group_comparisons == TRUE
paired_samples <- as.character(analysis_info$V3[4])
compareGroups <- as.character(analysis_info$V3[10])
results_Tab.output <- paste0(analysis_directory,"/results/tables")
results_Graph.output <-paste0(analysis_directory,"/results/graphs")

if (compareGroups == TRUE){

  summary_sig <- data.frame(Comparisons=character(), Significant_0.05=integer(), Significant_0.01 = integer(), n=integer(), stringsAsFactors=F)
  n_comp=0
  first_comp = TRUE
  
  for (i in 1:nrow(comparisons)){
    n_comp = n_comp+1
    groups <- c(as.character(comparisons[i,1]),as.character(comparisons[i,2]))
    results <- comparison_func(paired_samples, groups, rma.data, target)
    
    if (dir.exists(results_Tab.output)==FALSE){ # If the output directory does not exist, create it.
      dir.create(results_Tab.output,recursive=T)
    }
    if (dir.exists(results_Graph.output)==FALSE){ # If the output directory does not exist, create it.
      dir.create(results_Graph.output,recursive=T)
    } 
    
    cat("Writing results table\n")
    if (analysis_info$V3[7] == TRUE){
      write.table(results, file=paste0(results_Tab.output,"/Results_GCCN_SST_",as.character(comparisons[i,1]),"_VS_",as.character(comparisons[i,2]),".csv"), sep=",", row.names = F)
    }else{
      write.table(results, file=paste0(results_Tab.output,"/Results_non-GCCN_SST_",as.character(comparisons[i,1]),"_VS_",as.character(comparisons[i,2]),".csv"), sep=",", row.names = F)
    }
    
    # Call PCA function on each group
    PCA_plotting(rma.data, results_Graph.output, "norm", groups)
    volcano_plotting(results, results_Graph.output, groups)

    # Write significant summary table
    n_sig_0.05 <- length(subset(results, (results$adj.P.Val>0.05)))
    n_sig_0.01 <- dim(subset(results, (results$adj.P.Val<0.01)))
    comp<-paste0(groups[1],"_vs_",groups[2])
    #cat(paste0(n_comp," ",comp, " ",n_sig_0.05," ",n_sig_0.01,"\n"))

    if (first_comp){
      summary_sig[n_comp,]<-c(comp,n_sig_0.05,n_sig_0.01, nrow(results))
      first_comp=FALSE
    }else{
      summary_sig<-rbind(summary_sig, c(comp, n_sig_0.05, n_sig_0.01, nrow(results)))
    }
  }
  print(summary_sig)
  write.table(file=paste0(results_Tab.output, "/summary_comparisons.csv"),summary_sig, row.names = FALSE, sep=",", quote=FALSE)
}

  