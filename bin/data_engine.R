# Script to hold all the functions for the data handling
# Created by: scw65
# Created on: 19/04/2018

suppressMessages(require(oligo))
suppressMessages(require(pd.hugene.2.1.st))
suppressMessages(require(pd.mogene.2.1.st))
suppressMessages(require(genefilter))


#### Function to run apt-cel-transformer  #################################################
transform_data=function(cel_data, array_info_folder, output){

    if (dir.exists(output)==FALSE){ # If the output directory does not exist, create it.
        dir.create(output,recursive=T)
    }

    QCC <- list.files(path=array_info_folder, pattern=".qcc", full.names=TRUE)
    PGF <- list.files(path=array_info_folder, pattern=".pgf", full.names=TRUE)
    CLF <- list.files(path=array_info_folder, pattern=".clf", full.names=TRUE)

    # Print error message if any required files are not found
    error_raised = FALSE
    if (length(QCC)==0){
        cat("Error: QCC file not found.")
        error_raised = TRUE
    }
    if (length(PGF)==0){
        cat("Error: PGF file not found.")
        error_raised = TRUE
    }
    if (length(CLF)==0){
        cat("Error: CLF file not found.")
        error_raised = TRUE
    }

    if(!error_raised){
        cat("Running GCCN_SST transformation on the raw cel files.\n")
        command_txt <- paste0("apt-cel-transformer --pgf-file ",PGF," --clf-file ",CLF,
                                " --qc-probesets ",QCC," --cel-files ",cel_data,
                                " -c gc-correction,scale-intensities -o ",output)

        system(command=command_txt)  # Run the apt-cel-transformer (to apply the GCCN_SST correction)
        }
    }

######### END of function ################################################################
##########################################################################################



#### Function to load the raw cel files (transformed or not) #############################
load_data=function(cel_folder, sample_info){

  info <- read.csv(sample_info, header=T)

  cat("Retrieving raw data\n")
  raw.celfiles <- list.celfiles(cel_folder,full.names=T)
  raw.data<-read.celfiles(raw.celfiles)

  # Process sample names by removing the chip identifier
  chip_pattern <- "_\\(.+\\)_.*"  #.+ means any char - should be consistent for all chip names in CEL files

  sampleNames(raw.data)<- gsub(sampleNames(raw.data), pattern = chip_pattern, replacement = "")
  sampleNames(raw.data)<- gsub(sampleNames(raw.data), pattern = " ",replacement = "_")
  cat("Samples loaded:\n")
  cat(sampleNames(raw.data),"\n")

  # Check the IDs match between the sample info file and pData
  index_sample<-match(as.character(info$Sample.Name),sampleNames(raw.data))
  row.names(info)<-info$Sample.Name

  error_raised = FALSE
  if (NA %in% index_sample){
    cat("Sample names do not match\n")
      error_raised = TRUE
  }

  # Select samples present in the info file
  if(!error_raised){
      raw.data<-raw.data[,index_sample]
      pData(raw.data)=info
      cat("Sample info:\n")
      print(pData(raw.data))

      # Return the raw data
      return(raw.data)

  }
}

######### END of function ################################################################
##########################################################################################




#### Function to exclude samples #########################################################
exclude_samples=function(data, samples_to_remove){

  cat(paste0("Removing: ",samples_to_remove,"\n"))
  new_data <- data[, !(sampleNames(data) %in% samples_to_remove)]

  return(new_data)
}

######### END of function ################################################################
##########################################################################################




#### Function to apply RMA normalization #################################################
RMAnorm <- function(data, target="core"){

  cat("Applying RMA normalization\n")
  # Check class of the data
  if (class(data)[1]=="GeneFeatureSet"){
    rma.data<-rma(data, target=target)
  }else{                   # If using a Clariom array, the class is an expression feature set
    rma.data<-rma(data)
  }

  return(rma.data)
}

######### END of function ################################################################
##########################################################################################




#### Function to compare groups of data #################################################
comparison_func=function(paired_design, comparison_groups, data, gene_name_conv=TRUE, target="core"){

  # Set up the annotation target
  if (target == "core"){
    annot <- "transcript"
  }
  else{
    annot <- target
  }

  # Split the 2 groups and begin comparison
  group1 <- as.character(comparison_groups[1])
  group2 <- as.character(comparison_groups[2])
  cat(paste0("Comparing groups: ",group1," v ",group2,"\n"))

  index1<-row.names(pData(data)[pData(data)$Sample.Group%in%group1,])       # All the samples in the first group
  index2<-row.names(pData(data)[pData(data)$Sample.Group%in%group2,])       # All the samples in the second group
  index_all<-c(index1,index2)
  selected_samples<-(pData(data)[index_all,])
  selected_samples<-droplevels(selected_samples)
  # print(selected_samples)

  data.s <- data[,(index_all)]
  pData(data.s)<- droplevels(pData(data.s))
  pData(data.s)$Sample.Group=relevel(pData(data.s)$Sample.Group, ref=group1) # reorder based on group1


  if (class(data.s)=="GeneFeatureSet"){
    cat("Normalizing raw data 'GeneFeatureSet', using RMA, to retrieve annotation\n")
    data.s <- rma(data.s,target=target)
  }

  ### Changed to 25% quantile range - need to be here as rma data will be logged and raw will not
  median_exprs<-quantile(oligo::exprs(data.s),probs=0.5)

  # Filter data
  cat("Filtering data...\n")
  filter_setting<-kOverA(length(sampleNames(data.s))/2,median_exprs)
  filter_function<-filterfun(filter_setting)
  probe_filter<-genefilter(oligo::exprs(data.s),filter_function)
  data.s.f<-data.s[probe_filter,]
  #cat("Row name after probe filtering\n")
  #print(head(row.names(exprs(data.s.f))))

  # Get annotation
  cat("Retrieving annotation at the ",annot," level...\n")
  data.annot<-getNetAffx(data.s.f, type=annot)
  fData(data.s.f)<-pData(data.annot)

  # Filtering according to annotation
  data.s.f.a<-data.s.f[!rowSums(is.na(fData(data.s.f)[,c(8:11)]))==4,]
  #cat("Row names after NA filtering\n")
  #print(head(row.names(exprs(data.s.f.a))))

  cat(paste0("Fitting the linear model and performing differential testing between the groups: ",group1," v ",group2,"\n"))

  if (paired_design == TRUE){  # Include the SibShip data in the design matrix
    design <- model.matrix(~pData(data.s.f.a)$SibShip+pData(data.s.f.a)$Sample.Group)
  }else{
    design <- model.matrix(~pData(data.s.f.a)$Sample.Group)
  }

  colnames(design)<-gsub(colnames(design), pattern=".*Sample.Group", replacement="")
  fit <-lmFit(data.s.f.a, design)
  ##contrast.matrix=makeContrasts(paste(as.character(group2),as.character(group1),sep="-"),levels=design)
  ##fit.contrast<-contrasts.fit(fit, contrast.matrix)
  fit.ebayes<-eBayes(fit)   # Don't always need a contrast matrix for more simple experiments - can jump straight to here

  cat("Creating results table\n")
  results.table <-topTable(fit.ebayes, number=dim(oligo::exprs(data.s.f.a))[1], adjust.method="BH", coef = colnames(design)[length(colnames(design))]) # coef takes the last column in the design matrix (sample group)
  results.ids <- row.names(results.table)
  results<-cbind(results.table, oligo::exprs(data.s.f.a)[results.ids,])

  # Convert the gene assignment column to the given gene symbol if conversion is TRUE
  if (gene_name_conv == TRUE){
    ResultsTable <- GeneNameChange(results)
  }

  else{
    ResultsTable <- results
  }

  return(ResultsTable)
}

######### END of function ################################################################
##########################################################################################

######### END of function ################################################################
##########################################################################################



#### Function to convert results table gene assignment column to gene symbol #############
GeneNameChange <- function(x_table){
  x_table$geneassignment<-as.character(x_table$geneassignment)     # Change table column from factor to char so can be changed
  for (i in 1:nrow(x_table)){
    currentGene <- (as.character(x_table$geneassignment[i]))       # The row corresponding to each gene in the table
    splitNames <- strsplit(currentGene, "//")                      # Need character vector for strsplit, split on "//"
    GeneName <- splitNames[1]                                      # Hold the values of each split gene name
    if(is.na(GeneName)){                                           # If the GeneName is NA
      x_table$geneassignment[i]<-"NA"                              # Reassign it as 'NA'
    }  else{
      ExtractedName <- GeneName[[1]][2]
      x_table$geneassignment[i]<-ExtractedName                     # Replace the geneassignment column with the ExtractedName
    }
  }

  return(x_table)                                                  # Return the converted table at the end of the function
}

######### END of function ################################################################
##########################################################################################




