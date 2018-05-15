# Created by: scw65
# Created on: 09/05/2018


#### Function to hold the non-changing chunks of report #################################
report_template=function(array_type, title, groups, comparions){
  
  RMD_header <- paste0(
"---
title: ",title," analysis report
output: 
  html_document:
    number_section: yes
    toc: yes
    fig_caption: true
bibliography: biblio.bib
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = getwd())
```
  ")
  
  
  Materials_Methods <- readLines("bin/report_temps/Materials_Methods.txt", warn=F)
  Materials_Methods.s <- strsplit(Materials_Methods, "//")
  Materials_Methods.s.array <- paste0(Materials_Methods.s[[1]][1],"\n",Materials_Methods.s[[2]][1],
                                      "\n",Materials_Methods.s[[3]][1],array_type,Materials_Methods.s[[3]][2])
  
  Bioinfo <- readLines("bin/report_temps/Bioinfo.txt", warn=F)
  Normalization <- readLines("bin/report_temps/Normalisation.txt", warn=F)
  TechnicalQC <- readLines("bin/report_temps/TechnicalQC.txt", warn=F)
  
  # Need to call the functions to retrieve the appropriate plots
  NUSE <- readLines("bin/report_temps/NUSE.txt", warn=F)
  NUSE_plot <- report_custom("NUSE", "all_groups", NA)
  
  RLE <- readLines("bin/report_temps/RLE.txt", warn=F)
  RLE_plot <- report_custom("RLE", "all_groups", NA)
  
  BAC <- readLines("bin/report_temps/BAC.txt", warn=F)
  BAC_plot <- report_custom("bac", "all_groups", NA)
  
  PolyA <- readLines("bin/report_temps/PolyA.txt", warn=F)
  PolyA_plot <- report_custom("polyA", "all_groups", NA)
  
  PosNeg <- readLines("bin/report_temps/PosNeg.txt", warn=F)
  PosNeg_plot <- report_custom("positive_negative", "all_groups", NA)
  
  Density_raw <- readLines("bin/report_temps/Density_raw.txt", warn=F)
  Density_raw_plot <- suppressWarnings(report_custom("density_plot_raw_", groups, "raw"))
  
  Density_RMA <- readLines("bin/report_temps/Density_rma.txt", warn=F)
  Density_RMA_plot <- suppressWarnings(report_custom("density_plot_norm_",  groups, "rma"))
  
  MA_raw <- readLines("bin/report_temps/MA_raw.txt", warn=F)
  MA_raw_plot <- suppressWarnings(report_custom("MAplot_raw_", groups, "raw"))
  
  MA_RMA <- readLines("bin/report_temps/MA_rma.txt", warn=F)
  MA_RMA_plot <- suppressWarnings(report_custom("MAplot_norm_", groups, "rma"))
  
  GCcont <- readLines("bin/report_temps/GC_cont.txt", warn=F)
  NoGCplots <- list.files(path="Analysis/raw_data_plots", pattern="GCContent_.*")
  NoGCplots <- c(1:length(NoGCplots))             # Get the number of GC content plots
  GCcontent_plots <- suppressWarnings(report_custom("GCContent_", NoGCplots, "raw"))
  
  Exprs_raw <- readLines("bin/report_temps/ExprsBoxplot_raw.txt", warn=F)
  Exprs_raw_plot <- report_custom("expression_boxplot_raw",  "all_groups", "raw")
  
  Exprs_RMA <- readLines("bin/report_temps/ExprsBoxplot_rma.txt", warn=F)
  Exprs_RMA_plot <- report_custom("expression_boxplot_norm",  "all_groups", "rma")
  
  PCA_raw <- readLines("bin/report_temps/PCA_raw.txt", warn=F)
  PCA_raw_plot <- report_custom("pca_allgroups_raw", "all_groups", "raw")
  
  PCA_RMA <- readLines("bin/report_temps/PCA_rma.txt", warn=F)
  PCA_RMA_plot <- report_custom("pca_allgroups_norm", "all_groups", "rma")
  
  HClust_raw <- readLines("bin/report_temps/HClust_raw.txt", warn=F)
  HClust_raw_plot <- report_custom("hierarchical_clust_raw", "all_groups", "raw")
  
  HClust_RMA <- readLines("bin/report_temps/HClust_RMA.txt", warn=F)
  HClust_RMA_plot <- report_custom("hierarchical_clust_norm",  "all_groups", "rma")
  
  CorrHeatMap_raw <- readLines("bin/report_temps/CorrHeatMap_1.txt", warn=F)
  CorrHeatMap_raw_plots <- suppressWarnings(report_custom("Correlation_heatmap_raw_",  c("groups","samples"), "raw"))
  
  CorrHeatMap_rma <- readLines("bin/report_temps/CorrHeatMap_2.txt", warn=F)
  CorrHeatMap_rma_plots <- suppressWarnings(report_custom("Correlation_heatmap_norm_",  c("groups","samples"), "rma"))
  
  Comparisons_1 <- readLines("bin/report_temps/Comparisons_1.txt", warn=F)
  
  Individual_comps <- individualComps(comparisons)
  
  combined_templates <- c(RMD_header,"\n", Materials_Methods.s.array,"\n",Bioinfo,
                          "\n",Normalization,"\n",TechnicalQC,"\n",NUSE,"\n",NUSE_plot,
                          "\n",RLE,"\n",RLE_plot,"\n",BAC,"\n",BAC_plot,"\n", PolyA,
                          PolyA_plot,"\n",PosNeg,"\n",PosNeg_plot,"\n",Density_raw,
                          "\n", Density_raw_plot,"\n",Density_RMA,"\n",Density_RMA_plot,
                          "\n",MA_raw,"\n",MA_raw_plot,"\n",MA_RMA,"\n",MA_RMA_plot,"\n",GCcont,
                          "\n", GCcontent_plots,"\n",Exprs_raw,Exprs_raw_plot,"\n","\n",Exprs_RMA,"\n",Exprs_RMA_plot,
                          "\n",PCA_raw,"\n",PCA_raw_plot,"\n", PCA_RMA,"\n", PCA_RMA_plot,
                          "\n",HClust_raw,"\n", HClust_raw_plot,"\n",HClust_RMA,
                          "\n",HClust_RMA_plot,"\n",CorrHeatMap_raw,"\n",CorrHeatMap_raw_plots,
                          "\n",CorrHeatMap_rma,"\n",CorrHeatMap_rma_plots,"\n",
                          Comparisons_1,"\n", Individual_comps, "\n# References")
  
  return(combined_templates)
  
}

######### END of function ################################################################
##########################################################################################




#### Function to retrieve the plots/tables for each individual comparison group ###########
individualComps=function(comparisons_tab){
  
  combinedComps <- NULL
  for (i in 1:nrow(comparisons_tab)){
    groups <- c(as.character(comparisons[i,1]),as.character(comparisons[i,2]))
    header <- paste0("\n### Group ",groups[1]," v Group ", groups[2])
    
    code <- paste0("```{r echo=FALSE, results='asis'}
    library(knitr)
    results<-read.csv('results/tables/Results_GCCN_SST_",groups[1],"_VS_",groups[2],".csv')
    kable(results[1:10,1:8], caption='",groups[1],"v",groups[2],"')
    ```

    
[Link to the table](results/tables/Results_GCCN_SST_",groups[1],"_VS_",groups[2],".csv)")

   plots <- paste0("![Volcano plot](results/graphs/Volcano_",groups[1],"_vs_",groups[2],".png)
                   [Link to interactive volcano plot](results/graphs/interactive_plot_",groups[1],"_vs_",groups[2],"/Volcano_",groups[1],"_vs_",groups[2],"_interactive.html)
                   ![PCA plot](results/graphs/pca_group_",groups[1],"_vs_",groups[2],"_norm.png)")
    
    combinedComps<-c(combinedComps, header, code, plots)
  }
  return(combinedComps)
}


######### END of function ################################################################
##########################################################################################




#### Function to retrieve the correct plots for each section of the report ###############
report_custom=function(section, groups, dtype){
  
  if (groups == "all_groups"){
    if (is.na(dtype)){
      section_response = paste0("![",section," control](control_QC/control_",section,"_all.png)")
    }else
      section_response = paste0("![",section,"](",dtype,"_data_plots/",section,".png)")
    
  }
  else{
    section_response <- NULL
    
    for (i in 1:length(groups)){
      if (section == "GCContent_"){
        header <- paste0("\n### ",section,"GC content plot ",groups[i],"\n")
      }else if (grepl("Correlation", section) ==TRUE) {
        header <- paste0("\n####",section,"_",groups[i],"\n")
      }else{
        header <- paste0("\n#### ",section,"group_",groups[i],"\n")
      }
      
      if (dtype == "raw"){
        plot <- paste0("![", section," ",groups[i],"](raw_data_plots/",section,groups[i],".png)")
      }else{
        plot <- paste0("![", section," ",groups[i],"](rma_data_plots/",section,groups[i],".png)")
      }
      section_response <- c(section_response, header, plot)
    }
  }
    
  return(section_response)
}



######### END of function ################################################################
##########################################################################################



#### Function to generate all report files ###############################################
report_generation=function(report_name, array_type, project_name, project_folder, groups, comparisons){
  
  setwd(as.character(project_folder))
  
  # Create the analysis report script
  if (file.exists(report_name)==FALSE){ # If the output directory does not exist, create it.
    file.create(report_name,recursive=T)
  }
  
  combined_report <- report_template(array_type, project_name, groups, comparisons)
  write.table(combined_report, report_name, quote=F, row.names = F, col.names = F)
  
  # Compile the analysis script into a report.
  rmarkdown::render(report_name)
  
}

######### END of function ################################################################
##########################################################################################



