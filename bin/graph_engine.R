# Script to hold all the functions to plot graphs
# Created by: scw65
# Created on: 20/04/2018

suppressMessages(require(oligo))
suppressMessages(require(pd.hugene.2.1.st))
suppressMessages(require(pd.mogene.2.1.st))
suppressMessages(require(genefilter))
suppressMessages(require(limma))
suppressMessages(require(ggplot2))
suppressMessages(require(reshape))
suppressMessages(require(dendextend))
suppressMessages(require(colorspace))
suppressMessages(require(dplyr))
suppressMessages(require(gridExtra))
suppressMessages(require(gplots))
suppressMessages(require(plotly))


#### PCA no shape functions (single and paired) ##########################################
#### PAIRED
pca_plot_pairs_noshape=function(x, plot_name, column_col, log=FALSE){

  extension <- strsplit(plot_name, split="\\.")     # Get the extension/plot device
  extension <- extension[[1]][-1]                   # Take always the last element which will be the device name

  cat("\n")
  cat("PCA plot pairs no shape\n")
  if (log==TRUE){
    pca=prcomp(t(log2(oligo::exprs(x)+1)))  # If the data needs to be transformed by log2 and transposed, e.g. raw data
  } else {
    pca=prcomp(t(oligo::exprs(x)))          # If the data is already normalised, don't need to re-transform it, just transpose
  }

  percentVar<-pca$sdev^2/sum(pca$sdev^2)
  percentVar=round(percentVar*100)
  names(percentVar)=paste0('PC',1:length(percentVar))
  pca_points<-pca$x
  pca_points<-data.frame(SampleName=row.names(pca_points),pData(x),pca_points)

  pcs=c('PC1','PC2')
  plot1=pca_plot_single_noshape(pca_points,pcs,column_col, percentVar=percentVar)
  pcs=c('PC1','PC3')
  plot2=pca_plot_single_noshape(pca_points,pcs,column_col, percentVar=percentVar)
  pcs=c('PC1','PC4')
  plot3=pca_plot_single_noshape(pca_points,pcs,column_col, percentVar=percentVar)
  pcs=c('PC2','PC3')
  plot4=pca_plot_single_noshape(pca_points,pcs,column_col, percentVar=percentVar)
  pcs=c('PC2','PC4')
  plot5=pca_plot_single_noshape(pca_points,pcs,column_col, percentVar=percentVar)
  pcs=c('PC3','PC4')
  plot6=pca_plot_single_noshape(pca_points,pcs,column_col, percentVar=percentVar)
  g <- grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,ncol=3)
  saved_plot <- ggsave(plot_name, plot=g, device=extension, width=21, height=10)

}

#### SINGLE
   pca_plot_single_noshape=function(x, pcs, col, percentVar){
  plus_margin=diff(range(x[,pcs[1]]))/6
  plot=ggplot(as.data.frame(x), aes_string(x=pcs[1], y=pcs[2], col=col)) +
    geom_point(size=3) +
    xlab(paste0(pcs[1],' ', percentVar[pcs[1]], "%")) + ylab(paste0(pcs[2], ' ', percentVar[pcs[2]], "%")) +
    ggtitle(paste0(pcs[1],' vs ',pcs[2])) +
    theme(plot.margin= unit(c(1, 1, 1, 1), "cm"),
          panel.border = element_rect(fill=NA),
          panel.grid.major =element_line(colour = 'grey', linetype='dashed'),
          legend.text = element_text(size=15),
          legend.title= element_text(size=15),
          plot.title = element_text(size=15, face='bold'),
          axis.text = element_text(size=16),
          axis.title = element_text(size=17))  +
    guides(colour = guide_legend(override.aes = list(size=2.2)))  +
    xlim(min(x[,pcs[1]])-plus_margin,max(x[,pcs[1]])+plus_margin)
  return(plot)
}
##########################################################################################
##########################################################################################




#### Function to plot the PCA graphs for specific or all groups ##########################
PCA_plotting=function(data, output_dir, type_of_data, groups=NA){  # By default groups == NA

    if (dir.exists(output_dir)==FALSE){                 # If the output directory does not exist, create it.
       dir.create(output_dir,recursive=T)
    }

    if (type_of_data == "raw"){                         # Means the data will be logged if raw and left if normalized
        log_choice = TRUE
    }else{
        log_choice = FALSE
    }

    if(length(groups)>1){
        group1=as.character(groups[1])  # Take each row numbers' 1st column
        group2=as.character(groups[2])  # Take each row numbers' 2nd column
        cat(paste0("Plotting PCA of: ",group1," vs ",group2,"\n"))

        # extract just the samples in the groups being compared
        index1<-row.names(pData(data)[pData(data)$Sample.Group%in%group1,])       # All the samples in the first group
        index2<-row.names(pData(data)[pData(data)$Sample.Group%in%group2,])       # All the samples in the second group
        index_all<-c(index1,index2)
        selected_samples<-(pData(data)[index_all,])
        selected_samples<-droplevels(selected_samples)

        group_data <- data[,(index_all)]
        pData(group_data)<- droplevels(pData(group_data))
        pData(group_data)$Sample.Group=relevel(pData(group_data)$Sample.Group, ref=group1) # reorder based on group1

        pca_plot_pairs_noshape(group_data, paste0(output_dir, "/pca_group_",group1,"_vs_",group2,"_",type_of_data,".png"), column_col = 'Sample.Group', log=log_choice)
    }else{                                  # Plot all groups if none are given
         pca_plot_pairs_noshape(data, paste0(output_dir, "/pca_allgroups_",type_of_data,".png"), column_col = 'Sample.Group', log=log_choice)
    }
}

######### END of function ################################################################
##########################################################################################




#### Function to create volcano plots ####################################################
volcano_plotting=function(results_table, output_dir, groups){ # Results table will come from the comparisons function in data_engine (will be for single group comparison)

  if (dir.exists(output_dir)==FALSE){                       # If the output directory does not exist, create it.
    dir.create(output_dir,recursive=T)
  }

  group1=as.character(groups[1])
  group2=as.character(groups[2])
  cat(paste0("Plotting volcano plot of: ",group1," vs ",group2,"\n"))

  logFC<-results_table$logFC
  geneNames <- results_table$geneassignment
  log_adj.pvalue=-log(results_table$adj.P.Val,base = 10) #-log10 to give +ve number
  significant<-results_table$adj.P.Val<0.05
  volcano<-as.matrix(cbind(logFC,log_adj.pvalue))

  # Normal volcano plot
  g<-ggplot(data.frame(volcano), aes(x=as.numeric(logFC),as.numeric(log_adj.pvalue),col=significant))+geom_point()+
    ggtitle(paste0("Volcano_",group1,"_vs_",group2))+ylab("significance level (-log(adj.pval))")

  # Interactive volcano plot
  p <- plot_ly(data=data.frame(volcano), x=logFC, y=log_adj.pvalue, type="scatter", text=geneNames,
              mode="markers", color=~significant, colors=c("red", "cadetblue"))%>%
              layout(title=paste0("Volcano_",group1,"_vs_",group2))

  plot_name=paste0(output_dir,"/Volcano_",group1,"_vs_",group2)

  oldir=getwd()
  setwd("~")

  # Create a folder to hold all the comparison ggplots - then move 'mv' whole folder
  temp_plots <- paste0(getwd(),"/interactive_plot_",group1,"_vs_",group2)
  if (dir.exists(temp_plots)==FALSE){
    dir.create(temp_plots,recursive=T)
  }

  saved_plot <- ggsave(paste0(plot_name,".png"), plot=g, device="png", width=15, height=10)
  htmlwidgets::saveWidget(as_widget(p), file=paste0(temp_plots,"/Volcano_",group1,"_vs_",group2,"_interactive.html"), selfcontained=FALSE)
  system(paste0("mv ",temp_plots," ",output_dir),ignore.stdout = TRUE, ignore.stderr = TRUE) # Ignore any messages from the system
  setwd(oldir)
}

######### END of function ################################################################
##########################################################################################



#### Function to create MA plots ##########################################################
MA_plotting=function(data, output_dir, type_of_data){

    all_groups <- levels(factor(pData(data)$Sample.Group))

    for(i in 1:length(all_groups)){
      data.s<-data[,as.character(pData(data)$Sample.Group)%in%all_groups[i]]
      cat("Plotting MA data for samples in: ", all_groups[i], "\n")

        if (ncol(data.s)>3){   # If there are more than 3 samples
        png(file=paste0(output_dir,"/MAplot_", type_of_data,"_",all_groups[i],".png",sep=""),width = 1000, height=1800)
        par(mfrow=c(5,3))
      } else {
        png(file=paste0(output_dir,"/MAplot_", type_of_data,"_",all_groups[i],".png",sep=""),width = 1000, height=400)
        par(mfrow=c(1,3))
      }
      for (j in 1:length(sampleNames(data.s))){
        oligo::MAplot(data.s, which=j, cex.main=1.5, cex.axis=1.5, cex.lab=1.5,main=paste0("vs pseudo-median ", all_groups[i]))
      }
      dev.off()
    }
}
######### END of function ################################################################
##########################################################################################



#### Function to create density plots ####################################################
density_plotting=function(data, output_dir, type_of_data){

    all_groups <- levels(factor(pData(data)$Sample.Group))

    for(i in 1:length(all_groups)){

      data.s<-data[,as.character(pData(data)$Sample.Group)%in%all_groups[i]] # N.B. this will only work if rows of pData are Sample.Name

      cat(paste0("Plotting density plots for samples in: ", all_groups[i], "\n"))
      png(file=paste0(output_dir,"/density_plot_",type_of_data,"_", all_groups[i],".png"), width=780)
      par(mar=c(5.1,4.1,4.1,12.1),xpd=T)

        if (type_of_data == "norm"){  # rma data is an expression set
          hist(data.s, col=1:length(pData(data.s)$Sample.Name), lwd=1.3, lty = c(1:length(pData(data.s)$Sample.Name)))
          legend(legend=pData(data.s)$Sample.Name, lty = c(1:length(pData(data.s)$Sample.Name)),lwd=1.3, col=1:length(pData(data.s)$Sample.Name), "topright",inset=c(-0.2,0))
          dev.off()
        }
        else{  # raw data (gene feature set) requires the which argument to be filled in
          hist(data.s, col=1:length(pData(data.s)$Sample.Name), which="all",lwd=1.3, lty = c(1:length(pData(data.s)$Sample.Name)))
          legend(legend=pData(data.s)$Sample.Name, lty = c(1:length(pData(data.s)$Sample.Name)),lwd=1.3, col=1:length(pData(data.s)$Sample.Name), "topright",inset=c(-0.2,0))
          dev.off()

        }
    }


}
######### END of function ################################################################
##########################################################################################


#### Function to create GC content plots #################################################
GC_content_plots=function(data, output_dir){

    pmSeq<-pmSequence(data)                                     # Access the sequence information of data
    pmsLog2<-log2(pm(data))                                     # Create an intensity matrix from data
    coefs<-getAffinitySplineCoefficients(pmsLog2, pmSeq)
    counts<-Biostrings::alphabetFrequency(pmSeq,baseOnly=T)     #Base only = true as DNA=TGCA, RNA=UGCA
    GCcontent<-ordered(counts[,"G"]+counts[,"C"])
    total_graph=ceiling(length(sampleNames(data))/4)            # A way of dividing up the graphs and knowing the max numbers
    graph_n=0

    cat("Plotting GC content plots\n")
    for(i in 1:length(sampleNames(data))){
      current_n=ceiling(i/4)
      if(current_n==graph_n)                                    # If currently at 0
      {
        plot(y=pmsLog2[,i], x=GCcontent, main=paste("GC content ",sampleNames(data)[i]), xlab="GC frequency", ylab="Intensity" )
      }else{
        if(current_n>1){                                        # So only 4 GC plots in 1 figure
          dev.off()
        }
        graph_n=current_n
        png(file=paste0(output_dir,"/GCContent_",graph_n,".png"),width = 480, height=480)
        par(mfrow=c(2,2))
        plot(y=pmsLog2[,i], x=GCcontent, main=paste("GC content ",sampleNames(data)[i]), xlab="GC frequency", ylab="Intensity" )
      }
    }
    graphics.off()
}

######### END of function ################################################################
##########################################################################################



#### Function to create clustering and correlation heatmap plots ####################################
clust_and_heatmap_plots=function(data, output_dir, type_of_data){

  cat("Plotting the correlation heat maps\n")
  cor_mat<-cor(oligo::exprs(data))
  png(file=paste0(output_dir,"/Correlation_heatmap_",type_of_data,"_groups.png"), height=15, width=20, units="cm", res=100)
  heatmap.2(cor_mat, trace="none", labRow=pData(data)$Sample.Group, cexRow=1.2, cexCol=1.2, margin=c(6,8))
  dev.off()

  png(file=paste0(output_dir,"/Correlation_heatmap_",type_of_data,"_samples.png"), height=15, width=20, units="cm", res=100)
  heatmap.2(cor_mat, trace="none", cexRow=1.2, cexCol=1.2, margin=c(6,8))
  dev.off()

}

######### END of function ################################################################
##########################################################################################


#### Function to create hierarchical clustering plots ####################################
hierarch_cluster_plots=function(data, output_dir, type_of_data){

  cat("Plotting the hierarchical clustering plots\n")
  dist_mat<-dist(t(oligo::exprs(data))) # transpose the distance matrix data (by default calculated on rows)
  png(file=paste0(output_dir,"/hierarchical_clust_",type_of_data,".png"), width=780)
  hierarch_clust <- plot(hclust(dist_mat), xlab="Expression distance matrix")
  dev.off()

}
######### END of function ################################################################
##########################################################################################



#### Function to create expression box plot ##############################################
expression_boxplot=function(data, output_dir, type_of_data){

  cat("Plotting the expression box plots\n")
  png(file=paste0(output_dir,"/expression_boxplot_",type_of_data,".png"), width=780)
  if (class(data)[1]=="GeneFeatureSet"){
    plot <- boxplot(data, which="all", main="Boxplot of the expression  values (log2 intensities) for all samples", las=2)
  }else{                   # If using a Clariom array, the class is an expression feature set
     plot <- boxplot(data, main="Boxplot of the expression  values (log2 intensities) for all samples", las=2)
  }
  dev.off()
}
######### END of function ################################################################
##########################################################################################






#### Function to create all the control QC  plots ########################################
controlQC_plots=function(data, QCC_info, output){                # QCC_info = file path of QCC file

    cat("Running all control QC plots\n")

     if (dir.exists(output)==FALSE){                             # If the controls output directory does not exist, create it.
        dir.create(output,recursive=T)
    }

    # Calls all control QC plots, starting with get_controlLevels, as its returned data is needed as input for other control plots
    controls_info <- read.table(QCC_info, sep="\t", header=T)                         # Read the QCC file

    controlsLevels <- get_controlLevels(data, controls_info)                          # 1st function call: Get the normalized control expression for all samples from diff probes
    bac_controls(controlsLevels, plot_name=paste0(output,'/control_bac_all.png'))     # Use the controlLevels data to extract bacterial control probes expression
    polya_controls(controlsLevels, plot_name=paste0(output,'/control_polyA_all.png')) # Use the controlLevels data to extract RNA control probes expression
    posneg_controls(controlsLevels, plot_name=paste0(output,'/control_positive_negative_all.png'))
    nuse_controls(raw.data, plot_name = paste0(output,'/control_NUSE_all.png'))
    rle_controls(raw.data, plot_name=paste0(output,'/control_RLE_all.png'))

    rank_table=get_rank_table(controlsLevels)
    write.table(rank_table, file=paste0(output,"/rank_table.txt"), sep="\t", quote=FALSE)

}

######### END of function ################################################################
##########################################################################################





#########################################################################################
###### DEFINING ALL THE CONTROL PLOT FUNCTIONS ##########################################
#########################################################################################

#### Function to get the control info ###################################################
get_controlLevels <- function(data, QCC_info, target="core"){  # Required for PolyA, signal-to-noise rank table, positive and negative controls and BAC spike-in controls

  # Check the class of data before applying the rma (need 'target' arg for gene_st and not for clariom)
  if (class(data)[1]=="GeneFeatureSet"){
    rma.probe<-rma(data, target=target)
    match_control_probe<-match(QCC_info$probeset_id, row.names(rma.probe)) # Match control file data with the raw.data (match probe ids)
    controlsLevels<-oligo::exprs(rma.probe)[match_control_probe,]                                               # Take expression levels of only probes in the control file
    controlsLevels<-cbind(as.data.frame(controlsLevels), QCC_info$probeset_name, QCC_info$group_name)    # Combine expression of the samples, probe types/names and group name
    colnames(controlsLevels)[length(colnames(controlsLevels))]<-"probe_group"                             # Name the last column probe group
    colnames(controlsLevels)[length(colnames(controlsLevels))-1]<-"probe_type"                           # The second last column (all other colnames are sample names)

  }
  else{ # If using a Clariom array, the class is an expression feature set
    rma.probe<-rma(data)
    match_control_probe<-match(as.character(QCC_info$probeset_name), as.character(row.names(rma.probe))) # Match control file data with the raw.data (match probe ids)
    controlsLevels<-oligo::exprs(rma.probe)[match_control_probe,]                                               # Take expression levels of only probes in the control file
    controlsLevels<-cbind(as.data.frame(controlsLevels), QCC_info$probeset_name, QCC_info$group_name)    # Combine expression of the samples, probe types/names and group name
    colnames(controlsLevels)[length(colnames(controlsLevels))]<-"probe_group"                            # Name the last column probe group
    colnames(controlsLevels)[length(colnames(controlsLevels))-1]<-"probe_type"
  }


  ## DEBUG
  #write.table(controlsLevels, file="Control_check.csv", sep=",")
  return(controlsLevels)  # Returns the expression levels of the control samples.
}

######### END of function ################################################################
##########################################################################################





#### Function to plot the BAC controls ###################################################
bac_controls<-function(controlsLevels, plot_name){

  extension <- strsplit(plot_name, split="\\.")     # Get the extension/plot device
  extension <- extension[[1]][-1]                   # Take always the last element which will be the device name

  cat("Running and plotting BAC spike-in controls\n")
  BioControls<-controlsLevels[grep("[B|b]io.+5_at", controlsLevels$probe_type),]        # Get BioB,C,D spike controls (regex: . any char, + any length)
  CreControls<-controlsLevels[grep("cre-5_at",controlsLevels$probe_type),]              # Get the cre spike controls
  BioControls<-rbind(BioControls,CreControls)                                           # Combine into a table
  BioControls<-melt(BioControls[,-length(colnames(BioControls))], id="probe_type")      # Melt so 1 ID variable belongs to probe_type all others are melted together
  g<-ggplot(data=BioControls, aes(colour=probe_type,y=value, group=probe_type, x=as.character(variable)))+
    geom_point()+geom_line()
  g<-g+theme(axis.text.x=element_text(angle=90, vjust=.5, size=12), panel.grid.major =element_line(colour = 'grey', linetype='dashed'))+labs(x="Samples", y="Intensity", colour="Control")
  saved_plot <- ggsave(plot_name, plot=g, device=extension, width=15, height=10)
}

######### END of function ################################################################
##########################################################################################




#### Function to plot the polyA controls #################################################
polya_controls<-function(controlsLevels, plot_name){

  cat("Running and plotting Poly A controls\n")

  # Get the extension/plot device
  extension <- strsplit(plot_name, split="\\.")     # Get the extension/plot device
  extension <- extension[[1]][-1]                   # Take always the last element which will be the device name

  polyAControls<-controlsLevels[grep("polya", controlsLevels$probe_group),]          # Select the probe group (polyA_spike control probes)
  polyAControls<-polyAControls[grep("r2.+5.+st", polyAControls$probe_type),]         # Select from those, the r2...st probe type
  polyAControls<-melt(polyAControls[,-length(colnames(polyAControls))], id="probe_type")
  g<-ggplot(data=polyAControls, aes(colour=probe_type,y=value, group=probe_type, x=as.character(variable)))+
    geom_point()+geom_line()
  g<-g+theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), panel.grid.major =element_line(colour = 'grey', linetype='dashed'))+labs(x="Samples", y="Intensity", colour="Control")
  saved_plot <- ggsave(plot_name, plot=g, device=extension, width=15, height=10)
}

######### END of function ################################################################
##########################################################################################




#### Function to plot the positive and negative controls #################################
posneg_controls<-function(controlsLevels, plot_name){

  extension <- strsplit(plot_name, split="\\.")     # Get the extension/plot device
  extension <- extension[[1]][-1]                   # Take always the last element which will be the device name

  cat("Running and plotting positive and negative controls\n")
  posControls<-controlsLevels[grep("pos_control", controlsLevels$probe_group),]
  negControls<-controlsLevels[grep("neg_control", controlsLevels$probe_group),]
  pos_cont_avg<-colMeans(posControls[, -c(length(colnames(posControls)),length(colnames(posControls))-1)]) #Last col lists the probe groups, ignore it
  neg_cont_avg<-colMeans(negControls[, -c(length(colnames(negControls)),length(colnames(negControls))-1)])
  neg_cont_avg<-cbind(as.data.frame(neg_cont_avg), "neg_control")
  pos_cont_avg<-cbind(as.data.frame(pos_cont_avg), "pos_control")
  neg_cont_avg<-cbind(as.data.frame(neg_cont_avg),row.names(neg_cont_avg))
  pos_cont_avg<-cbind(as.data.frame(pos_cont_avg), row.names(pos_cont_avg))
  colnames(neg_cont_avg)<-c("Intensity","Control","SampleID")
  colnames(pos_cont_avg)<-c("Intensity","Control", "SampleID")
  posneg_contr<-rbind(pos_cont_avg,neg_cont_avg)
  g<-ggplot(data=posneg_contr, aes(colour=Control, y=Intensity, x=SampleID,group=Control))+
    geom_point()+geom_line()
  g<-g+theme(axis.text.x=element_text(angle=90, size=12, vjust=0.5), panel.grid.major =element_line(colour = 'grey', linetype='dashed'))+labs(x="Samples", y="Intensity", colour="Control")
  saved_plot <- ggsave(plot_name, plot=g, device=extension, width=15, height=10)
}

######### END of function ################################################################
##########################################################################################




#### Function to produce the signal noise rank table #####################################
get_rank_table<-function(controlsLevels){

  cat("Running and writing the signal to noise rank table\n")
  posControls<-controlsLevels[grep("pos_control", controlsLevels$probe_group),]
  negControls<-controlsLevels[grep("neg_control", controlsLevels$probe_group),]
  pos_cont_avg<-colMeans(posControls[, -c(length(colnames(posControls)),length(colnames(posControls))-1)])
  neg_cont_avg<-colMeans(negControls[, -c(length(colnames(negControls)),length(colnames(negControls))-1)])
  signal_to_noise<-pos_cont_avg-neg_cont_avg
  signal_to_noise_rank<-rank(-signal_to_noise)
  rank_table<-as.data.frame(signal_to_noise)
  rank_table<-cbind(rank_table, signal_to_noise_rank)
  colnames(rank_table)[1]<-"Signal to noise"
  colnames(rank_table)[2]<-"Signal to noise rank (higher better)"
  return(rank_table)
}

######### END of function ################################################################
##########################################################################################




#### Function to plot the NUSE controls ##################################################
nuse_controls<-function(data, plot_name){

  cat("Running and plotting NUSE controls\n")
  pmData<-fitProbeLevelModel(data)
  png(file=plot_name, width=1080)
  par(mar=c(8.1,4.1,4.1,2.1))
  NUSE(pmData, las=2,cex.axis=1.25, cex.lab=1.24)
  dev.off()
}

######### END of function ################################################################
##########################################################################################




#### Function to plot the RLE controls ###################################################
rle_controls<-function(data,plot_name){

  cat("Running and plotting RLE controls\n")
  pmData<-fitProbeLevelModel(data)
  png(file=plot_name, width=1080)
  par(mar=c(8.1,4.1,4.1,2.1))
  RLE(pmData,las=2,cex.axis=1.25, cex.lab=1.25)
  dev.off()
}

######### END of function ################################################################
##########################################################################################








