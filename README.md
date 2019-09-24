Microarray (Affymetrix) pipeline
================

-   [Getting Started](#getting-started)
-   [Running the pipeline](#running-the-pipeline)
    -   [Analysis info file](#analysis-info-file)
        -   [Format of the analysis info file](#format-of-the-analysis-info-file)
        -   [Filling in the analysis info file](#filling-in-the-analysis-info-file)
        -   [About the array files](#array-files)
    -   [Sample information and pairwise comparisons](#sample-information-and-group-comparisons)
-   [Description of the scripts used in the pipeline](#description-of-the-r-scripts)
-   [Standard pipeline command example](#standard-pipeline-command-example)
-   [Standard output directories and files](#standard-output-files-running-the-above-command-line)
    
    
## Getting Started
1. Create a new folder in cgs-fs3 or any other location with your chosen project name - referred to as the "project folder" in this README.   
2. Create a "bin/" folder in the project folder.   
3. Download and copy the scripts from the github "Microarray/bin" to the "bin/" folder. 
4. Create an "Analysis" folder in the project folder. 
5. In this "Analysis" folder, create two more folders, "Array_files" and "cels"
6. Copy the corresponding QCC, PGF, and CLF files for your array and species into the "Array_files" folder. 
7. Copy the **raw** CEL files for each sample into the "cels" folder. 
8. If you want to automatically generate a report for your analysis, you will also need to download the 'biblio.bib' file and place it inside the 'Analysis' folder, so the report references can be retrieved. 


## Running the pipeline

### Analysis info file
A central part of the pipeline is the **analysis info** file. It has information about the analysis directory location, the original CEL files folder, the array files folder and additional parameters used throughout the analysis. 

#### Format of the analysis info file
The analysis info file is a simple .txt file, with each line providing information. Different parameters on the same line are separated by the semicolon (i.e ";") character. These arguments are in a set order so please do not change the order of the lines. A template analysis_info.txt file can be downloaded in the github Microarray folder.    


The following is an example of the analysis info file:   

----------------------------------------------------------------------------------------------------
**analysis_directory =** /mnt/cgs-fs3/SERVICES/Microarrays/GX-Affy/GX-A.Tester-10235/Analysis

**rawCELs_folder =** /mnt/cgs-fs3/SERVICES/Microarrays/GX-Affy/GX-A.Tester-10235/Analysis/cels 

**Array_info =** /mnt/cgs-fs3/SERVICES/Microarrays/GX-Affy/GX-A.Tester-10235/Analysis/Array_files

**paired_samples =** FALSE 

**GCCN-SST_correction =** TRUE

**data_target =** core

**QC_raw =** TRUE

**QC_norm =** TRUE

**GEO_output =** TRUE

**group_comparisons =** TRUE

**report_gen =** TRUE

**functional_analysis =** TRUE

**samples_to_exclude =** None

**control_QC =** TRUE

**Batch_correction =** FALSE

------------------------------------------------------------------------------------------------------
<br>

The following is an explanation of the lines within the analysis info file:  

------------------------------------------------------------------------------------------------------
**analysis_directory =** *\<path to directory of the analysis\>* 

**rawCELs_folder =** *\<path to the raw CEL files folder\>*  

**Array_info =** *\<path to the array folder, containing the QCC, PGF and CLF files\>*   

**paired_samples=** *TRUE or FALSE whether the samples are paired and a SibShip column is required in the "comparisons.csv" and design matrix.*

**GCCN-SST_correction =** *TRUE or FALSE to perform GCCN-SST normalization on the raw CEL files. Default = TRUE*

**data_target =** *The level at which the RMA normalization is performed. Options are "core" or "probeset", which refer to the transcript and probeset level data, respectively. Default = core.*

**QC_raw =** *TRUE or FALSE to perform QC on the raw data. Default = TRUE.*

**QC_norm =** *TRUE or FALSE to perform QC on the RMA normalized data. Default = TRUE.*

**GEO_output =** *TRUE or FALSE to output the data required to produce a GEO submission. Default = TRUE.*

**group_comparisons =** *TRUE or FALSE to perform group comparisons on the data, using the comparisons from the "comparisons.csv" file. Default = TRUE.*

**report_gen =** *TRUE or FALSE to generate a report on the whole analysis, once it has been completed. Default = TRUE.*

**functional_analysis =** *TRUE or FALSE to perform functional analysis on the data. Default = TRUE.*

**samples_to_exclude =** *The name of any samples which need to be excluded. Default = None. Multiple sample names should be separated with "," like 'Sample1,Sample2,Sample3'.*

**control_QC =** *TRUE or FALSE to perform QC on the control data, using the corresponding QCC array file. Default = TRUE.*

**Batch_correction =** *TRUE or FALSE to perform a batch correction (using ComBat: sva) on the data. Default = FALSE. To do this you would need to add the batch information in a 'Batch' column in the sample_info file.*
 
 ------------------------------------------------------------------------------------------------------  
   
### Filling in the analysis info file 
Once you have downloaded the analysis_info.txt, fill in the correct arguments. **The arguments given on each line dictates which functions are run in the script.** Using the default arguments will allow you to run the pipeline from start to finish, which will mean running: GCCN-SST correction on the raw CELs, control QC, raw data QC, rma normalization, filtering (so that at least half the samples have over the median intensity value), normalized data QC, GEO submission folder creation, group comparisons and report generation. This should be suitable for most microarray analyses. In order to do this you will need to supply the correct array files, sample_info file and comparisons file (see below sections). 


### Array files   

#### QCC 
The QCC file is used to define group membership for probesets on a chip (typically different quality assessment groups), enumerate what individual probeset values should be written to CHP file headers and define textual labels for probesets for genotyping and expression arrays and should be specific to the chip/array you are using.   

#### PGF 
The probe group file (PGF) give information on what probes are contained within a probeset and information about the probes required for the analysis. 

#### CLF 
The CLF (cel layout file) maps probe IDs to a specific location in the CEL file.  


### Sample information and group comparisons
In the github Microarray folder, you will see two additional template .cvs files: "sample_info.csv" and "comparisons.csv". Download both templates and save them in the "Analysis folder", keeping their original names. 

**sample_info.csv:**
This file will have columns **SampleID** and **Group** (and **SibShip** for paired analysis, with **Batch** for batch effects). List the sample names and their corresponding group in the two columns and if the samples are paired, number them accordingly in the SibShip column (as shown in the example below).


**comparisons.csv:**
This file will have columns **baselineGroup**, **comparisonGroup**, and **Design**. The baselineGroup is normally the WT or control group. List the pairwise comparisons you want to make, using the corresponding group names. Also list the nature of the design as either  "**pairedSamples**" or "**non-pairedSamples**", depending on whether your samples are paired or un-paired. 

----------------------------------------------------------------
*The following is an example of the sample_info file:*


Sample.Name  | Sample.Group | SibShip | Batch
--- | --- | --- | ---
Sample1  | FirstGroup | 1 | 1
Sample2  | SecondGroup | 1 | 1
Sample3  | FirstGroup | 2 | 2
Sample4  | SecondGroup | 2 | 2
Sample5  | FirstGroup | 3 | 1
Sample6  | SecondGroup | 3 | 1


--------------------------------------------------------------
*The following is an example of the comparisons file:*


baselineGroup  | comparisonGroup | Design
--- | --- | ---
FirstGroup  | SecondGroup | pairedSamples

------------------------------------------------------------------



### Description of the R scripts

#### driver.R 
This script loads all the data processing and graph plotting functions from the other scripts and runs the pipeline according to the arguments given in the analysis_info.txt. By changing the arguments in the analysis info file, the driver.R should automatically run just the chosen analyses/functions.

This script should print messages on the command line to inform the user which stage of the analysis is being run and when plots and tables have been written. 

#### graph_engine.R
This is where the functions for plotting all the graphs are stored. All the arguments to run the functions should be given in the driver script. All the control and data QC plots, as well as volcano and expression plots (using the group comparison results) are created using this script.  

#### data_engine.R
This script holds the functions responsible for: reading in and loading the CEL files, transforming these using the GCCN-SST algorithm, excluding any samples specified in analysis_info.txt, applying RMA normalization, calculating batch effects and corrections and running the pairwise comparisons to give a list of regulated genes. The correct functions are called in the driver script, relating to the arguments given in the analyis info file.    

#### report_engine.R
If a report needs to be generated, this script is run and should produce an HTML output of an Rmarkdown report and an Rmd file, detailing the analysis steps, QC, data plots and results tables.

In order to send the report to a customer, who is not able to access our servers, you will need to copy the .Rmd report file, "biblio.bib" file, the "control_QC", "raw_data_plots", "results", and "rma_data_plots" directories into a separate folder e.g. ("N.Testers_reportfiles"). Then open the .Rmd file in Rstudio and press "Knit". This should produce an .html version of the analysis report into your newly created directory (N.Testers_reportfiles), using the plots within that directory.

The folder can now be zipped and sent to the customer, who should be able to see all the required tables and plots within the report and also in the separate folders, should they wish to look at them that way. 


### Standard pipeline command example

**Use the script**: bin/driver.R

**Command line example**:

```{r, engine = 'bash', eval = FALSE}
$ /usr/bin/Rscript bin/driver.R /mnt/cgs-fs3/SERVICES/Microarrays/GX-A-Tester-xxx/Analysis/analysis_info.txt
```

This will run the scripts using the arguments supplied in the "analysis_info.txt" file, which can be opened in a text editor and changed depending on which functions/outputs you require.

<br>   

### Standard output files (running the above command line)

When running the above command using the default settings in the analysis info file, the output folders produced should be those detailed below:

-------------------------------------------------------------------------------------------------------------------------
**control_QC =** *\< directory containing the RLE, NUSE, bacSpike, polyA and positive_negative control plots, as well as the signal rank table made using the control data\>*  

**GCCN_cels =** *\< directory holding the GCCN_SST transformed cels for each sample \>* 

**GEO_submission =** *\< directory containing the normalized expression matrix and a copy of the raw cels for all samples used in the analysis. This should be suitable for submission to the GEO.\>*

**raw_data_plots =** *\< directory containing the correlation heatmaps, hierarchical clustering, expression boxplots, all groups' PCA, GC content, MA and density plots produced when running the QC on the raw CEL files\>*  

**rma_data_plots =** *\< directory containing the correlation heatmaps, hierarchical clustering, expression boxplots, all groups' PCA, MA and density plots produced when running the QC on the RMA normalized CEL files\>*  

**results =** *\< directory containing the results from the pairwise comparisons, divided into sub-directories: "tables" and "graphs". The "tables" directory holds all the results tables from each pairwise comparison and an overall "summary_comparisons" .csv file, showing the number of differentially expressed genes across all comparisons. The "graphs" directory holds PCA and volcano plots for each set of pairwise comparison groups. Interactive volcano plots are also produced for each comparison group, in comparison-specific sub-directories.\>*

*\< an .Rmd and .html file will also be produced if report_gen is set to TRUE in analysis_info.txt, titled as the project name followed by "_analysisReport". \>* 

-----------------------------------------------------------------------------------------------------------------------------


