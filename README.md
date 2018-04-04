---
title: "Microarray pipeline"
output:
  github_document:
    html_preview: true
    toc: true
    toc_depth: 4
---

## Getting Started
1. Create a new folder in cgs-fs3 or any other location with your chosen project name - referred to as the "project folder" in this README.   
2. Create a "bin/" folder in the project folder.   
3. Download and copy the scripts to the "bin/" folder. 


## Running the pipeline

### 1. Analysis info file
A central part of the pipeline is the **analysis info** file. It has information about the project location (main working directory), the original CEL files folder, the QCC file folder and additional parameters used throughout the analysis. 

#### Format of the analysis info file
The analysis info file is a simple .txt file, with each line providing information. Different parameters on the same line are separated by the semicolon (i.e ";") character.   

The following is an example of the analysis info file:   

--------------------------------------------------------------------------------------------------   
**project_location =** /mnt/cgs-fs3/SERVICES/Microarrays/GX-Affy/GX-A.Tester-10235 
**rawCELs_folder =** /mnt/cgs-fs3/SERVICES/Microarrays/GX-Affy/GX-A.Tester-10235/cels 
**QCC_folder =** /mnt/cgs-fs3/SERVICES/Microarrays/GX-Affy/GX-A.Tester-10235/QCC_file
--------------------------------------------------------------------------------------------------   

The following is the explanation of the lines within the analysis info file:   

--------------------------------------------------------------------------------------------------   
**project_location =** *\<path to directory of the analysis\>*   
**rawCELs_folder =** *\<path to the raw CEL files folder\>*    
**QCC_folder =** *\<path to the QCC file folder\>*    
--------------------------------------------------------------------------------------------------    
   
<br>   

#### QCC file   
The QCC file is used to define group membership for probesets on a chip (typically different quality assessment groups), enumerate what individual probeset values should be written to CHP file headers and define textual labels for probesets for genotyping and expression arrays and should be specific to the chip/array you are using.    


#### How to create the analysis info file

**Use the script**: bin/create_analysisinfo_file.py

**Command line example**:

```{r, engine = 'bash', eval = FALSE}
$ python bin/create_analysisinfo_file.py
```


This will create an "analysis_info.txt" file, which you can open in a text editor and fill, with the project specific paths and parameters.

<br>   
