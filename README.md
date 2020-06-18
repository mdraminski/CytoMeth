---
output:
  pdf_document: default
  word_document: default
  html_document: default
---
# CytoMeth
<!--- 
CytoMeth is a tool that processes methylation data designed to deal with paired-end data sequencing from Illumina. 
-->
CytoMeth tool compiles a set of open source software named in the Roche pipeline guidelines to perform SeqCap Epi data analysis. The pipe includes read quality assessment, read filtering, mapping to a reference genome, removal of PCR duplicates, assessment of coverage statistics, analyse methylation and variant calling and filtering as well as some additional functionalities added to improve the process and facilitate obtaining the processed results. Here, to obtain methylomes for brain tumor samples we used SeqCap Epi CpGiant Methylation panel and performed bisulphite conversion followed by Illumina NGS sequencing and CytoMeth tool analysis.

## Table of contents
* [Installation](#installation)
  * [Environment Preparation](#environment-preparation)
  * [Installation of CytoMeth Components](#installation-of-cytometh-components)
  * [The Docker](#the-docker)
* [CytoMeth Usage](#usage)
* [Authors](#authors)
* [License](#license)
* [Acknowledgments](#acknowledgments)

# Installation
CytoMeth is implemented as a set of R scripts that run various tools in a specific sequence with specific set of parameters. It can be installed (and used) in two ways:

- as a set of scripts and third party required tools installed directly in your Linux OS
- as a docker image that can be run under any OS

If you prefer the docker installation please skip the section below and go to the section [The Docker]. If you prefer to install it directly on your Linux environment please go through all below steps.

## Environment Preparation
Notice that below steps refer to Linux OS. However more experienced user can also use this manual to install CytoMeth on OSX but it requires some slight adaptations that are not provided in the description.

To complete the installation process CytoMeth requires the following components installed on your OS:

- R and Rscript
- Conda - an open source package management system
- Python 2.x
- Java 8 (1.8) or above
- wget tool

If you are sure all of the above is working correctly (*R*,  *conda*, *Java*, *python 2.x*) on your system you can skip the next section and go to the section [Installation of CytoMeth Components].

### R and Rscript
Check if there is R installed on your machine. Type in a terminal window:

```bash
R --version
Rscript --version
```
If you dont have R or Rscript please install them. If missing R:
```bash
sudo apt install r-base
```
If still missing Rscript try to install:
```bash
sudo apt-get install littler
```

### Conda
Check if there is conda installed on your machine:
```bash
conda info
```
If conda *command not found* please install it. Download anaconda from the web:
```bash
wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
```
Install it in the following directory: '*/opt/anaconda*' and remove the installation file.
```bash
sudo bash Anaconda3-2020.02-Linux-x86_64.sh -b -p /opt/anaconda
rm Anaconda3-2020.02-Linux-x86_64.sh 
```
Create new users group 'anaconda' and give all required priviliges to that group.
```bash
sudo groupadd anaconda 
sudo chgrp -R anaconda /opt/anaconda
sudo chmod 770 -R /opt/anaconda
```
Add all CytoMeth users to the anaconda group. Please replace *\_user\_* with your username. These users will have an access to all tools installed by conda.
```bash
sudo adduser _user_ anaconda 
```
All CytoMeth users need to add path to anaconda to the *PATH* system variable in '*.bashrc*' file:
```bash
echo 'export PATH="/opt/anaconda/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```
Now, conda should be available from a terminal window. Open the new one and type again:
```bash
conda info
```
You should see all information about conda environment.

### Python 2.x
Check your python version. It is recommended to use python 2.x with CytoMeth.
```bash
python2 --version
```
If your OS lacks of python 2.x please install it with the following commnads:
```bash
sudo apt-get update
sudo apt-get install python2
```

### Java
Check Java version. It is recommended to use Java >=8 (1.8) with CytoMeth.
```bash
java -version
```
If your OS lacks of Java please install it with the following commnads:
```bash
sudo apt-get update
sudo apt-get install openjdk-8-jre-headless
```

## Installation of CytoMeth Components

To get CytoMeth from the github repository you may download it as a zip file or clone the project:
```bash
git clone https://github.com/mdraminski/CytoMeth.git
cd CytoMeth
```

### Installation script
To install or update all required *R* and conda packages, download required reference files and set up CytoMeth,
run '*install.sh*' and '*install.data.sh*' scripts both located in CytoMeth directory. The first one installs all required R and conda packages and the second one downloads all required reference files and basic example data. For the first time select '*y*' option to install all required by CytoMeth components. All required packages and files should be installed or updated automatically and if that succeeded there is no need of any manual installation presented later below. If there is any missing component and CytoMeth stops with appropriate warning you may take a look at a specific section 'Required Tools' or 'Reference Files'. In that case please also try to rerun the script in a terminal. Please notice that size of reference files is several GB and it can take a few minutes to download them, however the downloading time strongly depends on your internet connection speed.
```bash
./install.sh
./install.data.sh
```

### Required Tools
#### Required R packages
The script '*install.sh*' file should install the following *R* packages:

- data.table (CRAN package)
- ggplot2 (CRAN package)
- rjson (CRAN package)
- RColorBrewer (CRAN package)
- yaml (CRAN package)
- stringr (CRAN package)
- methylKit (Bioconductor package)
- GenomicRanges (Bioconductor package)
- genomation (Bioconductor package)

These packages can be also manually installed by typing the command in a terminal window:
```bash
Rscript R/install.packages.R
```

#### Required conda packages
The follwing *conda* packages are required by CytoMeth and these packages are automatically installed or updated by '*./install.sh*' command. The current version of CytoMeth was succesfully tested on versions presented below: 

- bsmap (ver. 2.90)
- bamtools (ver. 2.5.1)
- bamutil (ver. 1.0.14)
- bedtools (ver. 2.27)
- seqtk (ver. 1.3)
- fastqc (ver. 0.11.8)
- samtools (ver. 1.9)

These tools can be also manually installed by typing the command in the a terminal window:
```bash
conda update conda
conda update conda-build
conda install -y -c bioconda bsmap
conda install -y -c bioconda bamtools
conda install -y -c bioconda bamutil
conda install -y -c bioconda bedtools
conda install -y -c bioconda seqtk
conda install -y -c bioconda fastqc
conda install -y -c bioconda samtools
```

#### Required Java Tools
Java tools (in CytoMeth '*/tools/*' directory) are provided with CytoMeth with the following versions:

- Trimmomatic (ver. 0.36)
- GATK (ver. 3.8.1)
- picard (ver. 1.141) or picard2 (ver. 2.20.6)

#### Required 'conda.info' file 
CytoMeth also requires in its main directory '*conda.info*' file that can be manually created by typing the following command in a terminal window: 
```bash
conda info --json > conda.info
```
This file is also automatically created during installation process.


### Reference Files
Reference files required by CytoMeth are automatically installed by '*install.data.sh*' script. If you would like to download them manually plese run the following commands in a terminal window:

```bash
wget -c -O ./referenceData/CytoMethReferenceData.zip http://zbo.ipipan.waw.pl/tools/CytoMeth/referenceData/CytoMethReferenceData.zip;
unzip ./referenceData/CytoMethReferenceData.zip;
```
All reference files are located in */referenceData/* directory by default.

### Basic Example Data
Basic example data may be downloaded by wget command:
```bash
wget -c -O ./input/small_FAKE03_R1.fastq http://zbo.ipipan.waw.pl/tools/CytoMeth/input/small_FAKE03_R1.fastq;
wget -c -O ./input/small_FAKE03_R2.fastq http://zbo.ipipan.waw.pl/tools/CytoMeth/input/small_FAKE03_R2.fastq;
```

## The Docker
CytoMeth project is also available as a docker. The CytoMeth docker is a virtual machine that contains all the environment (apps and libraries) ready to run CytoMeth. To download and run CytoMeth docker please install Docker app from https://www.docker.com/. After successfull instalation of Docker app you may build your own CytoMeth docker from the sources or download ready to use CytoMeth docker from  Docker Hub [Downloading the docker from Docker Hub]. 

### Building your own docker locally
To build your own docker use *build* command and after the successful creation the docker is ready to run. Pleaese notice building of the docker may take tens of minutes because the proper environment must be created from the scratch, however it must be done only once.

To get CytoMeth from the github repository you may download it as a zip file or clone the project:
```bash
git clone https://github.com/mdraminski/CytoMeth.git
cd CytoMeth
```

To build your own docker run the command in the CytoMeth directory:
```bash
docker build -t cytometh .
```

### Downloading the docker from Docker Hub
The docker ready to go is also publicly available on Docker Hub and can be pulled to your system by the command below. The second command adds a new docker tag so the name of the local docker image is cytometh instead of inconvenient e.g. mdraminski/cytometh:2.
```bash
docker pull mdraminski/cytometh:2
docker tag mdraminski/cytometh:2 cytometh
```

### Running the docker
To run the docker that is already built in your system or pulled from Docker Hub type the command below. Please notice that you are running fresh virtual machine session and this image does not yet contain '*input*' and '*referenceData*' folders. However they may be created and filled by '*install.data.sh*' script (See below '*Reference data*' section).
```bash
docker run -it cytometh /bin/bash
```
Notice all data that you download or create (e.g. results) within the docker session is available until its shut down. Therefore it is highly recomennded to share the folder between the docker and the host system (for data and results transfer). To run the docker that shares the folder between host system and the docker it is needed to specify it right after '-v' parameter e.g. to share Desktop folder in your home folder run the command below: 
```bash
docker run -it -v ~/Desktop:/Desktop cytometh /bin/bash
```

Please remember if you want to share the folder '*Desktop*' between your docker and local system you need to modify the following '*config.yml*' paths accordingly:
```yaml
input_path: "/Desktop/input/"
results_path: "/Desktop/results/"
ref_data_path: "/Desktop/referenceData/"
```

### Quit from the docker
To shut down the virtual machine type command 'exit'. It is similar as quiting from ssh session.

### Reference data
Reference data is several Gigabytes big therefore it is not included in any parent docker. However after successful running of the docker you can download the data by running '*install.data.sh*' script in the CytoMeth main directory. 
```bash
./install.data.sh
```
After successful installation of the reference data CytoMeth docker is ready to use. All reference files are located in */referenceData/* directory by default.

# CytoMeth Usage

## Configuration

### File 'config.yml'
The file '*config.yml*' contains CytoMeth input parameters and before use of CytoMeth please define your processing. The default settings look like below:

```yaml
#General Params
verbose: TRUE
threads: 8
memory: 16G
overwrite_results: FALSE
clean_tmp_files: TRUE
plot_format: "pdf"

#in/out paths
input_path: "./input/"
results_path: "./results/"
#anaconda_bin_path: "/opt/anaconda/bin/"

### Reference Data - Path
ref_data_path: "./referenceData/"
### Reference Data - Files
reference_sequence_file: "hg38_phage.fa"
intervals_file: "SeqCap_EPI_CpGiant_hg38_custom_liftOver.bed"
ref_sequence_name: "NC_001416"

# Specific Tools params
trimmomatic_MINLEN: 50
sqtk_run: FALSE
sqtk_subset: 10000000
```

Input parameters:

- verbose - prints additional info and commands on the screen
- threads - defines number of threads used by tools. Most of the tools does not gain any processing speed for more than 10-12  threads.
- memory - amount of memory dedicated to Java and other tools. If you see Java 'out of memory' error or any sudden stop of the program please increase the parameter. The minimum amount that is recommended for human genome analysis is 6GB. Plese use one of the following sufixes: 'M', 'G', 'T' (case sensitive: mega, giga, tera).
- overwrite_results - if TRUE then all result files from the sample processed again will be overwritten. If FALSE CytoMeth will skip phases that related phase result file exists in apriopriate results_path.
- clean_tmp_files - if TRUE all useless temporary files will be removed after the processing of the sample.
- plot_format - set up plot format of report files. Available formats: 'pdf','png', 'eps', 'tiff', 'jpg'.
- input_path - defines path to input fastq R1/R2 files, all samples existing in this directory will be processed in batch process.
- results_path - the path to keep all temporary and result files.
- anaconda_bin_path - path to conda and conda packages. This parameter is retrieved from 'conda.info' file and it is commented out by default. If you want to specify specific path to conda/bin directory uncommend it and define. This parameter overwrites the setting from 'conda.info' file.
- ref_data_path - defines path to the reference data
- ref_data_sequence_file - additional control sequence file (see Input files section) by default it is set on 'hg38_phage.fa'.
- ref_data_intervals_file - panel file (see Input files section)  by default it is set on SeqCap_EPI_CpGiant_hg38_custom_liftOver_phage.bed'
- ref_control_sequence_name - name of control sequence (by default phage sequence)
- trimmomatic_MINLEN - MINLEN parameter of the trimmomatic tool.
- sqtk_run - if TRUE initial sqtk sampling is processed.
- sqtk_subset - size of the subset to select by the sqtk tool.

### File 'tools.conf.yml'
The file '*tools.conf.yml*' contains CytoMeth tools parameters and it is located in tools directory. The settings in the file configure paths and names of all tools needed by CytoMeth processing default values are highly recommended. The file by default is defined as below:

```yaml
### TOOLS - path and tools cfg file
tools_path: "./tools/"
tools_config: "tools.conf.csv"

### TOOLS Definition 
bedtools: "bedtools"
samtools: "samtools"
bamtools: "bamtools"
bamUtil: "bam"
bsmap: "bsmap"
methratio: "methratio/methratio.py"
trimmomatic: "Trimmomatic/trimmomatic-0.36.jar" 
picard: "Picard/picard.jar"
picard_ver: 1
gatk: "GATK/GenomeAnalysisTK.jar"
fastqc: "fastqc"
seqtk: "seqtk"
bisSNP:

### python2 path
python2: "python2"

### Reference Data - Remaining Files
ref_data_trimmomatic_adapter: "Trimmomatic/adapters/TruSeq3-PE-2.fa"
ref_data_CpgIslandAnnotation: "cpgIslandExt.hg38.bed"
ref_data_CpGGenomAnnotation: "geneAnnotationEnsemble.hg38.bed"
```

## Input files

Before you run the processing you need to: 

- Copy your nucleotide sequences FASTA R1/R2 files (both in *.fastq* format) named: 'SAMPLENAME\_R1.fastq' and 'SAMPLENAME\_R2.fastq' (where 'SAMPLENAME' is a unique name of your sequenced sample) to the '*./input/*' directory.
- If your files are compressed (*.gz* format) please decompress before use:
  ```bash
  gunzip -c SAMPLENAME_R1.fastq.gz > SAMPLENAME_R1.fastq
  gunzip -c SAMPLENAME_R2.fastq.gz > SAMPLENAME_R2.fastq
  ```
- Prepare reference FASTA (in .fa or .fasta format) file with additional control sequence. CytoMeth comes with 'hg38_phage.fa' reference file with an additional sequence used as control (phage DNA sequence) and the file 'hg38.fa' without that additional seqence. Any new reference '.fa' file requires corresponding '.fai' and '.dict' files that should be generated. The control is used to estimate bisulfite conversion efficiency. Remember to add the sequence of your control e.g. enterobacteria phage lambda genome to the reference genome file so that captured controls can be mapped to the lambda genome. Reassuming, the reference genome file must be extended with a control sequence. Notice that all reference files are located in */referenceData/* directory.
- Prepare panel file (in .bed format) with panel coordinates and control coordinates '*SeqCap\_EPI\_CpGiant\_hg38\_custom\_liftOver\_phage.bed*'. If you used different panel or performed whole genome analysis, please prepare the '.bed' file defining genomic regions covered by your design, to compute not biased statistics. **Important**: Check if your panel file (bed format) control sequence coordinates has the same name (header ID) as in reference fasta file.

## Running the CytoMeth Processing

To run entire CytoMeth processing and summary reporting please run '*CytoMethRun.sh*' bash script in a terminal window:
```bash
./CytoMethRun.sh
```

It is also possible to run the batch processing separately for all samples located in *'/input/'* directory. If it is required  please type in a terminal window:
```bash
Rscript R/CytoMeth.R
```

When above processing is finished create summary quality report on all results files located in *'/results/QC_report'* directory:
```bash
Rscript R/CytoMethQC.R
```
The script above creates summary csv file that aggregates quality measures values for all processed samples. It also creates two barplots: overall coverage report plot, CpG vs nonCpG frequency report. The methylation results can be also visualised in respect to specific genomic regions.
We annotate the level of methylation to CpG islands, promoters, intergenic regions, introns and exons and provide proper plots in '*results/QC_report*' directory.

It is also possible to define your own processing chain and run multiple experiments on different input and output folders or different set of input parameters. To set up the CytoMeth process manually plese edit '*CytoMeth.R*' file.
```R
source("./R/main.R")
#read default config from the config.yml file
conf <- readConfig()
#set up required parameters e.g. input path
conf$input_path <- "./myinputfolder/"
conf$overwrite_results <- F

# run batch processing of all files located in conf$input_path folder
CytoMeth(conf)

# For single sample processing it is required to define R1 and R2 files
CytoMeth(conf, file.path(conf$input_path,"small_FAKE03_R1.fastq"), 
  file.path(conf$input_path,"small_FAKE03_R2.fastq"))
```

## Output files

### Methylation beta values
The result files are located in *'/results/methyl_results'* directory. For each sample there are three output files:

- 'SAMPLENAME.methylation_results.bed.panel'
- 'SAMPLENAME.methylation_results.bed.panel.no_control.CpG_min10'
- 'SAMPLENAME.methylation_results.bed.panel.no_control.non_CpG_min10'

### FastQC report files
The result files are located in *'/results/QC/FastQC'* directory. For each sample there are four output files:

- 'SAMPLENAME_R1_fastqc.zip'
- 'SAMPLENAME_R2_fastqc.zip'
- 'SAMPLENAME_R2_fastqc.html'
- 'SAMPLENAME_R2_fastqc.html'

There is no need to unzip these files, FastQC report is available by opening html file in the browser.

### CytoMeth Quality Control report files
After the processing of each single sample CytoMeth creates a summary file associated with that sample. The default location for this summary file is '*results/QC_report/SAMPLENAME_QC_summary.yml*'. Example of that *yml* file is below (created for fake data):

```bash
Sample_ID: small_FAKE03
Input_read_pairs: 289087.0
Read_Pairs_Surviving_trimming: 289086.0
Prc_Read_Pairs_Survaving_trimming: 100.0
Prc_duplicated_reads_top: 7.5427
Prc_duplicated_reads_bottom: 7.7892
Number_of_reads_after_removing_duplicates: 533815.0
Number_of_reads_after_filtering: 474034.0
Prc_passed_filtering_step: 88.8011764
Number_of_on_target_reads: 533906.0
Prc_of_on_target_reads: 100.0170471
Mean_coverage: 0.41
Number_of_Cs_in_control: 0
Conversion_eff: .nan
Number_of_Cs_in_panel: 295448
Number_of_Cs_in_panel_CpG: 42450
Number_of_Cs_in_panel_non_CpG: 252998
Number_of_Cs_in_panel_non_CpG_cov_min10: 189937
Number_of_Cs_in_panel_CpG_cov_min10: 33856
Number_of_Cs_in_panel_CpG_cov_max9: 8594
Prc_of_Cs_in_panel_CpG_cov_min10: 79.7550059
Prc_of_Cs_in_panel_CpG_cov_max9: 20.2449941
```

### CytoMeth Quality Control summary file
The file '*results/QC_report/SummaryQC.csv*' contains all QC report files results for all samples.

### CytoMeth Quality Control plots 
The directory '*results/QC_report/*' contains set of plot files:

- 'SAMPLENAME_histCpGStats.pdf'
- 'BetaValuesSummary.pdf'
- 'SummaryCntCommonCpG.pdf'
- 'SummaryCpGCoverage.pdf'
- 'SummaryHypermethylatedCpGGenomAnnotation.pdf'
- 'SummaryHypermethylatedCpGIslandsAnnotation.pdf'
- 'SummaryHypomethylatedCpGGenomAnnotation.pdf'
- 'SummaryHypomethylatedCpGIslandsAnnotation.pdf'
- 'SummaryMethylationLevel.pdf'
- 'SummarySitesCovBy10CpGnonCpG.pdf'
- 'SummarySitesCovBy10.pdf'

# Version
For more information please see CHANGES.md

- Version: 0.9.16
- Date: 17.01.2020

# Authors
This tool has been created and implemented by:

- Michal Draminski [1] (author, developer, maintainer)
- Agata Dziedzic [1] (author, developer)
- Rafal Guzik [3] (author)
- Bartosz Wojtas [2] (author)
- Michal J. Dabrowski [1] (author)

1. Computational Biology Lab, Polish Academy of Science, Warsaw, Poland
2. Neurobiology Center, Nencki Institute of Experimental Biology, Warsaw, Poland
3. Andrzej Frycz Modrzewski Krakow University, Faculty of Medicine and Health Sciences, Department of Biochemistry, Poland.

# License
This program and the accompanying materials are made available under the terms of the GNU Public License v3.0 which accompanies this distribution, and is available at [http://www.gnu.org/licenses/gpl.html](http://www.gnu.org/licenses/gpl.html). For more information please see LICENSE file.

The set of tools, used by or provided with CytoMeth is under following licenses:

- fastqc - GNU General public license https://www.bioinformatics.babraham.ac.uk/projects/download.html
- trimmomatic - GNU GENERAL PUBLIC LICENSE https://github.com/timflutre/trimmomatic/blob/master/distSrc/LICENSE
- bsmap - GNU Public License (GPL) https://github.com/genome-vendor/bsmap/blob/master/README.txt
- picard - The Picard toolkit is open-source under the MIT license and free for all uses. https://broadinstitute.github.io/picard/
- bamtools - The MIT License https://github.com/pezmaster31/bamtools/blob/master/LICENSE
- bamUtil - GNU General Public License https://github.com/statgen/bamUtil/blob/master/src/Validate.h
- GATK - BSD 3-Clause "New" or "Revised" License https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT
- BisSNP GNU General Public License https://github.com/dnaase/Bis-tools/blob/master/Bis-SNP/src/main/java/edu/usc/epigenome/uecgatk/bissnp/BisSNP.java
- samtools - The MIT/Expat License https://github.com/samtools/htslib-plugins/blob/master/LICENSE

# Acknowledgments

- We would like to thank to everyone who helped and contributed to this work. 
- This work is supported by a grant from the Polish National Science Centre [DEC-2015/16/W/NZ2/00314].
