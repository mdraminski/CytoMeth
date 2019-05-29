# CytoMeth
<!--- 
CytoMeth is a tool that processes methylation data designed to deal with paired-end data sequencing from Illumina. 
-->
CytoMeth tool compiles a set of open source software named in the Roche pipeline guidelines to perform SeqCap Epi data analysis. The pipe includes read quality assessment, read filtering, mapping to a reference genome, removal of PCR duplicates, assessment of coverage statistics, analyse methylation and variant calling and filtering as well as some additional functionalities added to improve the process and facilitate obtaining the processed results. Here, to obtain methylomes for brain tumor samples we used SeqCap Epi CpGiant Methylation panel and performed bisulphite conversion followed by Illumina NGS sequencing and CytoMeth tool analysis.

## Table of contents
* [Installation](#installation)
* [Usage](#usage)
* [Authors](#authors)
* [License](#license)
* [Acknowledgments](#acknowledgments)

# Installation

CytoMeth is implemented as a set of R scripts that run various tools in a specific sequence with specific set of parameters. To complete the installation process CytoMeth requires the following components installed on your OS:

- R and Rscript
- Conda - an open source package management system
- Python 2.x
- Java 8 (1.8) or above
- wget tool

If you are sure all of the above is working correctly (*R*,  *conda*, *Java*, *python 2.x*) on your system you can skip the next section and go to the section CytoMeth installation.

## Preparation of the Environment

### R and Rscript
Check if there is R installed on your machine. Type in the console window:

```bash
R --version
Rscript --version
```
If you dont have R or Rscript please install them. If missing R:
```bash
sudo apt install r-base
```
If missing Rscript only:
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
wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
```
Install it in the following directory: '*/opt/anaconda*' and remove the installation file.
```bash
bash Anaconda3-2019.03-Linux-x86_64.sh -b -p /opt/anaconda
rm Anaconda3-2019.03-Linux-x86_64.sh 
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
Now, conda should be available from the console window. Open the new one and type again:
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

## Installation of CytoMeth components

### Installation script
To install or update all required *R* and conda packages, download required reference files and set up CytoMeth,
run the '*install.sh*' script located in CytoMeth directory. For the first time select '*y*' option to install all required by CytoMeth components. All required packages and files should be installed or updated automatically and if that succeeded there is no need of any manual installation presented later below. If there is any missing component and CytoMeth stops with appropriate warning you may take a look at a specific section 'Required Tools' or 'Reference Files'. In that case please also try to rerun the script in the console window. Please notice that size of reference files is several GB and it can take a few minutes, however the downloading time strongly depends on your internet connection speed.
```bash
./install.sh
```

### Required Tools
#### Required R packages
The script '*install.sh*' file should install the following *R* packages:

- yaml
- tools
- data.table
- rjson

These packages can be also manually installed by typing the command in the console window:
```bash
Rscript R/install.packages.R
```

#### Required conda packages
The follwing *conda* packages are required by CytoMeth and these packages are automatically installed or updated by '*./install.sh*' command. The current version of CytoMeth was succesfully tested on versions presented below: 

- bsmap (ver. 2.90)
- bamtools (ver. 2.5.1)
- bamutil (ver. 1.0.14)
- bedtools (ver. 2.27)
- samtools (ver. 1.9)

These tools can be also manually installed by typing the command in the console window:
```bash
conda update conda
conda install -y -c bioconda bsmap
conda install -y -c bioconda bamtools
conda install -y -c bioconda bamutil
conda install -y -c bioconda bedtools
conda install -y -c bioconda samtools
```

#### Required Java Tools
Java tools (in CytoMeth '*/tools/*' directory) are provided with CytoMeth in the following versions:

- Trimmomatic (ver. 0.36)
- GATK (ver. 3.8.1)
- picard (ver. 2.20.1)

#### Required 'conda.info' file 
CytoMeth also requires in its main directory '*conda.info*' file that can be manually created by typing the following command in the console window: 
```bash
conda info --json > conda.info
```
This file is also created during installation process.


### Reference Files
Reference files required by CytoMeth are automatically installed by '*install.sh*' script. If you would like to download them manually plese run the following commands in console window:

```bash
wget -c -O ./RefData/hg38.fa http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38.fa
wget -c -O ./RefData/hg38.fa.fai http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38.fa.fai
wget -c -O ./RefData/hg38.fa.dict http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38.dict
wget -c -O ./RefData/hg38_phage.fa http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38_phage.fa
wget -c -O ./RefData/hg38_phage.fa.fai http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38_phage.fa.fai
wget -c -O ./RefData/hg38_phage.dict http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38_phage.dict
wget -c -O ./RefData/SeqCap_EPI_CpGiant_hg38_custom_liftOver.bed http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/SeqCap_EPI_CpGiant_hg38_custom_liftOver.bed
```
All reference files by default are located in */RefData/* directory.


# Usage

## Configuration

The file '*config.yml*' contains CytoMeth input parameters and before use of CytoMeth please define your processing. The default settings look like below:

```
#General Params
threads: 12
java_mem: 12G
overwrite_results: TRUE
clean_tmp_files: TRUE

#in/out paths
input_path: "./input/"
results_path: "./results/"
#anaconda_bin_path: "/opt/anaconda/bin/"

### Reference Data Definition
reference_sequence_file: "hg38_phage.fa"
intervals_file: "SeqCap_EPI_CpGiant_hg38_custom_liftOver.bed"
trimmomatic_adapter: "Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa"
```

Input parameters: 

- threads - define number of threads used by tools. Most of them does not gain any processing speed for more than around 10 threads.
- java_mem - amount of memory dedicated to Java tools. If you face Java 'out of memory' error please increase the parameter.
- overwrite_results - if TRUE then all result files from the sample processed again will be overwritten. If FALSE CytoMeth will skip phases that related phase result file exists in apriopriate results_path.
- clean_tmp_files - if TRUE all useless temporary files will be removed after the processing of the sample.
- plot_format - set up plot format of report files. Available formats: 'pdf','png', 'eps', 'tiff', 'jpg'.
- input_path - defines path to input fastq R1/R2 files, all samples existing in this directory will be processed in batch process.
- results_path - the path to keep all temporary and result files.
- anaconda_bin_path - path to conda and conda packages. This parameter is retrieved from 'conda.info' file and it is commented out by default. If you want to specify specific path to conda/bin directory uncommend it and define. This parameter overwrites the setting from 'conda.info' file.
- reference_sequence_file - additional control sequence file (see Input files section) by default it is set on 'hg38_phage.fa'.
- intervals_file - panel file (see Input files section)  by default it is set on SeqCap_EPI_CpGiant_hg38_custom_liftOver.bed'
- trimmomatic_adapter - trimmomatic adapter file (see Input files section) by default it is set on 'Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa'

## Input files

Before you run the processing you need to: 

- Copy your nucleotide sequences FASTA R1/R2 files (both in *.fastq* format) named: '*samplename\_R1.fastq*' and '*samplename\_R2.fastq*' (where *samplename* is a unique name of your sequenced sample) to the '*./input/*' directory. 
- Prepare reference FASTA (in .fa or .fasta format) file with additional control sequence. CytoMeth comes with 'hg38_phage.fa' reference file with an additional sequence used as control (phage DNA sequence) and the file 'hg38.fa' without that additional seqence. Any new reference '.fa' file requires corresponding '.fai' and '.dict' files that should be generated.
- Prepare panel file (in .bed format) with panel coordinates and control coordinates 'SeqCap\_EPI\_CpGiant\_hg38\_custom\_liftOver.bed''.

If you used different panel or performed whole genome analysis, please prepare the '.bed' file defining genomic regions covered by your design, to compute not biased statistics.
The control is used to estimate bisulfite conversion efficiency. Remember to add the sequence of your control e.g. enterobacteria phage lambda genome to the reference genome file so that captured controls can be mapped to the lambda genome. Reassuming, the reference genome file must be extended with a control sequence.
Important: Check if your panel file (bed format) control sequence coordinates has the same name (header ID) as in reference fasta file.


## Running the CytoMeth Processing

To run the batch processing of all samples located in *'/input/'* directory please type in the console window:
```bash
Rscript R/CytoMeth.R
```

After processing is finished create quality report on all results files located in *'/results/QC_report'* directory:
```bash
Rscript R/CytoMethQC.R
```
The script above creates summary csv file that aggregates quality measures values for all processed samples. It also creates two barplots: overall coverage report plot, CpG vs nonCpG frequency report.

The methylation results can be visualised in respect to specific genomic regions.
We annotate the level of methylation to CpG islands, promoters, intergenic regions, introns and exons.
```bash
Rscript CytoMethAnnotate.R
```

## Output files
The result files are located in *'/results/methyl_results'* directory:

- '*samplename.methylation_results.bed.panel*'
- '*samplename.methylation_results.bed.panel.no_control.CpG_min10*'
- '*samplename.methylation_results.bed.panel.no_control.non_CpG_min10*'

After the processing of each single sample CytoMeth creates a summary file associated with that sample. The default location for this summary file is '*results/QC_report*'. Example of that file for is below (created for fake data):

```
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

# Authors
This tool has been created and implemented by:

- Michal Draminski [1] (author, developer, maintainer)
- Agata Dziedzic [1] (author, developer)
- Rafal Guzik (author)
- Bartosz Wojtas [2] (author)
- Michal J. Dabrowski [1] (author)

1. Computational Biology Lab, Polish Academy of Science, Warsaw, Poland
2. Neurobiology Center, Nencki Institute of Experimental Biology, Warsaw, Poland


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
- methratio - BSD 3-Clause License https://github.com/zyndagj/BSMAPz/blob/master/methratio.py
- BisSNP GNU General Public License https://github.com/dnaase/Bis-tools/blob/master/Bis-SNP/src/main/java/edu/usc/epigenome/uecgatk/bissnp/BisSNP.java
- samtools - The MIT/Expat License https://github.com/samtools/htslib-plugins/blob/master/LICENSE

# Acknowledgments