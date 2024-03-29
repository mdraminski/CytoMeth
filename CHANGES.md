CytoMeth - Changes
## Version 1.0.3
* fixed CytoMethQC script - now all data is loaded from rds results files

## Version 1.0.2
* samtools is used for sorting (multithreaded)

## Version 1.0.1
* more than one control is now available

## Version 1.0.0
* added BS-Snper tool that provides methylations and SNP
* adjusted start/end positions in the output .bed file
* fixed some minor issues
* updated documentation

## Version 0.9.21
* minor updates
* updated documentation

## Version 0.9.20
* updated reference files 
* updated documentation

## Version 0.9.19
* new parameter in config file methratio_processing
* bunch of fixes and improvements

## Version 0.9.18
* some rare issues in exporting to bed file fixed
* improved continuation of calculation after the stop
* logging on the screen improved

## Version 0.9.17
* param java_mem is now memory
* BSMAP uses I parameter (automatically set up based on available memory)
* manual updated

## Version 0.9.16
* bam files can be used as an input
* less memory usage during CalcMethylation Step

## Version 0.9.15
* reshape2 fix
* updated Dockerfile and install.packages.R files
* autodetect total memory size and cpu cores and fix the config if needed
* reference data directory change
* manual updated 

## Version 0.9.14
* new plot functions plotSitesCpG and plotSitesNonCpG allow for any min_coverage
* manual updated

## Version 0.9.13
* output files bed format - finalized structure
* new output files rds format

## Version 0.9.12
* generation of QCSummary is now independet process
* plot plotBetaValuesSummary supports data sampling and non-CpG context
* new config parameter 'ref_control_sequence_name'
* plots updated

## Version 0.9.11
* new fontsize parameter in plots

## Version 0.9.10
* picard updated
* tools.conf.yml updated
* readme file updated and extended
* output files section in readme file updated and extended

## Version 0.9.9
* picard updated
* checking of log files updated

## Version 0.9.8
* README updated
* QC extended
* install.packages.R extended by bioconductor packages
* Dockerfile added
* install.data.sh added

## Version 0.9.7
* README updated
* FastQC QC report creation added
* trimmomatic_MINLEN parameter added
* seqtk added

## Version 0.9.6
* GATK version 3.8.1
* improved error handling and missing file checking
* improved skipping of already calculated phases
* improved README
* methratio.py fixed and added (-X removed)

## Version 0.9.5
* first GitHub Version
* installation process install.sh file
* picard version 2.20.1 integrated
* conda integration - all non Java tools installed by conda
* tools.conf.yml added
* README.md and CHANGES.md files added
