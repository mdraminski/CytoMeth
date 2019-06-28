#!/bin/bash
echo "#############################"
echo "Geting conda environment info"
echo "#############################"
conda info --json > conda.info

echo "###################################"
echo "Installation of Required R packages"
echo "###################################"
#install R packages
while true; do
    read -p "Do you wish to install or update required R CRAN and BIOCONDUCTOR packages? (y/n)" yn
    case $yn in
        [Yy]* ) Rscript R/install.packages.R; break;;
        [Nn]* ) echo "Skipping installation"; break;;
        * ) echo "Please answer yes or no";;
    esac
done

echo "#######################################"
echo "Installation of Required conda packages"
echo "#######################################"
#install tools by conda 
while true; do
    read -p "Do you wish to install or update required conda packages? (y/n)" yn
    case $yn in
        [Yy]* ) conda update conda;
                conda install -y -c bioconda bsmap;
                conda install -y -c bioconda bamtools;
                conda install -y -c bioconda bamutil;
                conda install -y -c bioconda bedtools;
                conda install -y -c bioconda seqtk;
                conda install -y -c bioconda fastqc;
                conda install -y -c bioconda/label/cf201901 samtools;
                #conda install -y -c bioconda -c r samtools --override-channels;
                break;;
        [Nn]* ) echo "Skipping installation"; break;;
        * ) echo "Please answer yes or no";;
    esac
done

echo "#######################################"
echo "Downloading of Required Reference Files"
echo "#######################################"
mkdir -p ./RefData/
while true; do
    read -p "Do you wish to download required reference files? (y/n)" yn
    case $yn in
        [Yy]* ) 
                wget -c -O ./RefData/hg38.fa http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38.fa;
                wget -c -O ./RefData/hg38.fa.fai http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38.fa.fai;
                wget -c -O ./RefData/hg38.fa.dict http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38.dict;
                wget -c -O ./RefData/hg38_phage.fa http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38_phage.fa;
                wget -c -O ./RefData/hg38_phage.fa.fai http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38_phage.fa.fai;
                wget -c -O ./RefData/hg38_phage.dict http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/hg38_phage.dict;
                wget -c -O ./RefData/SeqCap_EPI_CpGiant_hg38_custom_liftOver.bed http://zbo.ipipan.waw.pl/tools/CytoMeth/RefData/SeqCap_EPI_CpGiant_hg38_custom_liftOver.bed;
                break;;
        [Nn]* ) echo "Skipping downloading"; break;;
        * ) echo "Please answer yes or no";;
    esac
done

echo "##################################"
echo "Downloading of Simple Example Data"
echo "##################################"
mkdir -p ./input/
while true; do
    read -p "Do you wish to download input example files? (y/n)" yn
    case $yn in
        [Yy]* ) 
                wget -c -O ./input/small_FAKE03_R1.fastq http://zbo.ipipan.waw.pl/tools/CytoMeth/input/small_FAKE03_R1.fastq;
                wget -c -O ./input/small_FAKE03_R2.fastq http://zbo.ipipan.waw.pl/tools/CytoMeth/input/small_FAKE03_R2.fastq;
                break;;
        [Nn]* ) echo "Skipping downloading"; break;;
        * ) echo "Please answer yes or no";;
    esac
done

echo "########################"
echo "Installation is Finished"
echo "########################"

