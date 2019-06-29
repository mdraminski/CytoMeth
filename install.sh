#!/bin/bash

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
                conda install -y -c bioconda -c r samtools --override-channels;
                break;;
        [Nn]* ) echo "Skipping installation"; break;;
        * ) echo "Please answer yes or no";;
    esac
done

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

echo "########################"
echo "Installation is Finished"
echo "########################"

