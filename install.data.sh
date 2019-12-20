#!/bin/bash
# zip -lv ./referenceData/CytoMethReferenceData.zip ./referenceData/*
echo "#######################################"
echo "Downloading of Required Reference Files"
echo "#######################################"
mkdir -p ./referenceData/
while true; do
    read -p "Do you wish to download required reference files? (y/n)" yn
    case $yn in
        [Yy]* ) 
                wget -c -O ./referenceData/CytoMethReferenceData.zip http://zbo.ipipan.waw.pl/tools/CytoMeth/referenceData/CytoMethReferenceData.zip;
                unzip ./referenceData/CytoMethReferenceData.zip;
                rm ./referenceData/CytoMethReferenceData.zip;
                break;;
        [Nn]* ) echo "Skipping downloading"; break;;
        * ) echo "Please answer yes or no";;
    esac
done

echo "##################################"
echo "Downloading of Basic Example Data"
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
echo "Downloading is Finished"
echo "########################"
