#!/bin/bash
#cd referenceData
#zip -v ./referenceData/CytoMethRefData.zip ./referenceData/hg38.fa ./referenceData/hg38.dict ./referenceData/hg38.fa.fai ./referenceData/SeqCap_EPI_CpGiant_hg38_custom_liftOver.bed ./referenceData/cpgIslandExt.hg38.bed ./referenceData/geneAnnotationEnsemble.hg38.bed
#zip -v ./referenceData/CytoMethRefDataNC_001416.zip ./referenceData/hg38_NC_001416.fa ./referenceData/hg38_NC_001416.dict ./referenceData/hg38_NC_001416.fa.fai ./referenceData/SeqCap_EPI_CpGiant_hg38_custom_liftOver_NC_001416.bed
echo "####################################"
echo "Downloading Required Reference Files"
echo "####################################"
mkdir -p ./referenceData/
while true; do
    read -p "Do you wish to download required reference files? (y/n)" yn
    case $yn in
        [Yy]* ) 
                wget -c -O ./referenceData/CytoMethRefData.zip http://zbo.ipipan.waw.pl/tools/CytoMeth/referenceData/CytoMethRefData.zip;
                unzip ./referenceData/CytoMethRefData.zip;
                rm ./referenceData/CytoMethRefData.zip;
                break;;
        [Nn]* ) echo "Skipping downloading"; break;;
        * ) echo "Please answer yes or no";;
    esac
done
echo "####################################"
echo "Downloading Optional Reference Files"
echo "####################################"
mkdir -p ./referenceData/
while true; do
    read -p "Do you wish to download optional (hg38 + NC_001416) reference files? (y/n)" yn
    case $yn in
        [Yy]* ) 
                wget -c -O ./referenceData/CytoMethRefDataNC_001416.zip http://zbo.ipipan.waw.pl/tools/CytoMeth/referenceData/CytoMethRefDataNC_001416.zip;
                unzip ./referenceData/CytoMethRefDataNC_001416.zip;
                rm ./referenceData/CytoMethRefDataNC_001416.zip;
                chmod 775 -R ./referenceData;
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
                chmod 775 -R ./input;
                break;;
        [Nn]* ) echo "Skipping downloading"; break;;
        * ) echo "Please answer yes or no";;
    esac
done
echo "########################"
echo "Downloading is Finished"
echo "########################"
