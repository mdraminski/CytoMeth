#!/bin/bash
echo "#############################"
echo " Running CytoMeth Processing "
echo "#############################"
Rscript R/CytoMeth.R

echo "#############################"
echo "  Running QC Summary Report"
echo "#############################"
Rscript R/CytoMethQC.R
