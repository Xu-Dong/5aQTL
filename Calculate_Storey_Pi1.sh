#!/bin/bash
dir=`pwd`

# prepare 5'aQTL file in 49 tissues
Rscript $dir/Scripts/prepare_data_for_pi1.R

# calculate pi1
for TISSUE in `cat $dir/GTEx_v8_tissues.txt`
do
    Rscript $dir/Scripts/calc_pi1_index.R -t $TISSUE
done
