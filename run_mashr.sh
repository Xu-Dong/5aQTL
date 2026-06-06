#!/bin/bash

# global settings
basedir=/media/london_B/zouxudong/2022-08-05-altTSS-QTL-Project/2026-03-21-mashr-5aQTL

main(){
    run_prepare_eQTL_in_mashr_format
    run_prepare_gene_list
    run_mashr
}

function run_prepare_eQTL_in_mashr_format(){
    aqtlDir=$basedir/input/QTL_fold
    mkdir -p $basedir/input/eQTL_mashr

    for f in `ls $aqtlDir/*.cis_QTL_all.txt.gz`
    do
            Tissue=`basename $f .cis_QTL_all.txt.gz`
            echo $Tissue
            zcat $f|awk 'BEGIN{OFS="\t"}{print $1,$8,$12,$14}'|bgzip -c > $basedir/input/eQTL_mashr/${Tissue}.txt.gz &
        done

}


function run_prepare_gene_list(){
        eqtlDir=$basedir/input/eQTL_mashr
        module load languages/R-4.1.3
        mkdir -p $basedir/input/gene_list
        for f in `ls $eqtlDir/*.txt.gz`
        do
                Tissue=`basename $f .txt.gz`
                echo $Tissue
                zcat $f|cut -f1|sort|uniq|awk -v TISSUE=$Tissue '{print $0"\t"TISSUE}' > $basedir/input/gene_list/${Tissue}.5aQTL_tested_gene.txt &
        done
        echo "get union gene list..."
        Rscript $basedir/src/get_union_genelist.R

}



function run_mashr(){
#      step 1
        bash mashr_step1.sh
#      step 2
        bash mashr_step2.sh
}


main
