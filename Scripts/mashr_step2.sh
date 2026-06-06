#!/bin/bash
#only use in mtcook
### docker -- Mashr erna------------------
baseDir=`pwd`
## step2 
docker run -it --name eQTL_mashr2 --security-opt label:disable -t -P -h MASH -w $PWD -v $HOME:/home/$USER -v $baseDir/:/tmp -v $PWD:$PWD -u $UID:${GROUPS[0]} -e HOME=/home/$USER -e USER=$USER gaow/mash-paper

indir=$baseDir
sos run /tmp/workflows/mashr_flashr_workflow.ipynb posterior \
        --vhat identity \
        --data $indir/fastqtl_to_mash_output/eQTL_file_list.mash.rds
