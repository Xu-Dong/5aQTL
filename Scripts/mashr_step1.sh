#!/bin/bash
#only use in mtcook
### docker -- Mashr------------------
## step1
baseDir=`pwd`
docker run -it --name aQTL_mashr --security-opt label:disable -t -P -h MASH -w $PWD -v $HOME:/home/$USER -v $baseDir/:/tmp_1 -v $PWD:$PWD -u $UID:${GROUPS[0]} -e HOME=/home/$USER -e USER=$USER gaow/hdf5tools
sos run /tmp_1/workflows/fastqtl_to_mash.ipynb \
        --data_list ./input/eQTL_file_list.txt \
        --gene_list ./input/Shared_5aQTL_genes.txt \
        --cols 4 3 \
        -j 16
