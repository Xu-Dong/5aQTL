#!/bin/bash
#only use in mtcook
### docker -- Mashr erna------------------
## step1
#docker run -it --name aQTL_mashr --security-opt label:disable -t -P -h MASH -w $PWD -v $HOME:/home/$USER -v /media/london_B/zouxudong/2024-10-21-aGTEx-main/2025-05-28-mashR-eQTL/:/tmp_1 -v $PWD:$PWD -u $UID:${GROUPS[0]} -e HOME=/home/$USER -e USER=$USER gaow/hdf5tools
sos run /tmp_1/workflows/fastqtl_to_mash.ipynb \
        --data_list ./input/eQTL_file_list.txt \
        --gene_list ./input/Shared_5aQTL_genes.txt \
        --cols 4 3 \
        -j 16
