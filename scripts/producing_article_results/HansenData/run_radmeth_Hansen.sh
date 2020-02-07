#!/bin/bash -l

cd /path/to/input/data/folder

#module load GCCcore/6.3.0
module load GSL/2.4-goolf-triton-2017a
SCRIPTFOLDER=/path/to/scripts
INPUTFOLDER=/path/to/input/data/folder
RADMETHPATH=/path/to/methpipe/software


"$RADMETHPATH"/methpipe-3.4.3/bin/radmeth regression -factor IsCase "$SCRIPTFOLDER"/design_matrix_fullrank_radmeth.txt "$INPUTFOLDER"/proportion_table_Hansen_wHeader.txt > "$INPUTFOLDER"/radmeth_cpgs_Hansen.bed
"$RADMETHPATH"/methpipe-3.4.3/bin/radmeth adjust -bins 1:200:1 "$INPUTFOLDER"/radmeth_cpgs_Hansen.bed > "$INPUTFOLDER"/radmeth_cpgs_Hansen.adjusted.bed
