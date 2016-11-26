#/bin/sh
export CSV=TI.csv
export BED=TI.bed
export WEIGHTS=weights.txt
export CONTIGS=pilon62x_fast_allreads_superscaffold.fasta
python -m jcvi.assembly.allmaps merge $CSV -o $BED -w $WEIGHTS
python -m jcvi.assembly.allmaps path $BED $CONTIGS -w $WEIGHTS 2> log.txt
#python -m jcvi.assembly.allmaps estimategaps $BED
#mv map.estimategaps.agp map.chr.agp
#python -m jcvi.assembly.allmaps build $BED $CONTIGS
