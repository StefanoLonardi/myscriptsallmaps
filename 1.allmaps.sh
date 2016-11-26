#/bin/sh
#export CSV="CB.csv LL.csv TI.csv CI.csv SV.csv WC.csv 5I.csv"
export CSV="CB.csv LL.csv TI.csv CI.csv SV.csv WC.csv"
export BED=map.bed
export WEIGHTS=weights.txt
export CONTIGS=canu91x_high_outcov_superscaffold.fasta

pypy pacbio_lg_CB_assignment.py $CONTIGS
mv ${CONTIGS}_CB_snps_contigs.csv CB.csv
rm -rf *.nhr *.nin *.nsq Cowpea_iSelect_Ref_seq_raw.* formatdb.log

pypy pacbio_lg_CI_assignment.py $CONTIGS
mv ${CONTIGS}_CI_snps_contigs.csv CI.csv
rm -rf *.nhr *.nin *.nsq Cowpea_iSelect_Ref_seq_raw.* formatdb.log

pypy pacbio_lg_LL_assignment.py $CONTIGS
mv ${CONTIGS}_LL_snps_contigs.csv LL.csv
rm -rf *.nhr *.nin *.nsq Cowpea_iSelect_Ref_seq_raw.* formatdb.log

pypy pacbio_lg_SV_assignment.py $CONTIGS
mv ${CONTIGS}_SV_snps_contigs.csv SV.csv
rm -rf *.nhr *.nin *.nsq Cowpea_iSelect_Ref_seq_raw.* formatdb.log

pypy pacbio_lg_WC_assignment.py $CONTIGS
mv ${CONTIGS}_WC_snps_contigs.csv WC.csv
rm -rf *.nhr *.nin *.nsq Cowpea_iSelect_Ref_seq_raw.* formatdb.log

#pypy pacbio_lg_5I_assignment.py $CONTIGS
#mv ${CONTIGS}_5I_snps_contigs.csv 5I.csv
#rm -rf *.nhr *.nin *.nsq Cowpea_iSelect_Ref_seq_raw.* formatdb.log

pypy pacbio_lg_TI_assignment.py $CONTIGS
mv ${CONTIGS}_TI_snps_contigs.csv TI.csv
rm -rf *.nhr *.nin *.nsq Cowpea_iSelect_Ref_seq_raw.* formatdb.log

python -m jcvi.assembly.allmaps merge $CSV -o $BED -w $WEIGHTS
python -m jcvi.assembly.allmaps path $BED $CONTIGS -w $WEIGHTS 2> log.txt
python -m jcvi.assembly.allmaps plotall $BED --figsize=12x8 --dpi=600 --format=eps
#python -m jcvi.assembly.allmaps estimategaps $BED
#mv map.estimategaps.agp map.chr.agp
#python -m jcvi.assembly.allmaps build $BED $CONTIGS
