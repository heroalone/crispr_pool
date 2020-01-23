# crispr_pool


FQfile="trimed_ABC.fastq"   ### trimed reads for plasmid pools
TARGETfile="example_sgRNA_targets.txt"   ### target information

perl step1_Get_sgRNA_fromPool.pl $FQfile
perl step1_Get_sgRNA_fromPool.pl $FQfile $TARGETfile


#### final outputs ####
example_fullymatch.txt # reads statistics that fully covered sgRNA seqs
## header:GeneID	sgRNA	covered_reads_number

example_mismatch.txt # reads that didn't fully covered sgRNA seqs, but meet the given mismatch conditions, 80% is set in the example script.
## GeneID	ID	sgRNA	identity(%)	mimatched_seq

example_nomatch.txt #  reads that didn't meet the given mismatch conditions, but the flanking sequences are fully matched to the barcodes. 
## readsID	seq_without_matched_sgRNA
