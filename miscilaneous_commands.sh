################################################################################
# count sites with a given depth of coverage in each sample

#(The output requires some post-processing to add column headers and sort etc).
samtools depth -f ../vcf/bam_list.txt  | cut -f 3- |  countThings.py > coverage_counts.txt

################################################################################

# make version of genome with only scaffolds >= 1 Mb

awk '($2 >= 1000000)' skylark.fa.fai | cut -f 1 > skylark.scafs1M.txt

python ~/Research/genomics_general/sequence.py -f skylark.scafs1M.txt < skylark.fa > skylark.scafs1M.fa

###############################################################################
### NGSrelate
ngsRelate -h ../vcf/raso26.HC.NRscafs.DP1ALL.vcf.gz > raso26.HC.NRscafs.DP1ALL.ngsrelate 

ngsRelate -h ../vcf/skyN13.HC.NRscafs.DP1ALL.vcf.gz > skyN13.HC.NRscafs.DP1ALL.ngsrelate

###############################################################################
### KGD

# counts of bases per site per individual
angsd -out raso26.NRscafs.min100reads.counts -doCounts 1 -dumpcounts 4 -setMinDepth 100 -bam raso.bam.NRscafsSep2018.list

angsd -out skyN13.NRscafs.min50reads.counts -doCounts 1 -dumpcounts 4 -setMinDepth 50 -bam skylarkNL.bam.NRscafsSep2018.list

#then run R script GBSRun.R

###############################################################################
