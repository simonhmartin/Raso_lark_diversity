################################################################################
### Normal Recombnination scaffolds, all individuals

# raso
angsd -gl 2 -dosaf 1 -bam raso.bam.NRscafsSep2018.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out raso.NRscafs.baq1MQ1Q20GL2
realSFS raso.NRscafs.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > raso.NRscafs.baq1MQ1Q20GL2.BS20.sfs

# Skylark Netherlands
angsd -gl 2 -dosaf 1 -bam skylarkNL.bam.NRscafsSep2018.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 10 -minQ 20 -out skylarkNL.NRscafs.baq1MQ1Q20GL2
realSFS skylarkNL.NRscafs.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > skylarkNL.NRscafs.baq1MQ1Q20GL2.BS20.sfs

# Skylark West Russia
angsd -gl 2 -dosaf 1 -bam skylarkWR.bam.NRscafsSep2018.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out skylarkWR.NRscafsSep2018.baq1MQ1Q20GL2
realSFS skylarkWR.NRscafsSep2018.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > skylarkWR.NRscafsSep2018.baq1MQ1Q20GL2.BS20.sfs


# Skylark East Russia
angsd -gl 2 -dosaf 1 -bam skylarkER.bam.NRscafsSep2018.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out skylarkER.NRscafs.baq1MQ1Q20GL2
realSFS skylarkER.NRscafs.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > skylarkER.NRscafs.baq1MQ1Q20GL2.BS20.sfs


# Oriental
angsd -gl 2 -dosaf 1 -bam oriental.bam.NRscafsSep2018.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out oriental.NRscafs.baq1MQ1Q20GL2
realSFS oriental.NRscafs.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > oriental.NRscafs.baq1MQ1Q20GL2.BS20.sfs

################################################################################
### All scafs all individual (for comparion with NR scaf results)

# raso
angsd -gl 2 -dosaf 1 -bam raso26.bam.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out raso26.baq1MQ1Q20GL2
realSFS raso26.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > raso26.baq1MQ1Q20GL2.BS20.sfs


# Skylark Netherlands
angsd -gl 2 -dosaf 1 -bam skylarkNL.bam.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 10 -minQ 20 -out skylarkNL.baq1MQ1Q20GL2
realSFS skylarkNL.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > skylarkNL.baq1MQ1Q20GL2.BS20.sfs


# Skylark East Russia
angsd -gl 2 -dosaf 1 -bam skylarkER.bam.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out skylarkER.baq1MQ1Q20GL2
realSFS skylarkER.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > skylarkER.baq1MQ1Q20GL2.BS20.sfs


# Oriental
angsd -gl 2 -dosaf 1 -bam oriental.bam.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out oriental.baq1MQ1Q20GL2
realSFS oriental.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > oriental.baq1MQ1Q20GL2.BS20.sfs


# Skylark West Russia
angsd -gl 2 -dosaf 1 -bam skylarkWR.bam.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out skylarkWR.baq1MQ1Q20GL2
realSFS skylarkWR.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > skylarkWR.baq1MQ1Q20GL2.BS20.sfs


################################################################################
### Drop-one-out tests for Raso to check consistency

# make bam lists
for i in {1..26}
do
awk -v XXX=$i 'NR!=XXX' raso.bam.NRscafsSep2018.list > drop_one_out/raso_DO$i.bam.NRscafsSep2018.list
done

# run angsd
for i in {1..26}
do
echo $i

angsd -gl 2 -dosaf 1 -bam drop_one_out/raso_DO$i.bam.NRscafsSep2018.list \
-ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  \
-baq 1 -C 50 -minMapQ 1 -P 10 -minQ 20 -out drop_one_out/raso_DO$i.NRscafs.baq1MQ1Q20GL2

realSFS drop_one_out/raso_DO$i.NRscafs.baq1MQ1Q20GL2.saf.idx \
-P 10 -maxIter 100 > drop_one_out/raso_DO$i.NRscafs.baq1MQ1Q20GL2.sfs

rm drop_one_out/raso_DO$i.NRscafs.baq1MQ1Q20GL2.saf.*
done


# concatenate them into one
for i in {1..12}
do
cat drop_one_out/raso_DO$i.NRscafs.baq1MQ1Q20GL2.sfs >> drop_one_out/raso_DOall.NRscafs.baq1MQ1Q20GL2.sfs
done




################################################################################
### male-only analysis for Raso, sky ER and Oriental, because they have good numbers of males

# raso (11 males)
angsd -gl 2 -dosaf 1 -bam raso11M.bam.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out raso11M.baq1MQ1Q20GL2
realSFS raso11M.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > raso11M.baq1MQ1Q20GL2.BS20.sfs


# skylark ER (5 males)

angsd -gl 2 -dosaf 1 -bam skylarkER5M.bam.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out skylarkER5M.baq1MQ1Q20GL2
realSFS skylarkER5M.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > skylarkER5M.baq1MQ1Q20GL2.BS20.sfs

# oriental (7 males)

angsd -gl 2 -dosaf 1 -bam oriental7M.bam.list -ref ../skylark_genome/skylark.fa -anc ../skylark_genome/skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out oriental7M.baq1MQ1Q20GL2
realSFS oriental7M.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > oriental7M.baq1MQ1Q20GL2.BS20.sfs
 
