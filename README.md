# Raso lark diversity

Commands and scripts from Dierickx et al. 'Neo-sex chromosomes, genetic diversity and demographic history in the Critically Endangered Raso lark'

Note that links o raw data will be provided in the paper.
Software installation instructions are not provided.
Some custom scripts listed here that are available at: https://github.com/simonhmartin/genomics_general

## Assembly

#### Run Trimmomatic
```sh
sh trim_fastq.sh Skylark_R1_Apollo.fastq.gz Skylark_trimmed_Apollo.fastq.gz
sh trim_fastq.sh Skylark_R1_3kb.fastq.gz Skylark_trimmed_3kb.fastq.gz
sh trim_fastq.sh Skylark_2nd_R1_Apollo.fastq.gz Skylark_2nd_trimmed_Apollo.fastq.gz
sh trim_fastq.sh Skylark_2nd_R1_3kb.fastq.gz Skylark_2nd_trimmed_3kb.fastq.gz
```

#### run allpaths assembly
```sh
PrepareAllPathsInputs.pl DATA_DIR=/PATH/TO/DATA PLOIDY=2 \
IN_GROUPS_CSV=assembly_data/in_groups.csv \
IN_LIBS_CSV=assembly_data/in_libs.csv |
tee prepare.out

RunAllPathsLG PRE=/PATH/TO/PARENT/DIR \ REFERENCE_NAME=skylark \
DATA_SUBDIR=data RUN=skylark_genome_2nd HAPLOIDIFY=TRUE OVERWRITE=True THREADS=64
```



#### MUMmer alignment with zebrafinch genome


Make version of skylark genome with scaffolds larger than 1Mb

```sh
awk '($2 >= 1000000)' skylark.fa.fai | cut -f 1 > skylark.scafs1M.txt

python genomics_general/sequence.py -f skylark.scafs1M.txt < skylark.fa > skylark.scafs1M.fa

```

Alignment
```sh
# run nucmer and filter output
MUMmer3.23/nucmer --prefix=Tgut_skylark_scaf1M Tgut_chroms.fa skylark.scafs1M.fa

# filter for hits longer than 5 kb
MUMmer3.23/delta-filter -l 5000  Tgut_skylark_scaf1M.delta > Tgut_skylark_scaf1M.len5k.delta
```

Plotting
```
# plot with layout adjusted to arrange scaffolds
MUMmer3.23/mummerplot Tgut_skylark_scaf1M.len5k.delta -l -t png -s large \
-R Tgut_chroms.fa \
-Q skylark.scafs1M.fa


# NOTE This creates a gnuplot script (out.gp). Now we modify the scale. Set the size to 3000x3000 and reduce point size (ps 0.1)
sed -e 's/tiny size 1400,1400/size 3000,3000/' \
-e 's/out.png/Tgut_skylark_scaf1M.len5k.ordered.png/' \
-e 's/ps 1/ps 0.1/g'  out.gp > out.fixed.gp

# Then run gnuplot to make the final plot
gnuplot out.fixed.gp 

```
---
## Read mapping and processing

#### bowtie2
Build genome index
```sh
bowtie2-build skylark.fa skylark_genome_indexing
```

Run bowtie2
```sh
bowtie2 -x skylark_genome_indexing \
-U trimmed_sample_$i.fq \
-S sample_$i.sam
```

#### bam processing
Remove singly aligned reads
```sh
sed '/XS:/d' sample_$i.sam > sample_$i_1alignmentonly.sam
```
Convert to bam
```sh
samtools view -bS sample_$i_1alignmentonly.sam > sample_$i_1alignmmentonly.bam | samtools sort -o sample_$i_1alignmmentonly.sort.bam
```

Add readgroup info (uses Picardtools)
```sh
java -jar AddOrReplaceReadGroups.jar \
      I=sample_$i_1alignmentonly.bam \
      O=sample_$i_1alignmentonly_final.bam \
      RGID=sample_$i_1alignmentonly \
      RGLB=none RGPL=illumina RGPU=none \
      RGSM=sample_$i_1alignmentonly
```


---
## Genotyping

#### GATK HaplotypeCaller
```sh
for b in sample*1alignmentonly.sort.bam
do
python run_hapCaller.py --threads 30 -b $b -r skylark.fa --runName HC \
--gatkDir GATK/current/ --javaRam 12g -ploidy 2 --heterozygosity 0.001 --outputMode EMIT_ALL_CONFIDENT_SITES
done
```



#### Genotype GVCFs by population
```sh
for prefix in RAS ESN ESE ESW OST OSQ CRE
do
#make sample list
samples=$(awk -v prefix="$prefix" -F "\t" '$5==prefix' lark75.data.txt | cut -f 2)

#make list of vcf files
vcf_list_with_dashV=$(for s in $samples; do (echo "-V $s.sort.HC.g.vcf.gz "); done)

# run GenotypeGVCFs
java -Xmx12g -Djava.io.tmpdir=/local -jar GATK/current/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-nt 30 -R skylark.fa --includeNonVariantSites \
 $vcf_list_with_dashV \
-o $prefix.vcf.gz
done
```
#### filter VCFs
```sh
# minimal filtering at first - remove sites without at least one called genotype
for prefix in RAS ESN ESE ESW OST OSQ CRE
do
bcftools view -U $prefix.HC.vcf.gz | bgzip > $prefix.HC.clean.vcf.gz
rm $prefix.HC.vcf.gz.tbi
mv $prefix.HC.clean.vcf.gz $prefix.HC.vcf.gz
done
```

## Density of heterozygous genotypes
#### make Geno files required as input

```
for prefix in RAS ESN ESE ESW OST OSQ CRE
do
echo $prefix
bcftools filter -e 'FORMAT/DP < 5' --set-GTs "." -O u $prefix.HC.vcf.gz |
python genomics_general/VCF_processing/parseVCF.py --skipIndels |
bgzip > $prefix.HC.DP5.geno.gz &
done

#### individual heterozygosity
```sh
for pref in RAS ESN ESE ESW OST OSQ CRE
do
python genomics_general/popgenWindows.py -g $pref.HC.DP5.geno.gz \
-f phased --analysis indHet popDist -o $pref.HC.DP5.het.w250m50s50.csv \
-w 250000 -m 50 -s 50000 --roundTo 5 -T 10
done
```

#### convert scaffold positions to genome positions
```
for pref in raso26 skyN10 skyWR9 skyER10 skyCH4 ori11 crest5
do
python genomics_general/tools/transferScafPos.py -i $pref.HC.DP5.het.w250m50s50.csv \
-o $pref.HC.DP5.het.w250m50s50.chrom.csv --endCol 3 --header --sep "," --failFile fails.txt \
--transfersFile Tgut_skylark_scaf1M.len5k.order_manual.transfers.tsv
done
```
---
## Relatedness

#### NgsRelate
Extract VCF for normal recombination scaffols only

```
scafs=$(cat Raso26_HC_DP5_NRscaffolds_250Kb_Sep2019.txt)

for prefix in RAS ESN
do
echo $prefix
tabix -p vcf $prefix.HC.vcf.gz
tabix -h $prefix.HC.vcf.gz $scafs | bgzip > $prefix.HC.NRscafs.vcf.gz
done
```

Highly filtered for sites with reads in all samples (lets see if this works...), and only Normal Rec. scaffolds
```
for prefix in RAS ESN
do
echo $prefix
bcftools view -e 'FORMAT/DP < 1' $prefix.HC.NRscafs.vcf.gz | bgzip > $prefix.HC.NRscafs.DP1ALL.vcf.gz
done
```
Run NgeRelate
```
for prefix in RAS ESN
do
ngsRelate -h $prefix.HC.NRscafs.DP1ALL.vcf.gz > $prefix.HC.NRscafs.DP1ALL.ngsrelate
done
```

#### KGD
Make input
```sh
angsd -out raso26.NRscafs.min100reads.counts \
-doCounts 1 -dumpcounts 4 -setMinDepth 100 -bam raso26.bam.NRscafs.list

angsd -out skyN10.NRscafs.min50reads.counts \
-doCounts 1 -dumpcounts 4 -setMinDepth 50 -bam skylarkNL10.bam.NRscafs.list
```
Run R script
```sh
R.script GBSRun.R
```
---
## Site Frequency Spectrum

```sh
for prefix in raso26 skyN10 skyWR skyER oriental
do
angsd -gl 2 -dosaf 1 -bam $prefix.bam.NRscafs.list -ref skylark.fa -anc skylark.fa \
-baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out $prefix.NRscafs.baq1MQ1Q20GL2

realSFS $prefix.NRscafs.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > $prefix.NRscafs.baq1MQ1Q20GL2.BS20.sfs
done
```
On full genome to see how far out we would ave been with neo-sex chromosomes included

```sh
for prefix in raso26 skyN10 skyWR skyER oriental
do
angsd -gl 2 -dosaf 1 -bam $prefix.bam.list -ref skylark.fa -anc skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out $prefix.baq1MQ1Q20GL2

realSFS $prefix.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > $prefix.baq1MQ1Q20GL2.BS20.sfs
done
```


Drop-one-out tests for Raso to check consistency

```sh
#make bam lists
for i in {1..26}
do
awk -v XXX=$i 'NR!=XXX' raso26.bam.NRscafs.list > drop_one_out/raso26DO$i.bam.NRscafs.list
done

#run angsd
for i in {1..26}
do
echo $i

angsd -gl 2 -dosaf 1 -bam drop_one_out/raso26DO$i.bam.NRscafs.list \
-ref skylark.fa -anc skylark.fa  \
-baq 1 -C 50 -minMapQ 1 -P 20 -minQ 20 -out drop_one_out/raso26DO$i.NRscafsSep2019.baq1MQ1Q20GL2

realSFS drop_one_out/RAS_DO$i.NRscafsSep2019.baq1MQ1Q20GL2.saf.idx \
-P 20 -maxIter 100 > drop_one_out/raso26DO$i.NRscafsSep2019.baq1MQ1Q20GL2.sfs

rm drop_one_out/raso26DO$i.NRscafsSep2019.baq1MQ1Q20GL2.saf.*
done


#concatenate them into one
for i in {1..26}
do
cat drop_one_out/raso26DO$i.NRscafs.baq1MQ1Q20GL2.sfs >> drop_one_out/raso26DOall.NRscafs.baq1MQ1Q20GL2.sfs
done
```
Males only

```sh
for prefix in raso11M skyER5M oriental7M
do
angsd -gl 2 -dosaf 1 -bam $orefix.bam.list -ref skylark.fa -anc skylark.fa  -baq 1 -C 50 -minMapQ 1 -P 30 -minQ 20 -out $prefix.baq1MQ1Q20GL2

realSFS $prefix.baq1MQ1Q20GL2.saf.idx -P 30 -maxIter 100 -bootstrap 20 > $prefix.baq1MQ1Q20GL2.BS20.sfs
done
```
---
## Demographic analyses
#### dadi

See the Jupyter notebook `dadi_raso_final.ipynb`

#### Stairwayplot
NOTE: first step is to make the blueprint file. That was done by manualling editing the example provided with the software. The averaged folded sfs from 20 bootstraps was used.

```sh
# Then generate the batch file
java -cp stairway_plot_v2/stairway_plot_es Stairbuilder raso_two-epoch_fold.blueprint

# then run the batch file
bash raso_two-epoch_fold.blueprint.sh
```
