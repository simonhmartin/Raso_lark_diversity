################################################################################
#SAMtools

samtools mpileup -f ../skylark_genome/skylark.fa --VCF -t DP --skip-indels -u --bam-list lark78.bam.list | bgzip > lark78.ST.vcf.gz

bcftools call -m lark78.ST.vcf.gz | bcftools filter -e 'FORMAT/DP<5' --set-GTs . | bcftools view -e 'AN<2' | bgzip > lark78.ST.DP5AN2.vcf.gz

python ~/Research/genomics_general/VCF_processing/parseVCF.py -i lark78.ST.DP5AN2.vcf.gz | bgzip > lark78.ST.DP5AN2.geno.gz


# To count sites and SNPs, keep only sites that have 100 reads across the dataset (this is arbitrary!)

bcftools filter -e 'DP < 100' -O u lark78.ST.vcf.gz | bcftools view -i 'TYPE=="snp" | TYPE=="ref"' | bgzip > lark78.ST.100reads.vcf.gz

# and a SNPs only version
bcftools filter -i 'MAC >= 1' lark78.ST.100reads.vcf.gz | bgzip > lark78.ST.100reads.VAR.vcf.gz

# get number of sites and number of SNPs

zcat lark78.ST.100reads.vcf.gz | grep -v "#" | wc -l
zcat lark78.ST.100reads.VAR.vcf.gz | grep -v "#" | wc -l


################################################################################
#### GATK hapcaller

for b in sample*1alignmentonly.sort.bam
do
echo $b
sleep 0.5
python ~/scripts/run_hapCaller.py --threads 30 -b $b --outDir ../vcf/ -r ../skylark_genome/skylark.fa --runName HC --gatkDir /home/shm45/programs/GATK/current/ --javaRam 12g -ploidy 2 --heterozygosity 0.001 --outputMode EMIT_ALL_CONFIDENT_SITES
done


#### Genotype GVCFs

#make lists of samples
samples=$(awk -F "\t" '$2=="RAS"' ../lark78.data.txt | cut -f 1)
prefix=raso26.HC

samples=$(awk -F "\t" '$2=="ESN"' ../lark78.data.txt | cut -f 1)
prefix=skyN13.HC


vcf_list_with_dashV=$(for s in $samples; do (echo "-V $s.sort.HC.g.vcf.gz "); done)

# run GenotypeGVCFs
java -Xmx12g -Djava.io.tmpdir=/local -jar /home/shm45/programs/GATK/current/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-nt 30 -R ../skylark_genome/skylark.fa --includeNonVariantSites \
 $vcf_list_with_dashV \
-o $prefix.vcf.gz


#### filter

# minimal filtering at first - remove sites without at least one called genotype

for prefix in raso26 skyN13
do
bcftools view -U $prefix.HC.vcf.gz | bgzip > $prefix.HC.clean.vcf.gz
rm $prefix.HC.vcf.gz.tbi
mv $prefix.HC.clean.vcf.gz $prefix.HC.vcf.gz
done



# Normal recombination scaffolds only

scafs=$(cat ../NR_scaffolds_250Kb_Sep2018.txt)

for prefix in raso26 skyN13
do
echo $prefix
tabix -p vcf $prefix.HC.vcf.gz
tabix -h $prefix.HC.vcf.gz $scafs | bgzip > $prefix.HC.NRscafs.vcf.gz
done


# highly filtered for sites with reads in all samples and only Normal Rec. scaffolds

for prefix in raso26 skyN13
do
echo $prefix
bcftools view -e 'FORMAT/DP < 1' $prefix.HC.NRscafs.vcf.gz | bgzip > $prefix.HC.NRscafs.DP1ALL.vcf.gz
done


