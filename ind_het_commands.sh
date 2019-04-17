# make Geno files required as input
for prefix in raso26 skyN13 skyWR9 skyER10 skyCH4 ori11 crest5
do
echo $prefix
bcftools filter -e 'FORMAT/DP < 5' --set-GTs "." -O u $prefix.HC.vcf.gz | bcftools view -U |
python ~/Research/genomics_general/VCF_processing/parseVCF.py --skipIndels | bgzip > $prefix.HC.DP5.geno.gz &
done

# individual heterozygosity
for pref in raso26 skyN13 skyWR9 skyER10 skyCH4 ori11 crest5
do
python popgenWindows.py -g ../vcf/$pref.HC.DP5.geno.gz -f phased --analysis indHet popDist -o $pref.HC.DP5.het.w250m50s50.csv -w 250000 -m 50 -s 50000 --roundTo 5 -T 10
done


# convert scaffold positions to geneome positions
for pref in raso26 skyN13 skyWR9 skyER10 skyCH4 ori11 crest5
do
python transferScafPos.py -i $pref.HC.DP5.het.w250m50s50.csv \
-o $pref.HC.DP5.het.w250m50s50.chrom.csv --endCol 3 --header --sep "," --failFile fails.txt \
--transfersFile Tgut_skylark_scaf1M.len5k.order_manual.transfers.tsv
done

