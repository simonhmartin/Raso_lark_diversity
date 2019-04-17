
# run nucmer
MUMmer3.23/nucmer --prefix=Tgut_skylark_scaf1M ../zebrafinch/Tgut_chroms.fa ../skylark_genome/skylark.scafs1M.fa

# filter for hits longer than 5 kb
MUMmer3.23/delta-filter -l 5000  Tgut_skylark_scaf1M.delta > Tgut_skylark_scaf1M.len5k.delta


# plot with layout adjusted to arrange scaffolds
MUMmer3.23/mummerplot Tgut_skylark_scaf1M.len5k.delta -l -t png -s large \
-R ../zebrafinch/Tgut_chroms.fa \
 -Q ../skylark_genome/skylark.scafs1M.fa


# NOTE This creates a gnuplot script (out.gp). Now we modify the scale. Set the size to 3000x3000 and reduce point size (ps 0.1)

sed -e 's/tiny size 1400,1400/size 3000,3000/' -e 's/out.png/Tgut_skylark_scaf1M.len5k.ordered.png/' -e 's/ps 1/ps 0.1/g'  out.gp > out.fixed.gp

# Then run gnuplot to make the final plot
gnuplot out.fixed.gp 
