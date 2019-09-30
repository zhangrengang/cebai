BAM=$1
REF=$2
export HOME=/share/home/nature
if [ ! -s $REF.fai ];then
	samtools faidx $REF &
fi

samtools depth $BAM -a -a | gzip -c > mapping.depth.gz
#	from xopen import xopen as open
python ~/src/depth_bins.py mapping.depth.gz &
samtools stats $BAM > reads.bam.stats && \
sh ~/src/mapping_rate.sh reads.bam.stats > reads.bam.stats.simple &

python ~/src/depth_plot.py mapping.depth.gz mapping.depth &
if [ "`grep chr $REF.fai -i`" ]; then
	python ~/src/depth_bins_plot.py mapping.depth.gz &
fi

python ~/src/depth_plot_busco.py mapping.depth.gz mapping.depth.busco run_ref.fa/full*tsv &
python ~/src/Busco_redundans.py run_ref.fa mapping.depth.gz > busco.depth &

python ~/src/depth_contig.py mapping.depth.gz > mapping.depth.ctgs &
python /share/home/nature/src/depth0_stat.py mapping.depth.gz > mapping_depth0.depth &

(
python ~/src/window_gc.py $REF > genome.fa.gc 
python ~/src/window_depth.py mapping.depth.gz $REF.fai > genome.fa.depth 

cat genome.fa.depth | cut -f 4 | sed '1i depth' | paste genome.fa.gc /dev/stdin > gc_depth.data && rm genome.fa.depth genome.fa.gc
Rscript ~/src/gc_depthPlot.R gc_depth.data gc_depth.data.png
) &

samtools mpileup --output-tags DP,AD -gf $REF $BAM |bcftools call -cv > samtools.vcf && \
python ~/src/vcf_het.py samtools.vcf > samtools.vcf.het &
wait;
[ "`grep chr $REF.fai -i`" ] && python /share/home/nature/src/het_bins_plot.py samtools.vcf
