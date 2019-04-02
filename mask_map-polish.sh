#ref=ref2.fa
ref=$1
depth=depth.$ref.fb1/mapping.depth
ncpu=20
tmpdir=tmpdir.pilon
python /share/home/nature/src/low_coverage.py $depth > $depth.low.bed
python /share/home/nature/src/mask_bed.py $ref $depth.low.bed > $ref.mask.fa

# map
outDir=$tmpdir/$ref.mapping-2
sortedBam=$outDir.mapped.bam
/share/home/app/src/map.sh $ref.mask.fa $tmpdir/reads.list $outDir $ncpu > $sortedBam
samtools index $sortedBam
python ~/src/filter_bam3.py $sortedBam $sortedBam.filter.bam

# vall
outDir=$tmpdir/$ref.calling-2
outVcf=$outDir.vcf
/share/home/app/src/call.sh $ref.mask.fa $sortedBam.filter.bam $outDir $ncpu  > $outVcf
python ~/src/filter_freebayes.py $outVcf | awk '$4 !~ "N"'> $outVcf.filter.vcf

# merge
cat $tmpdir/$ref.calling.1.freebayes.vcf.filter.vcf $outVcf.filter.vcf |  vcffirstheader | vcfstreamsort | vcfuniq > $tmpdir/$ref.calling.merge.vcf

[ -s $ref.changes ] && mv $ref.changes $ref.1.changes
python /share/home/nature/src/vcf2fasta.py $tmpdir/$ref.calling.merge.vcf $ref > $ref.mask-polished.fa

qqsub ~/line/assembly/quality-assess/busco_plant.sh $ref.mask-polished.fa


exit
ref=ref2.fa-2.polished.fa
lastBam=$tmpdir/ref2.fa.mapping.1.mapped.bam
samtools merge $sortedBam.p.bam $sortedBam $lastBam -pc -@ $ncpu -f
samtools flagstat $sortedBam.p.bam > $sortedBam.p.bam.flagstat &
samtools depth -aa $sortedBam.p.bam > $sortedBam.p.bam.depth
python ~/src/depth_plot.py $sortedBam.p.bam.depth $sortedBam.p.bam.depth &
python ~/src/depth_contig.py $sortedBam.p.bam.depth > $sortedBam.p.bam.depth.ctgs &
wait;

