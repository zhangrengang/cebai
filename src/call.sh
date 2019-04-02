#!/bin/bash
ploidy=2
fb_opts="--min-repeat-entropy 0 --min-alternate-count 2 --min-alternate-fraction 0.2 --ploidy $ploidy --use-best-n-alleles 4 --genotype-qualities "

function calling_parallel {
	ref=$1
    bam=$2
    outDir=$3
    ncpu=$4
    region_length=$5
    (fasta_generate_regions.py $ref.fai $region_length | \
        parallel -k -j $ncpu freebayes -f $ref $bam $fb_opts -r {} ) | \
         vcffirstheader | vcfstreamsort		# stdout
}
function calling_sge {
	ref=$1
    bam=$2
    outDir=$3
    ncpu=$4
    region_length=$5
    i=0
	src=$outDir/calling.sh
    for region in `fasta_generate_regions.py $ref.fai $region_length`
    do
        i=$[$i+1]
        echo "if [ \$SGE_TASK_ID -eq $i ]; then"
        echo "freebayes $fb_opts -f $ref $bam -r $region > $outDir/freebayes.$region.vcf"
        echo "fi"
    done > $src
    wait_exit `qsub -terse -cwd -j y -V -S /bin/bash -t 1-$i -tc $ncpu -o $src.out $src`
    cat $outDir/freebayes.*.vcf | vcffirstheader | vcfstreamsort	# stdout
}
function calling_wrapper {
	mode=$1
    ref=$2
    bam=$3
    outDir=$4
    ncpu=$5
    region_length=$6
    samtools faidx $ref
    samtools index $bam
    
	[ $mode = "parallel" ] && calling_parallel $ref $bam $outDir $ncpu $region_length	# stdout
	[ $mode = "sge" ]      && calling_sge	   $ref $bam $outDir $ncpu $region_length	# stdout
}

# SGE
function wait_exit {
    JID=$1
	i=0
    while :
    do
        STATS=`qacct -j $JID 2> /dev/null | grep exit_status | awk '{print $2}'`
        if [ x"$STATS" != x ]; then
            break
        fi
        QSTATS=`qstat | awk '$1=="'$JID'"{print $5}'`
        if [ x"$QSTATS" != x ]; then
            sleep 1m
        fi
        if [ x"$STATS" = x ] && [ x"$QSTATS" = x ]; then
            i=$[$i+1]
            if [ $i -gt 5 ];then
                break
            fi
        fi
        sleep 1m
    done
}

    # call variants
ref=$1
bam=$2
outDir=$3
ncpu=$4
mode=sge
region_length=2000000
#    region_length=$5
mkdir -p $outDir
calling_wrapper $mode $ref $bam $outDir $ncpu $region_length
	
