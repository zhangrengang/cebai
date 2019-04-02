#!/bin/bash
ref=$1
readsList=$2
tmpDir=$3
mode=$4    # [parallel, sge]
ncpu=$5

minimap2_opts="--sr"
#bwa_opts="-x ont2d"
#bwa_opts="-k14 -A1 -B1 -O1 -E1 -L0"
bwa_opts="-k15 -A1 -B2 -O1 -E1 -L1 -r5 -U 35"

function mapping_sge {
    refPrefix=$1
    readsList=$2
    outDir=$3
    ncpu=$4
    src=$outDir/mapping.sh
    i=0
    while read LINE
    do
        i=$[$i+1]
        arr=($LINE)
        r1=${arr[0]}
        r2=${arr[1]}
        cmd=""
        echo "if [ \$SGE_TASK_ID -eq $i ]; then"
#        echo "bwa mem $refPrefix $r1 $r2 -t 1 -R '@RG\tID:SR\tSM:SR' | samtools view -F 2304 -b | \
		echo "
if [ ! -f $outDir/samtools.$i.bam.ok ]; then 
#minimap2 -a $minimap2_opts -t 3 -R '@RG\tID:SR\tSM:SR' $refPrefix $r1 $r2 | samtools view -F 2304 -b | 
bwa mem $bwa_opts $refPrefix $r1 $r2 -t 1 -R '@RG\tID:SR\tSM:SR' | samtools view -F 2304 -b | \
            samtools sort -@ 1 -T /tmp/samtools.$i.$$ > $outDir/samtools.$i.bam && \
touch $outDir/samtools.$i.bam.ok
fi
"
        echo "fi"
    done < $readsList > $src
    wait_exit `qsub -terse -cwd -j y -V -S /bin/bash -t 1-$i -tc $ncpu -o $src.out $src`
    samtools merge - $outDir/samtools.*.bam -pc -@ $ncpu		#stdout
}
function mapping_parallel {
    refPrefix=$1
    readsList=$2
    outDir=$3
    ncpu=$4
    i=0
    while read LINE
    do
        i=$[$i+1]
        arr=($LINE)
        r1=${arr[0]}
        r2=${arr[1]}
        [ $i -eq 1 ] && samtools_opts="-h"
#        bwa mem $refPrefix $r1 $r2 -t $ncpu -R '@RG\tID:SR\tSM:SR' | samtools view -F 2304 -b | \
		minimap2 -a $minimap2_opts -t 1 -R '@RG\tID:SR\tSM:SR' $refPrefix $r1 $r2 | samtools view -F 2304 -b | \
            samtools view $samtools_opts -@ $ncpu
    done < $readsList 2> $outBam.err | samtools sort -@ $ncpu -T /tmp/samtools.$$	#stdout
}
function mapping_wrapper {
	mode=$1
    ref=$2
    readsList=$3
    outDir=$4
    ncpu=$5
#    refPrefix=$ref
    refPrefix=$outDir/ref
	if [ ! -f $refPrefix.ok ]; then
	    bwa index $ref -p $refPrefix && touch $refPrefix.ok
	fi
	[ $mode = "parallel" ] && mapping_parallel $refPrefix $readsList $outDir $ncpu	# stdout
	[ $mode = "sge" ]      && mapping_sge	   $refPrefix $readsList $outDir $ncpu	# stdout
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

# chunk reads
# main pipeline
# map reads
mode=sge
ref=$1
readsList=$2
outDir=$3
ncpu=$4
mkdir -p $outDir
mapping_wrapper $mode $ref $readsList $outDir $ncpu
