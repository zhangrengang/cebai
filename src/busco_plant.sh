
#$ -q midmem+.q
export PATH="/share/home/app/bin/hmmer-3.1b2/bin:$PATH"

LINEAGE=/share/home/yuanqyk/database/BUSCO_db/embryophyta_odb9


SEQS=$1
opts=($@)
opts=${opts[@]:1:${#opts[*]}} # for compatible
#echo $opts
#echo $@
bsname=$(basename $SEQS)
if [ -d run_$bsname ]; then
	echo "run_$bsname exists and skip" && exit 0
fi
ref0=`realpath $SEQS`
ref0base=$(basename $ref0)
ref0dir=$(dirname $ref0)
if [ -d $ref0dir/run_$ref0base ]; then
	ln $ref0dir/run_$ref0base run_$bsname -s && echo "use the real dir" && exit 0
fi

mode=geno
opts0=$(echo $opts | cut -b1)
if [ "x$opts" != "x" ] && [ "$opts0" != "-" ]; then
	mode=
fi

work_dir=$(pwd)
tmpdir=/io/node/$(hostname)/busco_$$
mkdir -p $tmpdir && cp $SEQS $tmpdir && cd $tmpdir
BUSCO.py -i $bsname -o $bsname -l $LINEAGE -c 8 -z --tmp /io/tmp/$(hostname)/ -m $mode $opts \
&& mv $tmpdir/run_$bsname $work_dir \
&& xrm $tmpdir
