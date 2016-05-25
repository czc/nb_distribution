#!/bin/bash

if [ $# != 5 -a $# != 6 ]; then
	echo $0 \<novoBreak_exe_dir\> \<ref\> \<tumor_bam\> \<normal_bam\> \<n_CPUs:INT\> \[outputdir:-PWD\]
	exit 1
fi

nbbin=`readlink -f $1`
ref=`readlink -f $2`
tumor_bam=`readlink -f $3`
normal_bam=`readlink -f $4`
n_cpus=$5
if [ $# == 6 ]; then
	output=`readlink -f $6`
fi
novobreak=$nbbin/novoBreak
bwa=$nbbin/bwa
samtools=$nbbin/samtools

lastdir=`pwd`

if [ $# == 6 ]; then
	mkdir $output
	cd $output
fi
$novobreak -i $tumor_bam -c $normal_bam -r $ref  -o kmer.stat 
#$samtools collate somaticreads.bam somaticreads.srt

mkdir group_reads
cd group_reads
#$samtools view -h ../somaticreads.srt.bam | perl $nbbin/fetch_discordant.pl - $tumor_bam > discordant.sam
$samtools bam2fq -1 read1.fq -2 read2.fq ../somaticreads.bam
perl $nbbin/group_bp_reads.pl ../kmer.stat read1.fq read2.fq  > bp_reads.txt
cls=`tail -1 bp_reads.txt | cut -f1`
rec=`echo $cls/$n_cpus | bc`
rec=$((rec+1))
mkdir split
cd split
awk -v rec=$rec '{print > int($1/rec)".txt"}' ../bp_reads.txt
for file in *.txt
do
	echo $file
	perl $nbbin/run_ssake.pl $file > /dev/null &
done
wait
cd ..
cd ..
mkdir ssake
cd ssake
#you can split the bp_reads.txt into multiple files to run them together
#perl $nbbin/run_ssake.pl ../group_reads/bp_reads.txt > /dev/null
awk 'length($0)>1' ../group_reads/split/*.ssake.asm.out > ssake.fa
$bwa mem -t $n_cpus -M $ref ssake.fa > ssake.sam
perl $nbbin/infer_sv.pl ssake.sam > ssake.vcf
grep -v '^#' ssake.vcf | sed 's/|/\t/g' | sed 's/read//' |  awk '{if(!x[$1$2]){y[$1$2]=$14;x[$1$2]=$0}else{if($14>y[$1$2]){y[$1$2]=$14; x[$1$2]=$0}}}END{for(i in x){print x[i]}}' | sort -k1,1 -k2,2n  | perl -ne 'if(/TRA/){print}elsif(/SVLEN=(\d+)/){if($1>100){print $_}}elsif(/SVLEN=-(\d+)/){if($1>100){print}}' > ssake.pass.vcf
#you can split the ssake.pass.vcf into multiple files to run them together
num=`wc -l ssake.pass.vcf | cut -f1 -d' '`
rec=`echo $num/$n_cpus | bc`
rec=$((rec+1))
mkdir split
cd split
split -l $rec ../ssake.pass.vcf # set proper split parameters when needed
for file in x??
do
	echo $file
	perl $nbbin/infer_bp_v4.pl $file $tumor_bam $normal_bam > $file.sp.vcf &
done
wait
cd ..

##below is a naive filter, pay attention to it
#perl $nbbin/filter_sv_icgc.pl nbasm.pass.sp.vcf > ../novoBreak.pass.flt.vcf
grep '^#' ssake.vcf > header.txt	
perl $nbbin/filter_sv_icgc.pl split/*.sp.vcf | cat header.txt - > ../novoBreak.pass.flt.vcf
cd ..
cd $lastdir
