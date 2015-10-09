#!/usr/bin/perl -w
#
#Auther: Zechen Chong
#
use strict;
use warnings;
use English;

my $vcf = shift or die "Usage: $0 <vcf> <tumor.bam> <normal.bam>\n";
my $tbam = shift or die "Usage: $0 <vcf> <tumor.bam> <normal.bam>\n";
my $nbam = shift or die "Usage: $0 <vcf> <tumor.bam> <normal.bam>\n";

open IN, $vcf or die $!;
while (<IN>) {
	print if /^#/;
	next if /^#/;
	chomp;
	my $line = $_;
	my @e = split /\s+/, $_;
	my $fh;
	my $beg = $e[1]-500;
	my $fin = $e[1]+500;
	my $chr1 = $e[0];
	my $end = 0;
	if ($line =~ /END=(\d+);/) {
		$end = $1;
	}
	my $chr2;
	if ($line =~ /CHR2=(\S+?);/) {
		$chr2 = $1;
	}
	my $beg2 = $end-500;
	my $fin2 = $end+500;
	my $len;
	if (/SVLEN=(\S+)/) {
		$len = $1>=0?$1:-$1;
	}
	system("samtools", "view", $tbam, $chr1.":".$beg."-".$fin, "-o", "$PID.regionsam2.tmp");
#	`samtools view $tbam $e[0]:$e[1]-$e[1] > regionsam2.tmp`;
	open $fh , "$PID.regionsam2.tmp" or die $!;
	my @pos1_tum = count_sp($fh, $e[1]);
	seek $fh, 0, 0;
	my $dr1_tum = count_dr(\$line, $fh, $len);
	close $fh;
	system("samtools", "view", $nbam, $chr1.":".$beg."-".$fin, "-o", "$PID.regionsam2.tmp");
#	`samtools view $nbam $e[0]:$e[1]-$e[1] > regionsam2.tmp`;
	open $fh , "$PID.regionsam2.tmp" or die $!;
	my @pos1_nor = count_sp($fh, $e[1]);
	seek $fh, 0, 0;
	my $dr1_nor = count_dr(\$line, $fh, $len);
	close $fh;
	system("samtools", "view", $tbam, $chr2.":".$beg2."-".$fin2, "-o", "$PID.regionsam2.tmp");
#	`samtools view $tbam $e[0]:$e[1]-$e[1] > regionsam2.tmp`;
	open $fh , "$PID.regionsam2.tmp" or die $!;
	my @pos2_tum = count_sp($fh, $end);
	seek $fh, 0, 0;
	my $dr2_tum = count_dr(\$line, $fh, $len);
	close $fh;
	system("samtools", "view", $nbam, $chr2.":".$beg2."-".$fin2, "-o", "$PID.regionsam2.tmp");
#	`samtools view $nbam $e[0]:$e[1]-$e[1] > regionsam2.tmp`;
	open $fh , "$PID.regionsam2.tmp" or die $!;
	my @pos2_nor = count_sp($fh, $end);
	seek $fh, 0, 0;
	my $dr2_nor = count_dr(\$line, $fh, $len);
	close $fh;
	print join("\t", ($line, @pos1_tum, @pos1_nor, @pos2_tum, @pos2_nor, $dr1_tum, $dr1_nor, $dr2_tum, $dr2_nor)), "\n";
	`rm $PID.regionsam2.tmp`;
}
close IN;

1;

sub count_dr {
	my $l = shift;
	my $fd = shift;
	my $len = shift;
	my $line = $$l;
	my ($del, $dup) = (0, 0);
	my %hash = ();
	while (<$fd>) {
		my @e = split /\s+/, $_;
		if ($len == 0) {
			if ($e[6] ne '=') {
				$hash{$e[0]} ++;
			}
		} else {
			my $len2 = $e[8]>0?$e[8]:-$e[8];
			if (abs($len-$len2) < 500) {
				$hash{$e[0]} ++;
				if ($line =~ /DUP/ or $line =~ /DEL/) {
					if (((($e[1]&144)==144) and $e[8]>0) or ((($e[1]&80)==80) and $e[8]>0) or ((($e[1]&96)==96) and $e[8]<0) or ((($e[1]&160)==160) and $e[8]<0)) { #duplication
						$dup ++;
					}
					if (((($e[1]&96)==96) and $e[8]>0) or ((($e[1]&160)==160) and $e[8]>0) or ((($e[1]&144)==144) and $e[8]<0) or ((($e[1]&80)==80) and $e[8]<0)) { # deletion
						$del ++;
					}
				}
			}
		}
		#last if scalar keys %hash > 100;
		#if ($e[8] > 1000 or $e[8] < -1000) {
		#next if $e[4] < 29;
		#if (($e[8] > 600 and $e[8] < 1000000) or ($e[8] < -600 and $e[8] > -1000000) or ($e[8] > -30 and $e[8] < 30)) {
		#	$hash{$e[0]} ++;
		#}
	}
	if ($line =~ /DUP/ or $line =~ /DEL/) { ##TODO
		if ($dup > $del+1) {
			$$l =~ s/DEL/DUP/g;
		} elsif ($del > $dup+1) {
			$$l =~ s/DUP/DEL/g;
		}
		return scalar keys %hash;
	}
	return scalar keys %hash;
}

sub count_sp {
	my ($fd, $target) = @_;
	my $cnt = 0;
	my $qual = 0;
	my $cnt2 = 0;
	my $qual2 = 0;
	my $total = 0;
	while (<$fd>) {
		chomp;
		my @e = split /\s+/, $_;
		$total++ if $e[3] < $target and $e[3]+100>$target;
		next if $e[1] == 4 or $e[5] !~ /[SH]/;
		my $bp;
		if ($e[5] =~ /(\d+)[SH].+M(\d+)[SH]/) {
			next if ($1>5 and $2>5);
		}
		if ($e[5] =~ /(\d+)M.+?(\d+)M/) {
			next if ($1>10 and $2>10);
		}
		
		my ($m1, $s1) = (0, 0);
		while ($e[5] =~ /(\d+)[SH]/g) {
				$s1 = $1 if $1 > $s1;
		}
		while ($e[5] =~ /(\d+)M/g) {
				$m1 = $1 if $1 > $m1;
		}
		if ($e[5] =~ /$m1[M].*?$s1[SH]/) {
			$bp = $e[3]+$m1-1;
		} else {
			$bp = $e[3];
		}
		

		my @scores = split //, $e[10];
		my $num = 0;
		foreach (@scores) {
			$num ++ if $_ le '#';
		}
		#	print $e[10], "\t", $num, "\t", $s1, "\n";
		next if (abs($num - $s1)<4);

		if ($bp >= $target-1 and $bp <= $target+1) {
			$cnt ++;
			$qual += $e[4];
			if ($e[4] >= 29) {
				$cnt2 ++;
				$qual2 += $e[4];
			}
		}
		#$cnt ++ if $bp == $target;
		#print $e[2], "\t", $bp, "\t", $e[4], "\t", $e[5], "\n";
	}
	$qual = $qual/$cnt if $cnt > 0;
	$qual2 = $qual2/$cnt2 if $cnt2 > 0;

	return ($total, $cnt, $qual, $cnt2, $qual2);
}
