#!/usr/bin/perl -w
#
use strict;
use warnings;
use English;

my $pre = -1;
my @reads = ();

my $file = shift or die "$0 <bp_reads.txt>\n";
open RD, $file or die $!;
open OUT, ">$file.localreads.fa" or die $!;
my $old_fh = select(OUT);
$| = 1;
select($old_fh);
while (<RD>) {
	chomp;
	my @e = split /\s+/, $_;
	if ($pre != $e[0]) {
		if (@reads > 0) {
			my $i = 0;
			seek OUT, 0, 0;
			truncate OUT, 0;
			my %seen = ();
			@reads = grep { ! $seen{$_}++ } @reads;
			foreach my $r (@reads) {
				print OUT ">read$i\n";
				print OUT $r, "\n";
				$i++;
			}
		my $consensus = "";
		my $name = "";
		#	system("cap3 localreads.fa -i 30 -j 31 -o 18 -s 300 > localreads.fa.log;");
		system("SSAKE -m 16 -f $file.localreads.fa -w 1 -z 50 -o 1 -b $file.ssake.asm");
		open OUT2, '>>', "$file.ssake.asm.out" or die $!;
		my $old_fh2 = select(OUT2);
		$| = 1;
		select($old_fh2);
		open IN, "$file.ssake.asm.contigs" or die $!;
		while (<IN>) {
			if (/^>(\S+)/) {
				if ($name ne "") {
					print OUT2 ">", $name, "\n", $consensus, "\n";
				}
				$name = "E".$pre."|".$1;
			} else {
				chomp;
				$consensus = $_;
			}
		}

		print OUT2 ">", $name, "\n", $consensus, "\n";
		close IN;
		close OUT2;
#		system("cat localreads.fa.cap.contigs >> cap3.asm.out");

		$pre = $e[0];
		@reads = ();
	}
		push @reads, ($e[1], $e[2]);
		$pre = $e[0];
	} else {
		push @reads, ($e[1], $e[2]);
	}
}

		if (@reads > 0) {
			my $i = 0;
			seek OUT, 0, 0;
			truncate OUT, 0;
			my %seen = ();
			@reads = grep { ! $seen{$_}++ } @reads;
			foreach my $r (@reads) {
				print OUT ">read$i\n";
				print OUT $r, "\n";
				$i++;
			}
		} 
		system("SSAKE -m 16 -f $file.localreads.fa -w 1 -z 50 -o 1 -b $file.ssake.asm");
		#system("cap3 localreads.fa -i 30 -j 31 -o 18 -s 300 > localreads.fa.log;");
		open OUT2, '>>', "$file.ssake.asm.out" or die $!;
		open IN, "$file.ssake.asm.contigs" or die $!;
		my $consensus = "";
		my $name = "";
		while (<IN>) {
			if (/^>(\S+)/) {
				if ($name ne "") {
					print OUT2 ">", $name, "\n", $consensus, "\n";
				}
				$name = "E".$pre."|".$1;
			} else {
				chomp;
				$consensus = $_;
			}
		}
		print OUT2 ">", $name, "\n", $consensus, "\n";
		close IN;
		close OUT2;
		#system("cat localreads.fa.cap.contigs >> cap3.asm.out");
close OUT;
close RD;
