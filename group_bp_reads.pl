#!/usr/bin/perl -w
# Auther: Zechen Chong
$| = 1;
use strict;
use warnings;

my $kmerfile = shift or die "Usage: $0 <kmer.stat> <germline_novo_kmer_read1.fq> <germline_novo_kmer_read2.fq> \n";

my $rd1file = shift or die "Usage: $0 <kmer.stat> <germline_novo_kmer_read1.fq> <germline_novo_kmer_read2.fq> \n";
my $rd2file = shift or die "Usage: $0 <kmer.stat> <germline_novo_kmer_read1.fq> <germline_novo_kmer_read2.fq> \n";

my ($start, $end);

my $id = 0;
my %kmer2id = ();
my %id2pair = ();
$start = time;
printf STDERR "begin kmer2id ...\n";
open IN, $kmerfile or die $!;
my $kmer;
while (<IN>) {
	next unless /SOMATIC/;
	my @e = split /\s+/,$_;
	next if defined $kmer2id{$e[0]};
	$kmer2id{$e[0]} = $id;
	$id++;
	$kmer = $e[0];
}
my $len = length $kmer;
my $kmer_num = $id;
$end = time;
printf STDERR "kmer2id takes %d seconds\n", $end - $start;
close IN;

my @ids = (); 
my @sz = ();
for (my $i = 0; $i < $id; $i++) {
	push @ids, $i;
	push @sz, 1;
}
printf STDERR "begin id2pair...\n";
$start = time;
open IN1, $rd1file or die $!;
open IN2, $rd2file or die $!;

my @relations = ();
my ($seq1, $seq2);
while (<IN1>) {
	if ($. % 4 == 2) {
		chomp;
		$seq1 = $_;
		{
			<IN2>;
			$_ = <IN2>;
			chomp;
			$seq2 = $_;
			$id2pair{$id} = $seq1."_".$seq2;
			push @ids, $id;
			push @sz, 1;
			for (my $i = 0; $i <= length($seq1) - $len; $i++) {
				my $kmer = substr($seq1, $i, $len);
				if (exists $kmer2id{$kmer}) {
					#print STDERR ($kmer, " ", $kmer2id{$kmer}, " ", $id, "\n");
					&union($kmer2id{$kmer}, $id);
				} elsif (exists $kmer2id{&rev_com($kmer)}) {
					#print STDERR ($kmer, " ", $kmer2id{&rev_com($kmer)}, " ", $id, "\n");
					&union($kmer2id{&rev_com($kmer)}, $id);
				}
			}
			for (my $i = 0; $i <= length($seq2) - $len; $i++) {
				my $kmer = substr($seq2, $i, $len);
				if (exists $kmer2id{$kmer}) {
					#print STDERR ($kmer, " ", $kmer2id{$kmer}, " ", $id, "\n");
					&union($kmer2id{$kmer}, $id);
				} elsif (exists $kmer2id{&rev_com($kmer)}) {
					#print STDERR ($kmer, " ", $kmer2id{&rev_com($kmer)}, " ", $id, "\n");
					&union($kmer2id{&rev_com($kmer)}, $id);
				}
			}
			$id++;
			<IN2>;
			<IN2>;
		}
	}
}

close IN1;
close IN2;
$end = time;
printf STDERR "id2pair takes %d seconds\n", $end - $start;

#printf STDERR "begin clustering ...\n";

#printf STDERR "clustering takes %d seconds\n", $end - $start;

printf STDERR "begin sorting ids...\n";
$start = time;
my $line = 0;
my @idmap = ();
foreach my $i (@ids) {
	push @idmap, { newid => $i, oldid => $line };
	$line ++;
}

@idmap = sort {$a->{newid} <=> $b->{newid}} @idmap;
$end = time;
printf STDERR "sorting ids takes %d seconds\n", $end - $start;
printf STDERR "begin output results...\n";

$start = time;
my $group = -1;
my $pre = -1;

my @holder = ();

for my $i (@idmap) {
	if ($i->{newid} != $pre) {
		if (scalar @holder >= 3) {
			$group ++;
			&print_pair($group, \@holder);
		}
		@holder = ();
		$pre = $i->{newid};
		if ($i->{oldid} >= $kmer_num) {
			my $pair = $id2pair{$i->{oldid}};
			push @holder, $pair;
		} 
	} else {
		if ($i->{oldid} >= $kmer_num) {
			my $pair = $id2pair{$i->{oldid}};
			push @holder, $pair;
		} 
	}
}
		if (scalar @holder >= 3) {
			$group ++;
			&print_pair($group, \@holder);
		}
$end = time;
printf STDERR "Outputting results takes %d seconds\n", $end - $start;
printf STDERR "Finished\n";

1;

sub print_pair {
	my $gid = shift;
	my $group = shift;
	
	foreach my $pair (@$group) {
		my ($seq1, $seq2) = split /_/, $pair;
		print join("\t", ($gid, $seq1, $seq2)), "\n";
	}
}

sub rev_com {
	my $seq = shift;

	$seq =~ tr/ACGT/TGCA/;
	return reverse $seq;
}
sub union {
	my $p = shift;
	my $q = shift;
	my ($i, $j);
	for ($i = $p; $i != $ids[$i]; $i = $ids[$i]) {
		$ids[$i] = $ids[$ids[$i]];
	}
	for ($j = $q; $j != $ids[$j]; $j = $ids[$j]) {
		$ids[$j] = $ids[$ids[$j]];
	}

	return if ($i == $j);

	if ($sz[$j] + $sz[$i] < 1000) {
		if ($sz[$i] < $sz[$j]) {
			$ids[$i] = $j;
			$sz[$j] += $sz[$i];
		} else {
			$ids[$j] = $i;
			$sz[$i] += $sz[$j];
		}
	}
}
