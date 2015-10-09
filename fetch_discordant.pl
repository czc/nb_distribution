#!/usr/bin/perl -w
#
use strict;
use warnings;

my $file = shift or die "Usage: $0 <sam> <original_bam>\n";
my $bam = shift or die "Usage: $0 <sam> <original_bam>\n";

$| = 1;

open IN, "samtools view $bam|" or die $!;

my @ins = ();
while (<IN>) {
	last if scalar @ins == 10000;
	next if /^@/;
	my @e = split;
	next unless $e[1] & 66;
	next if $e[5] =~ /S/;
	next if $e[4] < 20;
	next if $e[8] <= 100 or $e[8] > 1000;
	push @ins, $e[8];
}

my @ins1 = sort { $a <=> $b } @ins;
my $mean = &sum(@ins1)/scalar(@ins1);
my $sqsum = 0;
foreach (@ins1) {
	$sqsum += ($_-$mean)**2;
}
my $sd = sqrt($sqsum/scalar(@ins1));
close IN;


my $upper = $mean + 3*$sd;
my $lower = $mean - 3*$sd>0?$mean-3*$sd:0;
print STDERR "mean= ", $mean, " sd=", $sd, " upper=$upper", " lower=$lower", "\n";
open IN, $file or die $!;
my @pair = ();
while (<IN>) {
	if (/^@/) {
		print;
		next;
	}
	my @e = split;
	next if ($e[1] & 256);
	push @pair, $_;
	if (@pair == 2) {
		my $line1 = $pair[0];
		my @e = split /\s+/, $line1;
		my $line2 = $pair[1];
		my @e2 = split /\s+/, $line2;
		if ($e[8] > $upper or $e[8] < -$upper or $e[5] =~ /S/ or $e2[5] =~ /S/) { # insert size abnormal
			print $line1, $line2;
			@pair = ();
			next;
		}
		#unless (($e[1]==99 and $e2[1]==147) or ($e[1]==83 and $e2[1]==163) or ($e[1]==97 and $e2[1]==145) or ($e[1]==81 and $e2[1]==161) or $e[9] =~ /CGCTCTTCC/ or $e2[9]=~/CGCTCTTCC/ or $e[9]=~/AGATCGGAAG/ or $e2[9]=~/AGATCGGAAG/ or $e[9]=~/CTGTCTCTTAT/ or $e2[9]=~/CTGTCTCTTAT/ or $e[9]=~/AGATGTGTATAA/ or $e2[9]=~/AGATGTGTATAA/) {
		unless (($e[1]==99 and $e2[1]==147) or ($e[1]==83 and $e2[1]==163)) {
			print $line1, $line2;
			@pair = ();
		}
		@pair=();
	}
}

close IN;

1;
sub median {
	my @vals = @_;
	my $len = @vals;
	if ($len % 2) {
		return $vals[int($len/2)];
	} else {
		return ($vals[$len/2] + $vals[$len/2-1])/2;
	}
}

sub Q1 {
	my @vals = @_;
	my $len = @vals;
	if ($len % 2) {
		return &median(@vals[0..(int($len/2)-1)]);
	} else {
		return &median(@vals[0..($len/2-2)]);
	}
}

sub Q3 {
	my @vals = @_;
	my $len = @vals;
	if ($len % 2) {
		return &median(@vals[(int($len/2)+1)..($len-1)]);
	} else {
		return &median(@vals[($len/2+1)]);
	}
}
sub sum {
	my @vals = @_;
	my $ret = 0;
	foreach (@vals) {
		$ret += $_;
	}
	return $ret;
}
=head
	if ($e[8] > $upper or $e[8] < -$upper or ($e[8] > -$lower and $e[8] < $lower)) { # insert size abnormal
		print $_;
	} else {
		$_ = <IN>;
		my $line2 = $_;
		my @e2 = split;
		if ($e[2] eq "*" or $e2[2] eq "*") { # one end unmappable
			print $line1, $line2;
			#} elsif ((($e[1] & 16) and !($e2[1] & 16)) or ((!($e[1] & 16)) and ($e2[1] & 16)) or (($e2[1] & 16) and !($e[1] & 16) and $e2[3] < $e[3])) { # abnormal orientation
		} elsif (!($e[1])) { # abnormal orientation
			print $line1, $line2;
		}
	}
=cut 
