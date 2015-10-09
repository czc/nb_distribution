#!/usr/bin/perl -w
#
use strict;
use warnings;

while (<>) {
	my @e = split;
	my $len = abs $1 if /SVLEN=(\S+)/;
	if ($len < 1000 and $len > 100 and $e[18] + $e[28] > 2 and $e[21] < 2 and $e[31] < 2 and $e[13] > 3 and $e[5] > 10) { # small svs
		print;
	}
	#elsif ($e[35]>=2 and $e[36]/$e[35] < 0.1 and $e[37]>=2 and $e[38]/$e[37] < 0.1 or ($e[35]<=1 and $e[37]<=1 and $e[13]>2 and $e[23]<2 and $e[33]<2 and $e[18]+$e[28]>2)) { # this may be a better filter for real data
	elsif ($e[36] <= 1 and $e[38] <= 1 and $e[23] <= 1 and $e[33] <= 1) {
		if (!/TRA/ and $e[5]>=5) {
			if (!/DUP/) {
		   		if ($len > 100000) {
					print if $e[5] > 50 or ($e[35]>5 and $e[37]>5);
				} else {
					print;
				}
			} else {
				print;
			}
		}
	} elsif (/TRA/) {
		if ($e[35]>=3 and $e[36]/$e[35] < 0.1 and $e[37]>=3 and $e[38]/$e[37] < 0.1) {
			print if $e[5] >= 30;
		} elsif (($e[35]>=3 and $e[36]/$e[35] < 0.2) ) {
			#if ($e[5] >= 30 and $e[13]>=3) {
				my $end = $e[1]+1;
				$e[7] =~ s/END=(\d+);/END=$end;/;
				$e[7] =~ s/CHR2=(\S+?);/CHR2=$e[0];/;
				$e[7] =~ s/TRA/INS/;
				$e[4] =~ s/TRA/INS/;
				print join("\t", @e), "\n";
				#}
		} elsif (($e[37]>=3 and $e[38]/$e[37] < 0.2) ) {
			#if ($e[5] >= 30 and $e[13]>=3) {
				my $chr2 = $1 if /CHR2=(\S+?);/;
				my $end2 = $1 if /END=(\d+);/;
				$e[0] = $chr2;
				$e[1] = $end2-1;
				$e[7] =~ s/TRA/INS/;
				$e[4] =~ s/TRA/INS/;
				print join("\t", @e), "\n";
				#}
		}
	}
}
