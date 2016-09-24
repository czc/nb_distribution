#!/usr/bin/perl -w
#
use strict;
use warnings;

while (<>) {
	my @e = split;
	my $len = abs $1 if /SVLEN=(\S+)/;
	if ($len < 1000 and $len > 100 and $e[18] + $e[28] > 4 and $e[21] <= 1 and $e[31] <= 1 and $e[13] > 3 and $e[5] > 10) { # small svs
		print if $e[15]>0 and $e[16]/$e[15]>0.1 or ($e[25]>0 and $e[26]/$e[25]>0.1);
	} elsif ($e[35]>=3 and $e[36]/$e[35] < 0.1 and $e[37]>=3 and $e[38]/$e[37] < 0.1 and $e[36] <=2 and $e[38]<=2 or ($e[35]<=1 and $e[37]<=1 and $e[13]>4 and $e[21]<1 and $e[31]<1 and $e[18]>=2 and $e[28]>=2)) { # this may be a better filter for real data
#	elsif ($e[36] <= 1 and $e[38] <= 1 and $e[23] <= 1 and $e[33] <= 1) {
		if (/INV/) {
			#print if $len < 20000000;
			print;
		} elsif (/TRA/) {
			print if $e[35] >= 3 or $e[37] >= 3;
		} else {
			print;
		}
	} elsif (/TRA/) {
		print if ($e[13]>2 and $e[5] > 50 and $e[21]<=1 and $e[31]<=1 and $e[36]<=1 and $e[38]<=1 or ($e[13]==2 and $e[36]<=1 and $e[38]<=1 and $e[35]>=3 and $e[37]>=3));
	}
}
