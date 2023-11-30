#!/usr/bin/perl
#
# written by Kai Wang (kw2110@columbia.edu)
#
use warnings;
use strict;

my $work_dir = '/users/kw2110/project/_ARACNE';

print "Testing ARACNE parallel computation ... \n";

if ( !(-e $work_dir."/aracne") ) {
    print "Please compile the wokring ARACNE code first!\n";
	exit(1);
}

print "==> Reconstructing full matrix ... \n";
`$work_dir/aracne -i $work_dir/test/arraydata10x336.exp -k 0.16 -t 0.05 -e 0.1`;

print "Test I: reconstruction in parallel with a posterior processing.\n";
print "==> Reconstructing full matrix in parallel ";
for my $i (0 .. 9) {
    `$work_dir/aracne -i $work_dir/test/arraydata10x336.exp -h $i -k 0.16`;
    print ".";
}
print "\n";

my $tempout = $work_dir."/test/temp.adj";
open(OUT, "> $tempout");

print "==> Collecting results ";
for my $j (0 .. 9) {
	my $adjfile = $work_dir."/test/arraydata10x336_h".$j."_k0.16.adj";
	if ( -e $adjfile ) {
    	open(ADJ, "< $adjfile");
		while (<ADJ>) {
				if (!/^>/ && $_) {
					   print OUT $_;
			   	}
		}
        close(ADJ);
        print ".";

		# sleep 0.2;
        `rm $adjfile`;
	} else {
    	print "\nMissing file $adjfile ... \n";
        exit(0);
    }
}
close(OUT);
print "\n";

print "==> Applying thresholding & DPI a posterior ...";
`$work_dir/aracne -i $work_dir/test/arraydata10x336.exp -j $tempout -t 0.05 -e 0.1 -o $work_dir/test/arraydata10x336_k0.16_t0.05_e0.1_parallel_separate.adj`;
my $rslt1 = `diff -I "^>" $work_dir/test/arraydata10x336_k0.16_t0.05_e0.1.adj $work_dir/test/arraydata10x336_k0.16_t0.05_e0.1_parallel_separate.adj`;

if (!$rslt1) {
    print " Success!\n";
}
else {
	print " Fail!\n";
	exit(1);
}

print "Test II: subnetwork reconstruction.\n";
print "==> Reconstructing subnetwork ...";
`printf "G2\nG5\nG8\n" > $work_dir/test/list.dat`;
`$work_dir/aracne -i $work_dir/test/arraydata10x336.exp -s $work_dir/test/list.dat -k 0.16 -t 0.05 -e 0.1 -o $work_dir/test/arraydata10x336_k0.16_subnet.adj`;
`grep "^G2" $work_dir/test/arraydata10x336_k0.16_t0.05_e0.1.adj >> $work_dir/test/foo.adj`;
`grep "^G5" $work_dir/test/arraydata10x336_k0.16_t0.05_e0.1.adj >> $work_dir/test/foo.adj`;
`grep "^G8" $work_dir/test/arraydata10x336_k0.16_t0.05_e0.1.adj >> $work_dir/test/foo.adj`;
my $rslt2 = `diff -I "^>" $work_dir/test/foo.adj $work_dir/test/arraydata10x336_k0.16_subnet.adj`;

if (!$rslt2) {
    print " Success!\n";
}
else {
	print " Fail!\n";
	exit(1);
}
`rm $work_dir/test/list.dat`;
`rm $work_dir/test/*.adj`;

