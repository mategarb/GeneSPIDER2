#!/usr/bin/perl -w

######################################################
#                                                    #
# Program to change ADJ file format to 3 columns     #
#                                                    #
# Usage: perl adjix2id.pl <adj file>                 #
# Output: <adj>_3col.txt                             #
#                                                    #
######################################################

use Fcntl;
($adjfile) = @ARGV;
$adjfile =~ m/^(\S+)(.adj)/; 
$outfile = $1."_3col.txt";

sysopen(OUT, $outfile, O_WRONLY|O_TRUNC|O_CREAT, 0666) || die $!;

open(FH, $adjfile) or die "Error, could not open $adjfile\n";
while (<FH>) {
  chomp;
  @col = split("\t",$_);
  $numtar = $#col/2;
  foreach $i (1..$numtar) {
    print OUT $col[0],"\t",$col[(($i-1)*2)+1],"\t",$col[$i*2],"\n";
  }
}
close FH;

close OUT;

