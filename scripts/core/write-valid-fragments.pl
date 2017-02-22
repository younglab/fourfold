#!/usr/bin/perl

use strict;

if(scalar(@ARGV)<2) {
  die "write-valid-fragments.pl <fragment file> <output BED file>";
}

my ($infile,$outfile) = @ARGV;

die "Cannot find fragment file $infile!" unless -e $infile;

my @frags;

open(F,"<",$infile) or die "Cannot read: $!";
open(O,">",$outfile) or die "Cannot write: $!";

while(<F>) {
  chomp;
  
  my ($chr,$s,$e,$l,$r) = split /\t/;
  
  unless($l eq "NA") {
    my $beg = $s-$l;
    print O "$chr\t$beg\t$s\n";
  }
  
  unless($r eq "NA") {
    my $end = $e+$r;
    print O "$chr\t$e\t$end\n";
  }
}

close(F);
close(O);
