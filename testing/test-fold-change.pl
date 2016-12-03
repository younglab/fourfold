#!/usr/bin/perl

use strict;

die "arguments: <BED> <BAMs...>" unless scalar(@ARGV)>=2;

my ($bed,@bams) = @ARGV;

for(@bams) {
  die unless -e $_;
}

open(B,"<",$bed) or die "Cannot read $bed";

my $l = <B>;
chomp $l;

my ($chr,$s,$e) = split /\t/, $l;

close(B);

open(O,">","tad-temp.bed") or die;

for(my $i=0;$i<10000;$i++) {
  my $r = rand(1);
  my $p = $s+int(($e-$s-10000)*$r);
  
  print O "$chr\t$p\t" . ($p+10000) . "\n";
}

close(O);

for(@bams) {
  `bedtools intersect -c -wa -a tad-temp.bed -b $_ > tmp.bed`;
  `mv tmp.bed tad-temp.bed`;
}


