#!/usr/bin/perl

use strict;

my %reads;

sub procseq {
  my ($seq,$rl) = @_;
  
  my $slen = length($seq);
  
  for( my $i = 0; $i < ($slen-$rl); $i++ ) {
    my $s = substr($seq,$i,$rl);
    
    $reads{$s} = 1;
  }
}

if(scalar(@ARGV)<2) {
  die "Syntax: generate-mappability-track.pl <genome FASTA file> <read length> <output WIG file>";
}

my ($fasta,$readlen,$wigout) = @ARGV;

die "Cannot find FASTA file $fasta" unless -e $fasta;

die "Read length has to be at least 10" unless $readlen >= 10;

open(F,"<",$fasta) or die "Cannot read $fasta: $!";

my $seq = "";
my $curchr = "";

while(<F>) {
  chomp;
  
  if(/^>(\w+)/) {
    procseq($seq,$readlen);
    $seq =""
  } else {
    $seq .= $_;
  }
}

close(F);


procseq($seq,$readlen);

my $idx = 1;

open(T,">",".tmp.fasta") or die "Cannot write to .tmp.fasta: $!";

for(keys(%reads)) {
  print T ">read$idx\n$_\n";
}

close(T);
