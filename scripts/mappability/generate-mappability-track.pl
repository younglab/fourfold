#!/usr/bin/perl

use strict;


sub procseq {
  my ($chr,$seq,$rl) = @_;
  
  my $slen = length($seq);
  
  open(T,">>",".tmp.fasta") or die "Cannot write to .tmp.fasta! $!";
  for( my $i = 0; $i < ($slen-$rl); $i++ ) {
    my $s = substr($seq,$i,$rl);
    my $j = $i+$rl;
    
    print T ">$chr:$i-$j\n$s\n";
  }
  close(T);
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
    procseq($curchr,$seq,$readlen) unless $seq eq "";
    $seq ="";
    $curchr=$1;
    print "Reading $1...\n";
  } else {
    $seq .= $_;
  }
}

close(F);

procseq($curchr,$seq,$readlen);


