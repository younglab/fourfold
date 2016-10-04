#!/usr/bin/perl

use strict;

if(scalar(@ARGV)<3) {
  die "arguments <whole genome FA> <RE 1 seq> <RE 2 seq>";
}

my ($fasta,$r1,$r2) = @ARGV;

open(F,"<",$fasta) or die "Cannot read $fasta: $!";

my $chr = "";
my $seq = "";
my $lastpos = 1;

my $m1 = qr/$r1/;
my $m2 = qr/$r2/;

while(<F>) {
  chomp;
  
  if(/^>(\w+)/) {
    my $nextchr = $1;
    
    while( $seq =~ /$m1/ig ) {
      my $end = $-[0]-1;
      print "$chr\t$lastpos\t$end\n";
      $lastpos = $-[0];
    }
    
    $lastpos = 1;
    while( $seq =~ /$m2/ig ) {
      
    }
    
    $chr = $nextchr;
    $seq ="";
    $lastpos = 1;
    next;
  }
  
  $seq .= $_;
}

close(F);
