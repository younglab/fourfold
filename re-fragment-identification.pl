#!/usr/bin/perl

use strict;
use Spreadsheet::Read;


### this all assumes the REs are symmetric!


if(scalar(@ARGV) < 2) {
  die "arguments <XSLX> <whole genome FA>";
}

my ($sampletable,$fasta) = @ARGV;

die "Cannot find file $sampletable!" unless -e $sampletable;

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my %enzymepairs;
my %enzymecuts;

for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  next if $arr[0] =~ /^#/;
  
  my ($samplename,undef,undef,undef,undef,undef,undef,$re1,$re2,$re1seq,$re2seq) = @arr;
  
  my $pair = "$re1-$re2";
  
  if(defined($enzymepairs{$pair})) {
    my @seqs = @{$enzymepairs{$pair}};
    die "$pair sequence changed for $samplename!" unless ($seqs[0] eq $re1seq && $seqs[1] eq $re2seq);
  }
  else {
    $enzymepairs{$pair} = [$re1seq,$re2seq];
  }
}


open(F,"<",$fasta) or die "Cannot read $fasta: $!";

my $chr = "";
my $seq = "";
my $lastpos = 1;



while(<F>) {
  chomp;
  
  if(/^>(\w+)/) {
    my $nextchr = $1;
    
    for my $seqpair (keys(%enzymepairs)) {
      my ($r1,$r2) = @{$enzymepairs{$seqpair}};
      
      my $m1 = qr/$r1/;
      my $m2 = qr/$r2/;
      
      unless(defined($enzymecuts{$seqpair})) {
        $enzymecuts{$seqpair} = [];
      }
      
      $lastpos = 1;
      while( $seq =~ /$m1/ig ) {
        my ($start,$end) = ($-[0],$+[0]);
        push @{${$enzymecuts{$seqpair}}[0]},"$chr\t$start\t$end\n";
        $lastpos = $-[0];
      }
    
      $lastpos = 1;
      while( $seq =~ /$m2/ig ) {
        my ($start,$end) = ($-[0],$+[0]);
        push @{${$enzymecuts{$seqpair}}[1]},"$chr\t$start\t$end\n";
        $lastpos = $-[0];
      }
    

    }   
    
    $chr = $nextchr;
    $seq ="";
    $lastpos = 0;

    next;
  }
  
  $seq .= $_;
}

close(F);

for my $pair (keys(%enzymecuts)) {
  mkdir($pair);
  
  my ($re1,$re2) = split /-/, $pair;
  
  open(A,">","$pair/$re1.txt") or die "Cannot write to $pair/$re1.txt: $!";
  my @arr = @{${$enzymecuts{$pair}}[0]};
  print A for(@arr);
  close(A);
  
  open(B,">","$pair/$re2.txt") or die "Cannot write to $pair/$re2.txt: $!";
  @arr = @{${$enzymecuts{$pair}}[1]};
  print B for(@arr);
  close(B);
}


