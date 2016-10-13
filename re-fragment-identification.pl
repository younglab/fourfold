#!/usr/bin/perl

use strict;
use Spreadsheet::Read;


### this all assumes the REs are symmetric!


if(scalar(@ARGV) < 3) {
  die "arguments <XSLX> <whole genome FA> <min distance>";
}

my ($sampletable,$fasta,$mindistance) = @ARGV;

die "Cannot find file $sampletable!" unless -e $sampletable;
die "Cannot fine FASTA file $fasta!" unless -e $fasta;

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my %enzymepairs;
my %enzymecuts;

for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  next if $arr[0] =~ /^#/;
  
  my ($samplename,undef,undef,undef,undef,undef,undef,undef,undef,undef,$re1,$re2,$re1seq,$re2seq) = @arr;
  
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
      my ($re1,$re2) = split /-/, $seqpair;
      my ($r1,$r2) = @{$enzymepairs{$seqpair}};
      
      my $m1 = qr/$r1/i;
      my $m2 = qr/$r2/i;
      
      unless(defined($enzymecuts{$seqpair})) {
        $enzymecuts{$seqpair} = [];
      }
      
      $lastpos = 1;
      while( $seq =~ /$m1/g ) {
        my ($start,$end) = ($-[0],$+[0]);
        push @{${$enzymecuts{$seqpair}}[0]},"$chr\t$start\t$end\t$re1\n";
        $lastpos = $-[0];
      }
    
      $lastpos = 1;
      while( $seq =~ /$m2/g ) {
        my ($start,$end) = ($-[0],$+[0]);
        push @{${$enzymecuts{$seqpair}}[1]},"$chr\t$start\t$end\t$re2\n";
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
  
  open(C,">","$pair/tmp.bed") or die "Cannot write to $pair/tmp.bed: $!";
  open(A,">","$pair/$re1.bed") or die "Cannot write to $pair/$re1.bed: $!";
  print A "track name=\"$re1\" description=\"$re1\"\n";
  my @arr = @{${$enzymecuts{$pair}}[0]};
  print A for(@arr);
  print C for(@arr);
  close(A);
  
  open(B,">","$pair/$re2.bed") or die "Cannot write to $pair/$re2.bed: $!";
  print B "track name=\"$re2\" description=\"$re2\"\n";
  @arr = @{${$enzymecuts{$pair}}[1]};
  print B for(@arr);
  print C for(@arr);
  close(B);
  close(C);
  
  open(J,">","$pair/joint-re.bed") or die "Cannot write to $pair/joint-re.bed: $!";
  print J "track name=\"$re1-$re2\" description=\"$re1-$re2\"\n";
  close(J);

  `sort -k1,1 -k2,2n $pair/tmp.bed >> $pair/joint-re.bed`;
  unlink("$pair/tmp.bed");
  
  my @fragmentset;
  my @leftsecondre = ();
  my @rightsecondre = ();
  my $lastchr = "DEADBEEF";
  my $curstart = 0;
  my $curend = 0;
  
  open(O,"<","$pair/joint-re.bed") or die "Cannot read $pair/joint-re.bed: $!";
  <O>; ## skip header line
  
  while(<O>) {
    chomp;
    
    my ($c,$s,$e,$r) = split /\t/;
    
    if( $c ne $lastchr ) {
      @leftsecondre = (); ## dump second RE digests on last chromosome
      @rightsecondre = ();
      $lastchr = $c;
    }
    
    ### if hit RE2
    if( $r eq $re2 ) {
      unshift @rightsecondre, [$c,$s,$e];
    } else {  ## otherwise we hit a RE1
      if($curstart == 0 ) { ## the first one we've seen on the chromosome
        @leftsecondre = @rightsecondre;
        @rightsecondre = ();
        
        $curstart = $s;
        $curend = $e;
      } else {
        my $cut = "$lastchr\t$curstart\t$curend";
        
        my $left = "NA";
        for(@leftsecondre) {
          my @a = @{$_};
          if(($curstart-$a[1])>=$mindistance) {
            $left = $curstart-$a[1];
            last;
          }
        }
        my $right = "NA";
        for(@rightsecondre) {
          my @a = @{$_};
          if(($a[1]-$curstart)>=$mindistance) {
            $right = $curstart-$a[1];
            last;
          }
        }
        
        @leftsecondre = reverse @rightsecondre;
        push @fragmentset, "$cut\t$left\t$right\n";
        
        $curstart = $s;
        $curend = $e;
      }
    }
  }
  close(O);
  
  open(F,">","$pair/fragments.txt") or die "Cannot write to $pair/fragments.txt: $!";
  print F for(@fragmentset);
  close(F);
}


