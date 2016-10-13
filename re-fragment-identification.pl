#!/usr/bin/perl

use strict;
use Spreadsheet::Read;


sub findresites {
  my ($epairs,$chr,$seq,$ecuts) = @_;
  
  for my $seqpair (keys(%{$epairs})) {
    my ($re1,$re2) = split /-/, $seqpair;
    my ($r1,$r2) = @{$epairs->{$seqpair}};
    
    my $m1 = qr/$r1/i;
    my $m2 = qr/$r2/i;
      
    unless(defined($ecuts->{$seqpair})) {
      $ecuts->{$seqpair} = [];
    }
      
    while( $seq =~ /$m1/g ) {
      my ($start,$end) = ($-[0],$+[0]);
      push @{$ecuts->{$seqpair}->[0]}, [$chr,$start,$end,$re1];
    }
    
    while( $seq =~ /$m2/g ) {
      my ($start,$end) = ($-[0],$+[0]);
      push @{$ecuts->{$seqpair}->[1]},[$chr,$start,$end,$re2];
    }
  }   
}

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

while(<F>) {
  chomp;
  
  if(/^>(\w+)/) {
    my $nextchr = $1;
    
    unless($seq eq "") {
      findresites(\%enzymepairs,$chr,$seq,\%enzymecuts);
    }
    
    $chr = $nextchr;
    $seq ="";

    next;
  }
  
  $seq .= $_;
}

close(F);

unless($seq eq "") {
  findresites(\%enzymepairs,$chr,$seq,\%enzymecuts);
}

for my $pair (keys(%enzymecuts)) {
  mkdir($pair);
  
  my ($re1,$re2) = split /-/, $pair;
  
  open(A,">","$pair/$re1.bed") or die "Cannot write to $pair/$re1.bed: $!";
  print A "track name=\"$re1\" description=\"$re1\"\n";
  my @first = @{${$enzymecuts{$pair}}[0]};
  print A join("\t",@{$_}) . "\n" for(@first);
  close(A);
  
  open(B,">","$pair/$re2.bed") or die "Cannot write to $pair/$re2.bed: $!";
  print B "track name=\"$re2\" description=\"$re2\"\n";
  my @second = @{${$enzymecuts{$pair}}[1]};
  print B join("\t",@{$_}) . "\n" for(@second);
  close(B);

  `tail -n +2 $pair/$re1.bed > $pair/.re1.bed`;
  `tail -n +2 $pair/$re2.bed > $pair/.re2.bed`;
  `cat $pair/.re1.bed $pair/.re2.bed | sort -k1,1 -k2,2n > $pair/.tmp.bed`;
  `rm $pair/.re1.bed $pair/.re2.bed`;
  
  open(T,"<","$pair/.tmp.bed") or die "Cannot read $pair/.tmp.bed: $!";

  
  my %fragmentset;
  my @cuts;
  my @sites;
  my @past;
  
  @cuts = <T>;
  
  close(T);

  chomp @cuts;
  push @sites, [split( /\t/,$_)] for(@cuts);

  for(@sites) {
    my @cur = @{$_};
    my @l = ("NA",-1,-1);
    @l = @{$past[$#past]} if scalar(@past)>0;

    if($cur[3] eq $re2) {
      if($l[3] eq $re1 && $cur[0] eq $l[0] && $l[1] > 0 && ($cur[1]-$l[1]>=$mindistance) ) {
        my $site = join("\t",@l[0..2]);

        $fragmentset{$site}->[1] = $cur[1]-$l[1];
      }
    } else {
      my $site = join("\t",@cur[0..2]);

      $fragmentset{$site} = ["NA","NA"];
      
      if($l[3] eq $re2 && $cur[0] eq $l[0] && $l[1] > 0 && ($cur[1]-$l[1]>=$mindistance)) {
        $fragmentset{$site}->[0] = $cur[1]-$l[1];
      }
    }
    
    push @past, $_;
  }
  
  open(F,">","$pair/.fragments.txt") or die "Cannot write to $pair/.fragments.txt: $!";
  
  for my $site (keys(%fragmentset)) {
    my @dists = @{$fragmentset{$site}};
    
    print F "$site\t$dists[0]\t$dists[1]\n";
  }
  
  close(F);
  
  `sort -k1,1 -k2,2n $pair/.fragments.txt > $pair/fragments.txt`;
}


