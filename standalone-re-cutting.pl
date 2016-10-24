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


if(scalar(@ARGV) < 7) {
  die "arguments <RE1> <RE1 seq> <RE2> <RE2 seq> <organism> <organism database> <min distance>";
}

my ($e1,$e1s,$e2,$e2s,$organism,$organismdatabase,$mindistance) = @ARGV;

die "Cannot fine organismal database $organismdatabase!" unless -e $organismdatabase;

my %organisms;

open(D,"<","$organismdatabase") or die "Cannot read $organismdatabase: $!";

while(<D>) {
  chomp;
  
  my ($id,$bowtie,$fasta) = split /\t/;
  
  next if $id =~ /^$/;
  
  die "In organism database, cannot find bowtie index for $id!" unless -e "$bowtie.1.ebwt";
  die "In organism database, cannot fine FASTA file for $id!" unless -e $fasta;
  
  for(split(/,/,$id)) {
    $organisms{lc $_} = [$bowtie,$fasta];
  }
}
close(D);

my %enzymepairs;
my %enzymecuts;
  

my $pair = "$e1-$e2-$organism";

die "Cannot find organism $organism in organism database!" unless defined($organisms{lc $organism});

my @od = @{$organisms{lc $organism}};

$enzymepairs{$pair} = [$e1s,$e2s,$od[1]];


for my $kpair (keys(%enzymepairs)) {

  my $fasta = $enzymepairs{$kpair}->[2];
  
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

}
