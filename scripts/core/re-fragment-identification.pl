#!/usr/bin/perl

use strict;
use Spreadsheet::Read;


sub findresites {
  my ($re1,$re2,$r1,$r2,$chr,$seq,$ecuts) = @_;
  
  return if $seq eq "";
 
  my $m1 = qr/$r1/i;
  my $m2 = qr/$r2/i;

  while( $seq =~ /$m1/g ) {
    my ($start,$end) = ($-[0],$+[0]);
    push @{$ecuts->[0]}, [$chr,$start,$end,$re1];
  }
    
  while( $seq =~ /$m2/g ) {
    my ($start,$end) = ($-[0],$+[0]);
    push @{$ecuts->[1]},[$chr,$start,$end,$re2];
  }
  #}   
}

### this all assumes the REs are symmetric!


if(scalar(@ARGV) < 3) {
  die "arguments <XSLX> <organism database> <min distance>";
}

my ($sampletable,$organismdatabase,$mindistance) = @ARGV;

die "Cannot find file $sampletable!" unless -e $sampletable;
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



my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my %enzymepairs;
#my %enzymecuts;
  

for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  
  my($name,undef,undef,undef,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  next if $name =~ /^#/;

  my $pair = "$e1-$e2-$organism";
  
  die "Cannot find $organism in database!" unless defined($organisms{lc $organism});
  
  my @od = @{$organisms{lc $organism}};
  
  if(defined($enzymepairs{$pair})) {
    my @seqs = @{$enzymepairs{$pair}};
    die "$pair sequence changed for $name!" unless ($seqs[0] eq $e1s && $seqs[1] eq $e2s);
  }
  else {
    $enzymepairs{$pair} = [$e1s,$e2s,$od[1]];
  }
}


for my $kpair (keys(%enzymepairs)) {
  my ($re1,$re2,$organism) = split /-/, $kpair;
  my ($re1s,$re2s,$fasta) = @{$enzymepairs{$kpair}};
  
  my @enzymecuts = ([],[]);
  
  open(F,"<",$fasta) or die "Cannot read $fasta: $!";

  my $chr = "";
  my $seq = "";

  while(<F>) {
    chomp;
  
    if(/^>(\w+)/) {
      my $nextchr = $1;
    
      findresites($re1,$re2,$re1s,$re2s,$chr,$seq,\@enzymecuts);
    
      $chr = $nextchr;
      $seq ="";

      next;
    }
  
    $seq .= $_;
  }

  close(F);

  findresites($re1,$re2,$re1s,$re2s,$chr,$seq,\@enzymecuts);
  
  undef $seq; ### clean up some memory

  mkdir($kpair);
  
  open(A,">","$kpair/$re1.bed") or die "Cannot write to $kpair/$re1.bed: $!";
  print A "track name=\"$re1\" description=\"$re1\"\n";

  print A join("\t",@{$_}) . "\n" for(@{$enzymecuts[0]});
  close(A);
  
  open(B,">","$kpair/$re2.bed") or die "Cannot write to $kpair/$re2.bed: $!";
  print B "track name=\"$re2\" description=\"$re2\"\n";
  print B join("\t",@{$_}) . "\n" for(@{$enzymecuts[1]});
  close(B);

  `tail -n +2 $kpair/$re1.bed > $kpair/.re1.bed`;
  `tail -n +2 $kpair/$re2.bed > $kpair/.re2.bed`;
  `cat $kpair/.re1.bed $kpair/.re2.bed | sort -k1,1 -k2,2n > $kpair/.tmp.bed`;
  `rm $kpair/.re1.bed $kpair/.re2.bed`;
  
  open(T,"<","$kpair/.tmp.bed") or die "Cannot read $kpair/.tmp.bed: $!";

  
  my %fragmentset;
  my @cuts;
  my @sites;
  my @past;
  
  push @past, ["NA",-1,-1];
  
  @cuts = <T>;
  
  close(T);

  chomp @cuts;
  push @sites, [split( /\t/,$_)] for(@cuts);
  
  undef @cuts;

  for(@sites) {
    my @cur = @{$_};
    my @l = @{$past[$#past]};

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
  
  
  undef @past;
  
  open(F,">","$kpair/.fragments.txt") or die "Cannot write to $kpair/.fragments.txt: $!";
  
  for my $site (keys(%fragmentset)) {
    my @dists = @{$fragmentset{$site}};
  
    print F "$site\t$dists[0]\t$dists[1]\n";
  }
  
  close(F);
  
  `sort -k1,1 -k2,2n $kpair/.fragments.txt > $kpair/fragments.txt`;
  unlink("$kpair/.fragments.txt");
  unlink("$kpair/.tmp.bed");

  
  undef @enzymecuts;
}


