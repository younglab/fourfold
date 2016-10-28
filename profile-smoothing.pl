#!/usr/bin/perl

use strict;
use Spreadsheet::Read;

if(scalar(@ARGV)<7) {
  die "profile-smoothing.pl <sample table> <basedir> <organism databse> <output dir> <bin size> <step size> <samples...>";
}

my ($sampletable,$basedir,$organismdatabase,$outputdir,$binsize,$stepsize,@samples) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $outputdir!" unless -e $outputdir;
die "Cannot find $organismdatabase" unless -e $organismdatabase;

my %organisms;

open(D,"<","$organismdatabase") or die "Cannot read $organismdatabase: $!";

while(<D>) {
  chomp;
  
  my ($id,$bowtie,$fasta,$chromsizes) = split /\t/;
  
  next if $id =~ /^$/;
  
  die "In organism database, cannot find bowtie index for $id!" unless -e "$bowtie.1.ebwt";
  die "In organism database, cannot find FASTA file for $id!" unless -e $fasta;
  die "In organism database, cannot find chromosome sizes for $id!" unless -e $chromsizes;
  
  for(split(/,/,$id)) {
    $organisms{lc $_} = [$bowtie,$fasta,$chromsizes];
  }
}
close(D);

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my %samplegroups;

for( my $i = 0; $i < $nr; $i++ ) { 
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$celltype,$condition,$rep,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  next if $name =~ /^#/;
  
  unless( $samples[0] eq "all") {
    my $found = 0;
    
    for(@samples) {
      if(lc $_ eq lc $name) {
        $found = 1;
        last;
      }
    }
    
    next unless $found;
  }
  
  push @{$samplegroups{"$celltype:$condition"}}, ["bootstrap/$name.filtered.rpm.txt", "bootstrap/$name.filtered.rpm.bootstrap.txt"];
  
  my $chromsizes = $organisms{$organism}->[2];

  my $output = `Rscript $basedir/smooth-single-profile.r $binsize $stepsize $chromsizes bootstrap/$name.filtered.rpm.txt bootstrap/$name.filtered.rpm.bootstrap.txt $outputdir/$name.filtered.smoothed.rpm.txt`;
  
  die "Smoothing failed with an error: $output" unless( $? == 0 );
  
  ### generate a WIG file
  open(T,"<","$outputdir/$name.filtered.smoothed.rpm.txt") or die "Cannot read $outputdir/$name.filtered.smoothed.rpm.txt: $!";
  open(W,">","$outputdir/$name.filtered.smoothed.rpm.wig") or die "Cannot write $outputdir/$name.filtered.smoothed.rpm.wig: $!";
  
  print W "track type=wiggle_0 name=\"$name\" description=\"$name\"\n";
  my $lastchr = "";
  while(<T>) {
    chomp;
    
    my ($chr,$start,$end,$s,$lb,$ub) = split /\t/;
    
    unless($lastchr eq $chr) {
      print W "variableStep chrom=$chr span=$stepsize\n";
      $lastchr=$chr;
    }
    
    print W "$start\t$s\n";
  }
  
  close(W);
  close(T);
}

for my $group (keys(%samplegroups)) {
  my ($celltype,$condition) = split /:/, $group;
  
  ### code
}

