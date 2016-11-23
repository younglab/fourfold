#!/usr/bin/perl

use strict;
use Spreadsheet::Read;


sub in {
  my ($test,$ref) = @_;
  
  my $ret = 0;
  
  if(scalar(@{$ref})==1) {
    $ret = 1 if ($ref->[0] eq "all" || $test =~ m/$ref->[0]/);
  } else {
  
    for(@{$ref}) {
      if($test eq $_) {
        $ret = 1;
        last;
      }
    }
  }
  
  return $ret;
}

#print "DEBUG: @ARGV\n";
if(scalar(@ARGV) < 5) {
  die "process-quantile-norm.pl <sample table> <basedir> <organism database> <output dir> <samples...>";
}

my ($sampletable,$basedir,$organismdatabase,$outputdir,@samples) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $outputdir!" unless -e $outputdir;

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

my @norm;
my %samplegroups;

for( my $i = 0; $i < $nr; $i++ ) { 
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$celltype,$condition,$rep,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  next if $name =~ /^#/;
  
  if(in($name,\@samples)) {
    push @norm, [$name,"bootstrap/$name.filtered.rpm.txt"];
    
    my $samplekey = "$celltype-$condition";
    $samplegroups{$samplekey} = [] unless defined($samplegroups{$samplekey});
    
    push @{$samplegroups{$samplekey}}, "$outputdir/$name.filtered.rpm.txt";
  }
}

die "not enough samples" if scalar(@norm) < 2;

my $exe = "";
$exe .= join(" ",@{$_}) . " " for(@norm);

#print "DEBUG: $exe\n";

my $output = `Rscript $basedir/quantile-normalization.r $outputdir $exe`;

die "Failed to normalize: $output" unless $? == 0;

open(N,"<","$outputdir/quantile-normalized-samples.txt") or die "Cannot read normalized samples: $!";

<N>; ## skip header

my @output;

for(@norm) {
  my $ref = [$_->[0],[]]; ### place sample name and start array to build output wig file
  push @{$ref->[1]}, "track type=wiggle_0 name=\"$_->[1]QuantNorm\" description=\"$_->[1] Quantile Normalized\"";
  push @output, $ref;
}

my $lastchr = "";

while(<N>) {
  chomp;
  
  my ($chr,$pos,@vals) = split /\t/;
  
  if($chr ne $lastchr) {
    my $s = "variableStep chrom=$chr";
    
    push @{$_->[1]},$s for(@output);
    $lastchr = $chr;
  }
  
  for( my $i = 0; $i <= $#vals; $i++ ) {
    push @{$output[$i]->[1]}, "$pos\t$vals[$i]";
  }
}

close(N);

for(@output) {
  my $outputfile ="$outputdir/$$_[0].filtered.rpm.wig";
  open(O,">",$outputfile) or die "Cannot write to $outputfile: $!";
  print O "$_\n" for(@{$_->[1]});
  close(O);
}

unlink("$outputdir/quantile-normalized-samples.txt");

for my $skey (keys(%samplegroups)) { 
  print "DEBUG: $skey\n";
  my @files = @{$samplegroups{$skey}};
  
  next unless scalar(@files)>=2;
  
  my $ofile = "$outputdir/$skey.filtered.rpm.txt";
  my $wfile = "$outputdir/$skey.filtered.rpm.wig";

  my $f = join(" ",@files);

  
  my $output = `Rscript $basedir/combine-profiles.r $ofile $f`;
  
  open(T,"<","$ofile") or die "Cannot read $ofile: $!";
  open(W,">",$wfile) or die "Cannot write to $wfile: $!";
  
  print W "track type=wiggle_0 name=\"$skey-QuantNorm\" description=\"$skey-QuantNorm\"\n";
  
  my $lastchr = "";
  
  while(<T>) {
    chomp;
    
    my ($chr,$pos,$s) = split /\t/;
    
    unless($lastchr eq $chr) {
      print W "variableStep chrom=$chr\n";
      $lastchr =$chr;
    }
    
    print W "$pos\t$s\n";
  }
  
  close(T);
  close(W);
}
  
