#!/usr/bin/perl

use strict;
use Spreadsheet::Read;


die "Arguments: <template file> <organism database> <basedir> <genomic coordinates> <shading> <input dir> <output dir>" unless scalar(@ARGV)>=6;

my ($sampletable,$organismdatabase,$basedir,$genomecoord, $shading,$inputdir,$outputdir) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $organismdatabase" unless -e $organismdatabase;
die "Cannot find $inputdir" unless -e $inputdir;
die "Cannot find $outputdir" unless -e $outputdir;

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

if( !defined($sheet) || !defined($nr) ) {
  print "ERROR: Cannot seem to read $sampletable as an Excel file\n";
  exit 1;
}


for( my $i = 0; $i < $nr; $i++ ) { ## row 1 (index 0) is the header line

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$celltype,$condition,$rep,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  next if $name =~ /^#/;
  
  my $output = `Rscript $basedir/plot-4c-signal.r $inputdir/$name.filtered.rpm.txt $inputdir/$name.filtered.rpm.bootstrap.txt $shading $genomecoord $outputdir/$name.pdf $outputdir/$name.png 2>&1`;
  
  die "Failed to generate output for $name, messages: $output" unless $? == 0;
}

