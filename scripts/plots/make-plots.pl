#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use FourCOpts::Utils qw(issamplein);


die "Arguments: <template file> <organism database> <basedir> <genomic coordinates> <shading> <input dir> <output dir> <files 1> [files 2...]" unless scalar(@ARGV)>=8;

my ($sampletable,$organismdatabase,$basedir,$genomecoord, $shading,$inputdir,$outputdir,@files) = @ARGV;

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

my %samplegroups;

for( my $i = 0; $i < $nr; $i++ ) { ## row 1 (index 0) is the header line

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$celltype,$condition,$rep,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  next if $name =~ /^#/;
  
  next unless issamplein($name,\@files);
  
  my $skey = "$celltype-$condition";
  
  my $signalfile = "$inputdir/$name.filtered.rpm.txt";
  my $bootstrapfile = "";
  
  die "Cannot find signal file $signalfile" unless -e $signalfile;
  
  my $fullbootstrap = "$inputdir/$name.filtered.rpm.bootstrap.txt";
  my $statsbootstrap = "$inputdir/$name.filtered.rpm.bootstrap.stats.txt.gz";
  
  my $statsfile = 0;
  
  if( -e $statsbootstrap ) {
    $statsfile = 1;
    $bootstrapfile = $statsbootstrap;
  } elsif( -e $fullbootstrap ) {
    $bootstrapfile = $fullbootstrap;
  } else {
    die "Cannot find bootstrap files, either $statsbootstrap or $fullbootstrap";
  }
  
  my $output = `Rscript $basedir/plot-4c-signal.r $signalfile $bootstrapfile $statsfile $shading $genomecoord $outputdir/$name.pdf $outputdir/$name.png 2>&1`;
  
  die "Failed to generate output for $name, messages: $output" unless $? == 0;
  
  push @{$samplegroups{$skey}}, [$signalfile,$bootstrapfile,$statsfile];
}

for my $skey (keys(%samplegroups)) {
  my @fset = @{$samplegroups{$skey}};
  
  my $statsfile = $fset[0]->[2]; ## assume that we are not mixing and matching stats files here
  
  my @samples;
  
  push @samples, join(" ",@{$_}[0..1]) for(@fset);
  
  my $files = join(" ",@samples);
  my $pdfoutput = "$outputdir/$skey.pdf";
  my $pngoutput = "$outputdir/$skey.png";
  
  my $output = `Rscript $basedir/plot-joint-4c-signal.r $statsfile $shading $genomecoord $pdfoutput $pngoutput $files 2>&1`;
  
  die "Failed to generate output for $skey, messages: $output" unless $? == 0;
  
}

