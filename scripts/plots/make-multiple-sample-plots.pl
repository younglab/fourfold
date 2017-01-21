#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use FourCOpts::Utils qw(issamplein);


die "Arguments: <template file> <multiple sample template file> <organism database> <basedir> <genomic coordinates> <shading> <input dir> <output dir> <ylim low> <ylim high> <enhancer file> <promoter file> <vertlines>" unless scalar(@ARGV)>=11;

my ($sampletable,$multisampletable,$organismdatabase,$basedir,$genomecoord,$shading,$inputdir,$outputdir,$ylimlow,$ylimhigh,$enhancerfile,$promoterfile,$vertlines) = @ARGV;

## transfer into adjustable parameters someday
#my ($linecolor,$shadingcolor,$transparencyperc) = ("black","red",50);

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $multisampletable!" unless -e $multisampletable;
die "Cannot find $organismdatabase" unless -e $organismdatabase;
die "Cannot find $inputdir" unless -e $inputdir;
die "Cannot find $outputdir" unless -e $outputdir;
die "Cannot find $enhancerfile" unless ($enhancerfile eq "NA" || -e $enhancerfile);
die "Cannot find $promoterfile" unless ($promoterfile eq "NA" || -e $promoterfile);


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
  
  my $skey = "$celltype-$condition";
  next if defined($samplegroups{$skey});

  my $statsfile = 0;
  
  my $ssf = "$inputdir/$skey.filtered.rpm.txt";
  my $sbf = "$inputdir/$skey.filtered.rpm.bootstrap.txt";
  my $altsbf = "$inputdir/$skey.filtered.rpm.bootstrap.stats.txt.gz";
  if( -e $altsbf ) {
    $statsfile = 1;
    $sbf = $altsbf;
  }
  
  #die "Cannot find $ssf and $sbf for group plots for sample $skey" unless( -e $ssf && -e $sbf );
  next unless( -e $ssf && -e $sbf );

    
  $samplegroups{$skey} = [$ssf,$sbf,$statsfile];
}

print "Making plots...\n";

my $database = ReadData($multisampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

if( !defined($sheet) || !defined($nr) ) {
  print "ERROR: Cannot seem to read $sampletable as an Excel file\n";
  exit 1;
}

for( my $i = 0; $i < $nr; $i++ ) {
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($st1,$c1,$st2,$c2,$lc1,$sc1,$at1,$lc2,$sc2,$at2) = @arr;
  
  next if $st1 =~ /^#/;
  
  my $skey1 = "$st1-$c1";
  my $skey2 = "$st2-$c2";
  my $jointkey = "$skey1-$skey2";
  
  die "Cannot find cell type $st1 and condition $c1 in the sample template file!" unless defined($samplegroups{$skey1});
  die "Cannot find cell type $st2 and condition $c2 in the sample template file!" unless defined($samplegroups{$skey2});
  
  my ($signalfile1, $bootstrapfile1, $statsfile) = @{$samplegroups{$skey1}};
  my ($signalfile2, $bootstrapfile2) = @{$samplegroups{$skey2}};

  
  my $pdfoutput = "$outputdir/$jointkey.pdf";
  my $pngoutput = "$outputdir/$jointkey.png";
  
  print "\tPlotting $jointkey...\n";
  
  my $fileargs = "$signalfile1 $bootstrapfile1 $lc1 $sc1 $at1 $signalfile2 $bootstrapfile2 $lc2 $sc2 $at2";
  
  my $output = `Rscript $basedir/plot-4c-signal.r $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $fileargs 2>&1`;

  die "Failed to generate output for $jointkey (arguments $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $fileargs), messages: $output" unless $? == 0;
}

