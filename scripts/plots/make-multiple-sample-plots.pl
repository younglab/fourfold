#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use FourCOpts::Utils qw(issamplein convertcoordinatestring);

sub runfromtemplate {
  my ($multisampletable,$basedir,$genomecoord,$shading,$inputdir,$outputdir,$ylimlow,$ylimhigh,$enhancerfile,$promoterfile,$vertlines,$cialpha,$usecairopng,$sgref) = @_;
  
  my %samplegroups = %{$sgref};
  my $database = ReadData($multisampletable);
  my $sheet = $database->[1]; ## get first spreadsheet
  my $nr = ${$sheet}{"maxrow"};

  if( !defined($sheet) || !defined($nr) ) {
    print "ERROR: Cannot seem to read $multisampletable as an Excel file\n";
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
  
    my $fileargs = "$skey1 $signalfile1 $bootstrapfile1 $lc1 $sc1 $at1 $skey2 $signalfile2 $bootstrapfile2 $lc2 $sc2 $at2";
  
    my $output = `Rscript $basedir/plot-4c-signal.r $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $cialpha $usecairopng $fileargs 2>&1`;

    die "Failed to generate output for $jointkey (arguments $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $cialpha $fileargs), messages: $output" unless $? == 0;
  }

}

sub runallsets {
  my ($basedir,$genomecoord,$shading,$inputdir,$outputdir,$ylimlow,$ylimhigh,$enhancerfile,$promoterfile,$vertlines,$cialpha,$usecairopng,$sgref) = @_;
  
  my %samplegroups = %{$sgref};
  
  my @groups = keys(%samplegroups);
  my $ngroups = scalar(@groups);
  
  for( my $i = 0; $i < $ngroups; $i++ ) {
    for( my $j = $i+1; $j < $ngroups; $j++ ) {
      my ($g1,$g2) = ($groups[$i],$groups[$j]);
      
      my $jointkey = "$g1-$g2";
      
      my ($signalfile1, $bootstrapfile1, $statsfile) = @{$samplegroups{$g1}};
      my ($signalfile2, $bootstrapfile2) = @{$samplegroups{$g2}};
      
      my $pdfoutput = "$outputdir/$jointkey.pdf";
      my $pngoutput = "$outputdir/$jointkey.png";
  
      print "\tPlotting $jointkey...\n";
  
      my $fileargs = "$g1 $signalfile1 $bootstrapfile1 red red 50 $g2 $signalfile2 $bootstrapfile2 blue blue 50";
  
      my $output = `Rscript $basedir/plot-4c-signal.r $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $cialpha $usecairopng $fileargs 2>&1`;

      die "Failed to generate output for $jointkey (arguments $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $cialpha $usecairopng $fileargs), messages: $output" unless $? == 0;

    }
  }
  
  ### need to simplify this...
  
  for( my $i = 0; $i < $ngroups; $i++ ) {
    for( my $j = $i+1; $j < $ngroups; $j++ ) {
      for( my $k = $j+1; $k < $ngroups; $k++) {
        my ($g1,$g2,$g3) = ($groups[$i],$groups[$j],$groups[$k]);
      
        my $jointkey = "$g1-$g2-$g3";
      
        my ($signalfile1, $bootstrapfile1, $statsfile) = @{$samplegroups{$g1}};
        my ($signalfile2, $bootstrapfile2) = @{$samplegroups{$g2}};
        my ($signalfile3, $bootstrapfile3) = @{$samplegroups{$g3}};

      
        my $pdfoutput = "$outputdir/$jointkey.pdf";
        my $pngoutput = "$outputdir/$jointkey.png";
  
        print "\tPlotting $jointkey...\n";
  
        my $fileargs = "$g1 $signalfile1 $bootstrapfile1 red red 50 $g2 $signalfile2 $bootstrapfile2 blue blue 50 $g3 $signalfile3 $bootstrapfile3 green green 50";
  
        my $output = `Rscript $basedir/plot-4c-signal.r $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $cialpha $usecairopng $fileargs 2>&1`;

        die "Failed to generate output for $jointkey (arguments $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $cialpha $usecairopng $fileargs), messages: $output" unless $? == 0;

      }
    }
  }
}


die "Arguments: <template file> <multiple sample template file> <organism database> <basedir> <genomic coordinates> <shading> <input dir> <output dir> <ylim low> <ylim high> <enhancer file> <promoter file> <vertlines> <CI alpha> <Cairo PNG?>" unless scalar(@ARGV)>=15;

my ($sampletable,$multisampletable,$organismdatabase,$basedir,$genomecoord,$shading,$inputdir,$outputdir,$ylimlow,$ylimhigh,$enhancerfile,$promoterfile,$vertlines,$cialpha,$cairopng) = @ARGV;
my ($usecairopng) = ($cairopng eq "no" ? "false" : "true");

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $multisampletable!" unless (-e $multisampletable || $multisampletable eq "all");
die "Cannot find $organismdatabase" unless -e $organismdatabase;
die "Cannot find $inputdir" unless -e $inputdir;
die "Cannot find $outputdir" unless -e $outputdir;
die "Cannot find $enhancerfile" unless ($enhancerfile eq "NA" || -e $enhancerfile);
die "Cannot find $promoterfile" unless ($promoterfile eq "NA" || -e $promoterfile);

$genomecoord = convertcoordinatestring($genomecoord);
die "$genomecoord is an invalid genomic coordinate string!" unless $genomecoord ne "";


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

if( $multisampletable eq "all") {
  runallsets($basedir,$genomecoord,$shading,$inputdir,$outputdir,$ylimlow,$ylimhigh,$enhancerfile,$promoterfile,$vertlines,$cialpha,$usecairopng,\%samplegroups);
} else {
  runfromtemplate($multisampletable,$basedir,$genomecoord,$shading,$inputdir,$outputdir,$ylimlow,$ylimhigh,$enhancerfile,$promoterfile,$vertlines,$cialpha,$usecairopng,\%samplegroups);
} 
