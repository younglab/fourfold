#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use FourCOpts::Utils qw(issamplein convertcoordinatestring);


die "Arguments: <template file> <organism database> <basedir> <genomic coordinates> <shading> <input dir> <output dir> <group only> <ylim low> <ylim high> <enhancer file> <promoter file> <verticle line coordinates> <CI alpha> <files 1> [files 2...]" unless scalar(@ARGV)>=13;

my ($sampletable,$organismdatabase,$basedir,$genomecoord, $shading,$inputdir,$outputdir,$grouponly,$ylimlow,$ylimhigh,$enhancerfile,$promoterfile,$vertlines,$cialpha,@files) = @ARGV;

## transfer into adjustable parameters someday
my ($linecolor,$shadingcolor,$transparencyperc) = ("black","red",50);

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $organismdatabase" unless -e $organismdatabase;
die "Cannot find $inputdir" unless -e $inputdir;
die "Cannot find $outputdir" unless -e $outputdir;
die "Cannot find $enhancerfile" unless ($enhancerfile eq "NA" || -e $enhancerfile);
die "Cannot find $promoterfile" unless ($promoterfile eq "NA" || -e $promoterfile);

die "CI is not an appropriate fractional number in (0,1]!" unless ($cialpha eq "1" || $cialpha =~ m/^0?[.]\d+$/);

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

print "Making individual sample plots...\n" unless $grouponly;

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
  
  unless( defined($samplegroups{$skey})) { ## assume this hasn't changed from one iteration to the next if already there
    my $ssf = "$inputdir/$skey.filtered.rpm.txt";
    my $sbf = "$inputdir/$skey.filtered.rpm.bootstrap.txt";
    $sbf = "$inputdir/$skey.filtered.rpm.bootstrap.stats.txt.gz" if $statsfile;
  
    die "Cannot find $ssf and $sbf for group plots for sample $skey" unless( -e $ssf && -e $sbf );
    
    $samplegroups{$skey} = [$ssf,$sbf,$statsfile];
  }
  
  unless( $grouponly ) { ## skip plots if only grouponly is turned on
    print "\tPlotting $name...\n";
    
    my $pdfoutput = "$outputdir/$name.pdf";
    my $pngoutput = "$outputdir/$name.png";
    
    my $output = `Rscript $basedir/plot-4c-signal.r $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $cialpha $signalfile $bootstrapfile $linecolor $shadingcolor $transparencyperc 2>&1`;
  
    die "Failed to generate output for $name, messages: $output" unless $? == 0;
  }
}

if(scalar(keys(%samplegroups))==0) {
  die "No samples matched the list of sample IDs passed!";
}

print "Plotting combined sample group plots...\n";

for my $skey (keys(%samplegroups)) {
  my ($signalfile, $bootstrapfile, $statsfile) = @{$samplegroups{$skey}};
  
  my $pdfoutput = "$outputdir/$skey.pdf";
  my $pngoutput = "$outputdir/$skey.png";
  
  print "\tPlotting $skey...\n";
  my $output = `Rscript $basedir/plot-4c-signal.r $genomecoord $shading $statsfile $ylimlow $ylimhigh $enhancerfile $promoterfile $vertlines $pdfoutput $pngoutput $cialpha $signalfile $bootstrapfile $linecolor $shadingcolor $transparencyperc  2>&1`;

  die "Failed to generate output for $skey, messages: $output" unless $? == 0;
}

