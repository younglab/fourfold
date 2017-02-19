#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use File::Copy;

sub makeerror {
  print join(" ",@_) . "\n";
  exit 1;
}

die "Arguments: <base dir> <run all?> <data template file> <output dir> <template file>" unless scalar(@ARGV) >= 3;

my ($basedir,$runall,$samplefile,$outputdir,$templatefile) = @ARGV;

die "Cannot find base directory $basedir" unless -e $basedir;
die "Cannot find data file $samplefile" unless -e $samplefile;
die "Cannot find output directory $outputdir" unless -e $outputdir;
die "Cannot find template file $templatefile" unless -e $templatefile;

my $database = ReadData($templatefile);
my $sdatabase = ReadData($samplefile);

### Step 1

print "Step 1: Reading basic options and sample groups... ";

my $sheet = $database->[1];
my $nr = ${$sheet}{"maxrow"};

my %basicoptions;

for( my $i = 0; $i < $nr; $i++ ) {
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($key,$value) = @arr;
  
  next if $key =~ /^#/; ## skip commented lines
  next if $key =~ /^$/;
  
  $basicoptions{$key} = $value;
}

$sheet = $database->[2];
$nr = ${$sheet}{"maxrow"};

my %samplegroups;

for( my $i = 0; $i < $nr; $i++ ) {
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($groupid,$value) = @arr;
  
  next if $groupid =~ /^#/; ## skip commented lines
  next if $groupid =~ /^$/;
  
  $samplegroups{$groupid} = $value;
}

unless(scalar(keys(%samplegroups))>0) {
  print "error!\nNo sample groups detected!\n";
  exit 1;
}

print "completed\n";

### Step 2: read processing


print "Step 2: Processing 4C-seq reads... ";
unless( !$runall && -e "bootstrap") { ## skip if already present 
  unless($outputdir eq '.') {
    mkdir($outputdir);
    copy($samplefile,$outputdir) or die "Cannot copy 4C template file to new output directory!";
    chdir($outputdir);
  }
  
  my $args = "";
  
  $args .= "-n " . $basicoptions{"Bowtie -n Option"} . " " if defined($basicoptions{"Bowtie -n Option"});
  $args .= "-k " . $basicoptions{"Bowtie -k Option"} . " " if defined($basicoptions{"Bowtie -k Option"});
  $args .= "-m " . $basicoptions{"Bowtie -m Option"} . " " if defined($basicoptions{"Bowtie -m Option"});
  $args .= "-r " . $basicoptions{"Number of Bootstrap Iterations"} . " " if defined($basicoptions{"Number of Bootstrap Iterations"});
  $args .= "--lsf-queue=" . $basicoptions{"LSF Queue"} . " " if defined($basicoptions{"LSF Queue"});

  my $output = `$basedir/4c-read-processing.sh $args $samplefile 2>&1`;

  makeerror "failed to trim 4c-seq reads properly (arguments $args $samplefile): $output" unless $? == 0;
  print "finished\n";
} else {
  print "skipping (reads already processed)\n";
}

#### Step 3

print "Step 3: Normalizing sample groups... ";

my %normtypes;

$sheet = $database->[3];
$nr = ${$sheet}{"maxrow"};

my $tot = 0;
my $skipped = 0;

for( my $i = 0; $i < $nr; $i++ ) {
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($groupid,$normtype) = @arr;
  
  next if $groupid =~ /^#/; ## skip commented lines
  next if $groupid =~ /^$/;
  
  $tot++;
  
  my $outputdir = "$groupid-$normtype-norm";
  
  push @{$normtypes{$groupid}}, $normtype;

  if( -e $outputdir && !$runall ) {
    $skipped++;
    next;
  }
  
  makeerror "Unknown sample group $groupid!" unless( defined($samplegroups{$groupid}));
  
  my $samples = $samplegroups{$groupid};
  
  my $output = `$basedir/4c-normalize-samples.sh $samplefile $normtype $outputdir $samples`;
  
  makeerror "Failed to normalize group $groupid\nError messages: $output" unless $? == 0;
}

print "done ($skipped/$tot, " . sprintf("%.2f%%",$skipped/$tot*100) . ", skipped)\n";

####

print "Step 4: Smoothing samples... ";

$sheet = $database->[4];
$nr = ${$sheet}{"maxrow"};

my %smoothtypes;


$tot = 0;
$skipped = 0;

for( my $i = 0; $i < $nr; $i++ ) {
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($groupid,$smoothtype,$windowsize,$stepsize) = @arr;
  
  next if $groupid =~ /^#/;
  next if $groupid =~ /^$/;
  

  makeerror "unknown group type $groupid\n" unless defined($normtypes{$groupid});
  
  my $samples = $samplegroups{$groupid};
  
  for my $normtype (@{$normtypes{$groupid}}) {
    $tot++;
    
    my $ndir = "$groupid-$normtype-norm";
    

    makeerror "missing directory\nCannot find directory for associated normalization $normtype and $groupid\n" unless -e $ndir;

    my $outputdir = "$groupid-$normtype-$smoothtype-$windowsize-$stepsize-smoothing";
    
    push @{$smoothtypes{$groupid}}, $outputdir;

    if( !$runall && -e $outputdir ) {
      $skipped++;
      next;
    }
    

    my $output = `$basedir/4c-smooth-profiles.sh --inputdir=$ndir --mode=$smoothtype $samplefile $outputdir $windowsize $stepsize $samples`;
    
    makeerror "error in smoothing\nSee error messages: $output\n" unless $? == 0; 
  }
}

print "done ($skipped/$tot, " . sprintf("%.2f%%",$skipped/$tot*100) . ", skipped)\n";


####


print "Step 5: Drawing plots... ";

$sheet = $database->[5];
$nr = ${$sheet}{"maxrow"};

$tot = 0;
$skipped = 0;

for( my $i = 0; $i < $nr; $i++ ) {
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($groupid,$region,$shading,$ylow,$yhigh,$vpos,$efile) = @arr;
  
  next if $groupid =~ /^#/;
  next if $groupid =~ /^$/;
  
  makeerror "unknown group id $groupid\n" unless defined($smoothtypes{$groupid});
  
  my $samples = $samplegroups{$groupid};
  
  for my $sdir (@{$smoothtypes{$groupid}}) {
    $tot++;
    
    makeerror "missing directory\nCannot find directory for smoothed data $sdir\n" unless -e $sdir;
    
    my $outputdir = "$sdir-row$i-plots";
    
    #    push @{$smoothtypes{$groupid}}, $outputdir;

    #if( !$runall && -e $outputdir ) {
    #  $skipped++;
    #  next;
    #}
    
    my $extraopt = "";
    
    $extraopt .= "--add-enhancers=$efile " if (defined($efile) && $efile ne "NA");
    $extraopt .= "--add-vertical-lines=$vpos " if (defined($vpos) && $vpos ne "NA");
    
    my $output = `$basedir/4c-plots.sh --inputdir=$sdir --shading=$shading --ylim-low=$ylow --ylim-high=$yhigh $extraopt $samplefile $region $outputdir $samples`;
    
    makeerror "error in generating single profile plots\nSee error messages: $output\n" unless $? == 0;
    
    $output = `$basedir/4c-multiple-sample-plots.sh --inputdir=$sdir --ylim-low=$ylow --ylim-high=$yhigh $extraopt $samplefile all $region $outputdir`;
    
    makeerror "error in generating multi profile plots\nSee error messages: $output\n" unless $? == 0; 
  }
}

print "done\n";
