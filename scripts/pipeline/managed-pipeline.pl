#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use File::Copy;

die "Arguments: <base dir> <run all?> <template file>" unless scalar(@ARGV) >= 3;

my ($basedir,$runall,$templatefile) = @ARGV;

die "Cannot find base directory $basedir" unless -e $basedir;
die "Cannot find template file $templatefile" unless -e $templatefile;

my $database = ReadData($templatefile);

### Step 1: read processing

my $sheet = $database->[1];
my $nr = ${$sheet}{"maxrow"};

my %basicoptions;

for( my $i = 0; $i < nr; $i++ ) {
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($key,$value) = @arr;
  
  next if $key =~ /^#/; ## skip commented lines
  
  $basicoptions{$key} = $value;
}

die "Need to have template file information!" unless defined($basicoptions{"Template File"});
die "Need to specify the output directory" unless defined($basicoptions{"Output Directory"});

my $samplefile = $basicoptions{"Template File"};
my $outputdir = $basicoptions{"Output Directory"};

die "Cannot find 4C template file $samplefile!" unless -e $samplefile;

print "Step 1: Processing 4C-seq reads... ";
unless( !$runall && -e "bootstrap") { ## skip if already present 
  unless($outputdir eq '.') {
    mkdir($outputdir);
    copy($samplefile,$outputdir) or die "Cannot copy 4C template file to new output directory!";
    chdir($outputdir);
  }

  my $output = `$basedir/4c-read-processing.sh $samplefile`;

  die "Failed to trim 4c-seq reads properly: $output" unless $? == 0;
  print "finished\n";
} else {
  print "skipping\n";
}

### Step 2

print "Step 2: Reading sample groups... ";

$sheet = $database->[1];
$nr = ${$sheet}{"maxrow"};

my %samplegroups;

for( my $i = 0; $i < nr; $i++ ) {
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($groupid,$value) = @arr;
  
  next if $key =~ /^#/; ## skip commented lines
  
  $samplegroups{$groupid} = $value;
}

unless(scalar(keys(%samplegroups))>0) {
  print "error!\nNo sample groups detected!\n";
  exit 1;
}

print "done\n";

#### Step 3

print "Step 3: Normalizing sample groups... ";


