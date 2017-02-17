#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use File::Copy;

sub makeerror {
  print "$_\n";
  exit 1;
}

die "Arguments: <base dir> <run all?> <data template file> <output dir> <template file>" unless scalar(@ARGV) >= 3;

my ($basedir,$runall,$samplefile,$outputdir,$templatefile) = @ARGV;

die "Cannot find base directory $basedir" unless -e $basedir;
die "Cannot find data file $samplefile" unless -e $samplefile;
die "Cannot find output directory $outputdir" unless -e $outputdir;
die "Cannot find template file $templatefile" unless -e $templatefile;

my $database = ReadData($templatefile);

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

  my $output = `$basedir/4c-read-processing.sh $args $samplefile`;

  makeerror "failed to trim 4c-seq reads properly (arguments $args $samplefile): $output" unless $? == 0;
  print "finished\n";
} else {
  print "skipping\n";
}

#### Step 3

print "Step 3: Normalizing sample groups... ";

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
  
  if( -e $outputdir && !$runall ) {
    $skipped++;
    next;
  }
  
  makeerror "Unknown sample group $groupid!" unless( defined($samplegroups{$groupid}));
  
  my $samples = $samplegroups{$groupid};
  
  my $output = `$basedir/4c-normalize-samples.sh $templatefile $normtype $outputdir $samples`;
  
  makeerror "Failed to normalize group $groupid\nError messages: $output" unless $? == 0;
}

print "done ($skipped/$tot, " . sprintf("%.3f%%",$skipped/$tot) . ", skipped)\n";

####

print "Step 4\n";

$sheet = $database->[4];
$nr = ${$sheet}{"maxrow"};



####


print "Step 5\n";

$sheet = $database->[5];
$nr = ${$sheet}{"maxrow"};
