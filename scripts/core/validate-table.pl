#!/usr/bin/perl

use strict;
use File::Basename;
use Spreadsheet::Read;
use FourCOpts::OrganismDatabase qw(loadorgdatabase);

if(scalar(@ARGV)<3) {
  die "Syntax: ./validate-table.pl <sample table> <organism database> <validate only>";
}

my ($sampletable,$organismdatabase,$validateonly) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $organismdatabase!" unless -e $organismdatabase;

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

if( !defined($sheet) || !defined($nr) ) {
  print "ERROR: Cannot seem to read $sampletable as an Excel file\n";
  exit 1;
}

my %organisms = %{loadorgdatabase($organismdatabase)};

my %names;
my %viewpoints;
my %repairs;

for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  next if $name =~ /^#/;
  next if $name =~ /^$/;
  
  print "\tChecking $name...\n";
  
  die "The sample identifier ($name) cannot have any whitespace characters within it" if $name =~ /\s/;
  die "The sample identifier ($name) cannot have any forward or backslash characters within it" if $name =~ /[\/\\]/;

  
  die "$name is duplicated in the table, line " . ($i+1) . "!" if defined($names{$name});
  $names{$name} = 1;
  
  die "$organism for $name doesn't match any entry within the organism database" unless defined($organisms{lc $organism});
  
  die "Please set viewpoint information for sample $name to \"NA\" if not in use" if( $viewpointid =~ /^$/ || $viewpointchrom =~ /^$/ || $viewpointstart =~ /^$/ || $viewpointend =~ /^$/ || $readstart =~ /^$/ || $readend =~ /^$/ );
  
  if(defined($viewpoints{$viewpointid})) {
    my @d = @{$viewpoints{$viewpointid}};
    my @s = ($viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend);
    
    die "Differences detected for $viewpointid for sample $name across the samples in the table" unless $d[0] eq $s[0];

    
    for(my $i = 1; $i < $#d; $i++ ) {
      die "Differences detected for $viewpointid for sample $name across the samples in the table" unless $d[$i] == $s[$i];
    }
    
  } else {
    if( $viewpointid eq "NA") {
      warn "For sample $name, viewpoint consistentcy checking is turned off!";
    } else {
      $viewpoints{$viewpointid} = [$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend];
    }
  }
  
  my $origfastq = $fastq;
  if( $fastq =~ /^ftp:\/\//i ) {
    my $basename = basename($fastq);
    
    my $cmd = "wget --spider $fastq 2>&1";
    $cmd = "wget --spider $fastq >logs/$basename.wget.validate.txt 2>&1" unless $validateonly;
    
    my $output = `$cmd`;
    die "$name does not seem to exist at $fastq (output: $output)" unless $? == 0;
  } else {
    $fastq =~ s/\\\\WI-FILES/\/nfs/;
    $fastq =~ s/\\\\WI-HTDATA/\/lab/;
    $fastq =~ s/\\/\//g;
    die "Cannot find file $fastq (original name $origfastq)!" unless -e $fastq;
    die "$fastq appears to be a directory, not a FASTQ file (original name $origfastq)!" if -d $fastq;
  }
  
  die "Barcode of $name doesn't only contain DNA characters or NA" unless ($barcode =~ /^[ACTG]+$/i || $barcode eq "NA");
  die "Primer sequence of $name doesn't only contain DNA characters" unless ($primer =~ /^[ACTG]+$/i);
  
  my $pair = "$e1-$e2";
  my $s = "$e1s-$e2s";

  
  if(defined($repairs{$pair})) {
    die "Restriction enzyme sequences not the same in $name compared with previous samples" unless $s eq $repairs{$pair};
  } else {
    $repairs{$pair} = $s;
  }
}


