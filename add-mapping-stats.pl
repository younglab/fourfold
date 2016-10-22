#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use File::Basename;

if(scalar(@ARGV)<1) {
  die "need sample table";
}

my ($sampletable) = @ARGV;

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

if( !defined($sheet) || !defined($nr) ) {
  print "ERROR: Cannot seem to read $sampletable as an Excel file\n";
  exit 1;
}


for( my $i = 0; $i < $nr; $i++ ) { ## row 1 (index 0) is the header line

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  open(S,">>","stats/$name.out") or die "Cannot appened to stats/$name.out: $!";
  open(L,"<","logs/$name.align.log") or die "Cannot read logs/$name.align.log: $!";
  
  <L>; # skip first line
  my $line = <L>;
  chomp $line;
  my $mapped = (split /: /,$line)[1];
  $line = <L>;
  chomp $line;
  my $unmapped = (split /: /,$line)[1];
  $line = <L>;
  chomp $line;
  my $discarded = (split /: /,$line)[1];
  
  close(L);
  
  print S "Mapped reads: $mapped\n";
  print S "Unmapped reads: $unmapped\n";
  print S "Discarded repetitively mapped reads: $discarded\n";
  
  close(S);
}