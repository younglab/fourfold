#!/usr/bin/perl

use strict;
use Spreadsheet::Read;

if(scalar(@ARGV)<1) {
  die "map-to-fragments.pl <sample table> <basedir>";
}

my ($sampletable,$basedir) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my %stats;

for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;

  next if $name =~ /^#/;
  
  my $fragmentfile = "$e1-$e2-$organism/fragments.txt";
  die "Cannot find $fragmentfile! Make sure to run re-fragment-identification.pl first!" unless( -e $fragmentfile);
  
  `$basedir/mapping-from-bam-file $fragmentfile bamfiles/$name.sorted.bam wigfiles/$name.raw.wig wigfiles/$name.filtered.wig stats/$name.stats.txt $viewpointchrom $readstart $readend`;
}

