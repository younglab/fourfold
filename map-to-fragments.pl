#!/usr/bin/perl

use strict;
use Spreadsheet::Read;

sub torpm {
  my ($in,$out,$mapped) = @_;
  
  open(I,"<",$in) or die "Cannot read $in: $!";
  open(O,">",$out) or die "Cannot write to $out: $!";
  
  while(<I>) {
    if(/^track/) {
      print O;
      next;
    }
    
    if(/^variableStep/) {
      print O;
      next;
    }
    
    chomp;
    
    my ($pos,$val) = split /\t/;
    
    print O "$pos\t" . sprintf("%.3f",$val/$mapped*1e6) . "\n";
  }
  
  close(I);
  close(O);
}

sub executebootstrap {
  my ($basedir,$infile,$outfile,$num) = @_;
  
  my $output = `Rscript $basedir/bootstrap.r $infile $outfile $num`;
  die "Failed to run bootstrap script: $output" unless $? == 0;
}

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
  
  $viewpointchrom = "chr$viewpointchrom" unless $viewpointchrom =~ /^chr/;
  
  my $fragmentfile = "$e1-$e2-$organism/fragments.txt";
  die "Cannot find $fragmentfile! Make sure to run re-fragment-identification.pl first!" unless( -e $fragmentfile);
  
  my $output = `$basedir/mapping-from-bam-file $fragmentfile bamfiles/$name.sorted.bam wigfiles/$name.raw.wig wigfiles/$name.filtered.wig $name bootstrap/$name.raw.counts.txt bootstrap/$name.filtered.counts.txt bootstrap/$name.raw.rpm.txt bootstrap/$name.filtered.rpm.txt stats/$name.txt $viewpointchrom $readstart $readend`;
  
  die "Error in mapping fragments, output is: $output" unless $? == 0;
  
  open(O,"<","stats/$name.txt") or die "Cannot read stats file! $!";
  <O>; <O>; <O>; <O>; <O>; <O>; <O>; <O>; <O>; # skip the 9 lines
  my $l = <O>; ## only need 10th line (# of mapped reads)
  my (undef,$num) = split /: /, $l;
  close(O);
  
  torpm("wigfiles/$name.raw.wig","wigfiles/$name.raw.rpm.wig",$num);
  torpm("wigfiles/$name.filtered.wig","wigfiles/$name.filtered.rpm.wig",$num);
  
  executebootstrap($basedir,"bootstrap/$name.raw.counts.txt","bootstrap/$name.raw.counts.bootstrap.txt",1);
  executebootstrap($basedir,"bootstrap/$name.filtered.counts.txt","bootstrap/$name.filtered.counts.bootstrap.txt",1);
  executebootstrap($basedir,"bootstrap/$name.raw.counts.txt","bootstrap/$name.raw.rpm.bootstrap.txt",$num);
  executebootstrap($basedir,"bootstrap/$name.filtered.counts.txt", "bootstrap/$name.filtered.rpm.bootstrap.txt", $num);

  ### compress
  `gzip wigfiles/$name.raw.wig; rm wigfiles/$name.raw.wig`;
  `gzip wigfiles/$name.filtered.wig; rm wigfiles/$name.filtered.wig`;
  `gzip wigfiles/$name.raw.rpm.wig; rm wigfiles/$name.raw.rpm.wig`;
  `gzip wigfiles/$name.filtered.rpm.wig; rm wigfiles/$name.filtered.rpm.wig`;
}

