#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
#use Parallel::ForkManager;
use File::Basename;


if( scalar(@ARGV) < 2 ) {
  die "Need sample table and organism index!";
}

my ($sampletable,$organismdatabase) = @ARGV;

die "Cannot find file $sampletable!" unless -e $sampletable;
die "Cannot find organism database $organismdatabase!" unless -e $organismdatabase;

my %organisms;

open(D,"<","$organismdatabase") or die "Cannot read $organismdatabase: $!";

while(<D>) {
  chomp;
  
  my ($id,$bowtie,$fasta) = split /\t/;
  
  next if $id =~ /^$/;
  
  die "In organism database, cannot find bowtie index for $id!" unless -e "$bowtie.1.ebwt";
  die "In organism database, cannot fine FASTA file for $id!" unless -e $fasta;
  
  for(split(/,/,$id)) {
    $organisms{lc $_} = [$bowtie,$fasta];
  }
}
close(D);

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

if( !defined($sheet) || !defined($nr) ) {
  print "ERROR: Cannot seem to read $sampletable as an Excel file\n";
  exit 1;
}


my @tmpfiles;
#my $pm = Parallel::ForkManager->new(50);

### unpack sequences first
for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,$fastq) = @arr;
  
  
  next if $name =~ /^#/;
  next if $name =~ /^$/;
  
  
  my $origfastq = $fastq;
  my $basefastqname = basename($fastq);
  my $tmpseqfile = ".tmp.$basefastqname.seq.fq";
  
  next if -e $tmpseqfile;
  if( $fastq =~ /^ftp:\/\// ) {
    my $output = `wget $fastq`;
    die "wget failed: $output" if $? != 0;
    $fastq = basename($fastq);
    $output = `fastq-dump --gzip $fastq`;
    die "fastq-dump failed: $output" if $? != 0;
    $fastq =~ s/[.]sra/.fastq.gz/;
  } else {
    $fastq =~ s/\\\\WI-FILES/\/nfs/;
    $fastq =~ s/\\\\WI-HTDATA/\/lab/;
    $fastq =~ s/\\/\//g;
  }

  
  die "Cannot find file $fastq (original name $origfastq)!" unless -e $fastq;
  
  ### extract sequence
  if( $fastq =~ /[.]gz$/i && $fastq !~ /[.]tar[.]gz$/i ) {
    `zcat $fastq > $tmpseqfile`;
  } elsif( $fastq =~ /[.]bz2$/i ) {
    `bzcat $fastq > $tmpseqfile`;
  } elsif( $fastq =~ /[.]tar[.]gz$/i ) {
    `tar --strip-components=5 -xzOf $fastq > $tmpseqfile`;
  } else {
    `cat $fastq > $tmpseqfile`;
  }
  
  push @tmpfiles, $tmpseqfile;
}
  


for( my $i = 0; $i < $nr; $i++ ) { ## row 1 (index 0) is the header line

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  
  next if $name =~ /^#/;
  next if $name =~ /^$/;
  print "\tProcessing sample $name...\n";


  my $basefastqname = basename($fastq);
  my $tmpseqfile = ".tmp.$basefastqname.seq.fq";

  
  die "Cannot find file $fastq (temporary name $tmpseqfile)!" unless -e $tmpseqfile;

  
  open(D,"<","$tmpseqfile") or die "Cannot read $tmpseqfile: $!";
  open(P,">",".tmp.primer.fq") or die "Cannot write .tmp.primer.fq: $!";
  
  my $test;
  my $trimlength;
  
  if($barcode ne "NA") {
    $test = qr/^$barcode$primer/i;
    $trimlength = length($barcode) + length($primer) - length($e1s);
  } else {
    $test = qr/^$primer/i;
    $trimlength = length($primer) - length($e1s);
  }
  
  while(<D>) {
      my $line1 = $_;
      my $line2 = <D>;
      my $line3 = <D>;
      my $line4 = <D>;
      
      chomp $line1;
      chomp $line2;
      chomp $line3;
      chomp $line4;
      
      next unless $line2 =~ $test;
      
      my $np2 = substr($line2,$trimlength);
      my $np4 = substr($line4,$trimlength);
      
      print P "$line1\n$np2\n$line3\n$np4\n";
  }
  
  close(P);
  close(D);
  
  rename(".tmp.primer.fq","$name.trimmed.fq");
}

### cleanup tmp files
for my $file (@tmpfiles) {
  unlink($file);
}

