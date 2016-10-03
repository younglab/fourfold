#!/usr/bin/perl

use strict;
use Spreadsheet::Read;

if( scalar(@ARGV) < 3 ) {
  die "Need sample table and format type and index!";
}

my ($sampletable,$format,$bowtieidx) = @ARGV;

die "Cannot find file $sampletable!" unless -e $sampletable;

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

if( !defined($sheet) || !defined($nr) ) {
  print "ERROR: Cannot seem to read $sampletable as an Excel file\n";
  exit 1;
}

die "Old format not implemented yet" if $format eq "old";

open(S,">",".tmp.sampletable") or die "Cannot write to .tmp.sampletable: $!";

for( my $i = 1; $i < $nr; $i++ ) { ## row 1 (index 0) is the header line

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$chr,$start,$end,$fastq,$barcode,$primer,$e1,$e2) = @arr;
  
  print "\tProcessing sample $name...\n";
  
  print S "$name\t$chr\t$start\t$end\t$e1\t$e2\n";
  
  next if $name =~ /^$/;
  
  my $origfastq = $fastq;
  $fastq =~ s/\\\\WI-FILES/\/nfs/;
  $fastq =~ s/\\\\WI-HTDATA/\/lab/;
  $fastq =~ s/\\/\//g;

  
  die "Cannot find file $fastq (original name $origfastq)!" unless -e $fastq;
  
  ### extract sequence
  if( $fastq =~ /[.]gz$/ && $fastq !~ /[.]tar[.]gz$/ ) {
    `zcat $fastq > .tmp.seq.fq`;
  } elsif( $fastq =~ /[.]bz2$/ ) {
    `bzcat $fastq > .tmp.seq.fq`;
  } elsif( $fastq =~ /[.]tar[.]gz$/ ) {
    `tar --strip-components=5 -xzOf $fastq > .tmp.seq.fq`;
  } else {
    `cat $fastq > .tmp.seq.fq`;
  }
  
  
  open(D,"<",".tmp.seq.fq") or die "Cannot read .tmp.seq.fq: $!";
  open(P,">",".tmp.primer.fq") or die "Cannot write .tmp.primer.fq: $!";
  
  my $test;
  my $trimlength;
  
  if($barcode ne "NA") {
    $test = qr/^$barcode$primer/;
    $trimlength = length($barcode) + length($primer);
  } else {
    $test = qr/^$primer/;
    $trimlength = length($primer);
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
  
  `bowtie -n 1 $bowtieidx -k 1 -m 1 -S --chunkmbs 256 --best --strata $name.trimmed.fq > $name.sam`;
  `samtools view -Sb $name.sam > $name.bam`;
  `samtools sort -@ 6 -Ttmp $name.bam > $name.sorted.bam`;
  unlink("$name.sam");
  unlink("$name.bam");
}

close(S);

