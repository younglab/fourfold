#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use File::Basename;


if( scalar(@ARGV) < 5 ) {
  die "Need sample table and organism index and bowtie options!";
}

my ($sampletable,$organismdatabase,$bowtien,$bowtiek,$bowtiem) = @ARGV;

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
my %tmpfilemap;
my %tmpfilesize;

### unpack sequences first
for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,undef,undef,undef,undef,undef,undef,undef,$fastq) = @arr;
  
  
  next if $name =~ /^#/;
  next if $name =~ /^$/;
  
  
  my $origfastq = $fastq;
  my $basefastqname = basename($fastq);
  $basefastqname =~ s/\\/./g;
  $basefastqname =~ s/\//./g;
  $basefastqname =~ s/:/./g;
  my $tmpseqfile = ".tmp.$basefastqname.seq.fq";
  $tmpfilemap{$name} = $tmpseqfile;

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
  my $nlines = `wc -l $tmpseqfile`;
  chomp $nlines;
  $tmpfilesize{$tmpseqfile} = (split /\s+/, $nlines)[0];
}
  


for( my $i = 0; $i < $nr; $i++ ) { ## row 1 (index 0) is the header line

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  
  next if $name =~ /^#/;
  next if $name =~ /^$/;
  print "\tProcessing sample $name...\n";


  #my $basefastqname = basename($fastq);
  #my $tmpseqfile = ".tmp.$basefastqname.seq.fq";
  die "Somehow cannot find any temporary file for $name" unless defined($tmpfilemap{$name});
  my $tmpseqfile = $tmpfilemap{$name};
  die "Somehow the temporary file $tmpseqfile for $name disappeared" unless -e $tmpseqfile;

  
  open(D,"<","$tmpseqfile") or die "Cannot read $tmpseqfile: $!";
  open(P,">",".tmp.primer.fq") or die "Cannot write .tmp.primer.fq: $!";
  
  my $test;
  my $barcodetest = qr/NA/;
  my $trimlength;
  my ($nbarcode,$nprimer) = (0,0);
  
  if($barcode ne "NA") {
    $test = qr/^$barcode$primer/i;
    $barcodetest = qr/^$barcode/;
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
      
      $nbarcode++ if $line2 =~ $barcodetest;
      next unless $line2 =~ $test;
      $nprimer++;
      
      my $np2 = substr($line2,$trimlength);
      my $np4 = substr($line4,$trimlength);
      
      print P "$line1\n$np2\n$line3\n$np4\n";
  }
  
  close(P);
  close(D);
  
  ### right now die if a sample had no reads in it -- may want to handle this a bit more
  ### elegantly in the future to keep processing the other samples in the mean time
  if( $nprimer == 0 ) {
    if($nbarcode == 0 ) {
      die "Sample $name did not have any reads for primer $primer! Additionally, no reads with barcode $barcode were detected in the FASTQ file $fastq as well";
    } else {
      die "Sample $name did not have any reads for primer $primer! But the barcode $barcode was found on some reads in the FASTQ file $fastq";
    }
  }
  
  rename(".tmp.primer.fq","$name.trimmed.fq");
  
  my $bowtieidx = $organisms{lc $organism}->[0];
  my $bowtiecmd = "bowtie -n $bowtien -p 8 -k $bowtiek -m $bowtiem -S --chunkmbs 256 --best --strata $bowtieidx $name.trimmed.fq > bamfiles/$name.sam";
  open(S,">","$name.align.sh") or die "Cannot write to shell script: $!";
  print S "$bowtiecmd\n";
  close(S);
  
  open(O,">","stats/$name.txt") or die "Cannot write to stats/$name.out";
  print O "Sample ID: $name\n";
  print O "Number of sequenced reads: $tmpfilesize{$tmpseqfile}\n";
  print O "Number of reads with barcode $barcode: $nbarcode\n";
  print O "Number of reads with barcode and primer $barcode$primer: $nprimer\n";
  print O "Bowtie command: $bowtiecmd\n";
  close(O);
}

### cleanup tmp files
for my $file (@tmpfiles) {
  unlink($file);
}

