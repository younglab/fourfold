#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use File::Basename;
use FourCOpts::OrganismDatabase qw(loadorgdatabase);


if( scalar(@ARGV) < 6 ) {
  die "Need sample table and organism index and bowtie options!";
}

my ($sampletable,$organismdatabase,$bowtien,$bowtiek,$bowtiem,$geostr) = @ARGV;

my $geo = $geostr eq "yes";

die "Cannot find file $sampletable!" unless -e $sampletable;
die "Cannot find organism database $organismdatabase!" unless -e $organismdatabase;

my %organisms = %{loadorgdatabase($organismdatabase)};

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

print "\tRetrieving raw FASTQ files...\n";

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
    my $bfastq = basename($fastq);
    my $output = `wget $fastq 2>logs/$bfastq.log`;
    die "wget failed: $output" if $? != 0;
    $fastq = basename($fastq);
    $output = `fastq-dump --gzip $bfastq`;
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
  
  die "Failed to extract FASTQ data from $fastq, see error messages" unless $? == 0;
  
  push @tmpfiles, $tmpseqfile;
  my $nlines = `wc -l $tmpseqfile`;
  chomp $nlines;
  my $tl = (split /\s+/, $nlines)[0];
  $tmpfilesize{$tmpseqfile} = $tl/4;
}
  


for( my $i = 0; $i < $nr; $i++ ) { ## row 1 (index 0) is the header line

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  
  next if $name =~ /^#/;
  next if $name =~ /^$/;
  print "\tProcessing reads for sample $name...\n";

  die "Somehow cannot find any temporary file for $name" unless defined($tmpfilemap{$name});
  my $tmpseqfile = $tmpfilemap{$name};
  die "Somehow the temporary file $tmpseqfile for $name disappeared" unless -e $tmpseqfile;

  
  open(D,"<","$tmpseqfile") or die "Cannot read $tmpseqfile: $!";
  open(P,">",".tmp.primer.fq") or die "Cannot write .tmp.primer.fq: $!";
  if($geo) {
    open(G,">",".tmp.geo.fq") or die "Cannot write to .tmp.geo.fq: $!";
  }
  
  my $test;
  my $barcodetest = qr/INVALID1234/;
  my $trimlength;
  my $geotrimlength;
  my ($nbarcode,$nprimer) = (0,0);
  
  if($barcode ne "NA") {
    $test = qr/^$barcode$primer/i;
    $barcodetest = qr/^$barcode/i;
    $trimlength = length($barcode) + length($primer) - length($e1s);
    $geotrimlength = length($barcode);
  } else {
    $test = qr/^$primer/i;
    $trimlength = length($primer) - length($e1s);
    $geotrimlength = 0;
  }
  
  my @fastqlines = <D>;
  chomp @fastqlines;
  
  while(scalar(@fastqlines)>0) {
      my $line1 = shift @fastqlines;
      my $line2 = shift @fastqlines;
      my $line3 = shift @fastqlines;
      my $line4 = shift @fastqlines;
      
      $nbarcode++ if $line2 =~ $barcodetest;
      next unless $line2 =~ $test;
      $nprimer++;
      
      my $np2 = substr($line2,$trimlength);
      my $np4 = substr($line4,$trimlength);
      
      print P "$line1\n$np2\n$line3\n$np4\n";
      
      if($geo) {
        $np2 = substr($line2,$geotrimlength);
        $np4 = substr($line4,$geotrimlength);
        
        print G "$line1\n$np2\n$line3\n$np4\n";
      }
  }
  
  close(P);
  close(D);
  close(G) if $geo;
  
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
  if( $geo ) {
    rename(".tmp.geo.fq","$name.geo.fq");
    `gzip $name.geo.fq`;
    `md5sum $name.geo.fq.gz > $name.geo.fq.gz.md5sum`;
  }
  
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

