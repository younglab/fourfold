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
  
my %filegroups;

for( my $i = 0; $i < $nr; $i++ ) { ## row 1 (index 0) is the header line

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,undef,undef,undef,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  
  next if $name =~ /^#/;
  next if $name =~ /^$/;
  print "\tProcessing reads for sample $name...\n";

  die "Somehow cannot find any temporary file for $name" unless defined($tmpfilemap{$name});
  my $tmpseqfile = $tmpfilemap{$name};
  die "Somehow the temporary file $tmpseqfile for $name disappeared" unless -e $tmpseqfile;
  
    
  my $test;
  my $barcodetest = qr/INVALID1234/;
  my $trimlength;
  my $geotrimlength;
  my ($nbarcode,$nprimer) = (0,0);
  my $ofh;
  my $gfh;
  
  
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
  
  open($ofh,">","$name.trimmed.fq") or die "Cannot write $name.trimmed.fq: $!";
  open($gfh,">","$name.geo.fq") or die "Cannot write to $name.geo.fq: $!" if $geo;
  
  my $bowtieidx = $organisms{lc $organism}->[0];
  my $bowtiecmd = "bowtie -n $bowtien -p 8 -k $bowtiek -m $bowtiem -S --chunkmbs 256 --best --strata $bowtieidx $name.trimmed.fq > bamfiles/$name.sam";
  open(S,">","$name.align.sh") or die "Cannot write to shell script: $!";
  print S "$bowtiecmd\n";
  close(S);
  

  
  push @{$filegroups{$tmpseqfile}}, { Name =>$name,
                                      Barcode => $barcode,
                                      Primer => $primer,
                                      Enzyme1Seq => $e1s,
                                      TestRegex => $test,
                                      BarcodeRegex => $barcodetest,
                                      TrimLength => $trimlength,
                                      GeoTrimLength => $geotrimlength,
                                      BowtieCmd => $bowtiecmd,
                                      StatsFile => "stats/$name.txt",
                                      NBarcode => 0,
                                      NPrimer => 0,
                                      OutputFileHandle => \$ofh,
                                      GeoFileHandle => \$gfh };
}

for my $tmpseqfile (keys(%filegroups)) {
  my @samples = @{$filegroups{$tmpseqfile}};

  open(D,"<","$tmpseqfile") or die "Cannot read $tmpseqfile: $!";


  while(<D>) {
    my $line1 = $_;
    my $line2 = <D>;
    my $line3 = <D>;
    my $line4 = <D>;
    
    chomp $line1;
    chomp $line2;
    chomp $line3;
    chomp $line4;


    for my $sref (@samples) { 
        $sref->{NBarcode}++ if $line2 =~ $sref->{BarcodeRegex};
        next unless $line2 =~ $sref->{TestRegex};
        $sref->{NPrimer}++;
      
        my $fh = ${$sref->{OutputFileHandle}};
        
        my $np2 = substr($line2,$sref->{TrimLength});
        my $np4 = substr($line4,$sref->{TrimLength});
        
        print $fh "$line1\n$np2\n$line3\n$np4\n";
        
        if($geo) {
          $fh =  ${$sref->{GeoFileHandle}};
          $np2 = substr($line2,$sref->{GeoTrimLength});
          $np4 = substr($line4,$sref->{GeoTrimLength});
          
          print $fh "$line1\n$np2\n$line3\n$np4\n";
        }
    }
  }
  
  close(D);

  
  for my $sref (@samples) {
    my $fh = ${$sref->{OutputFileHandle}};
    close($fh);
    
    ### right now die if a sample had no reads in it -- may want to handle this a bit more
    ### elegantly in the future to keep processing the other samples in the mean time
    if( $sref->{NPrimer} == 0 ) {
      if($sref->{NBarcode} == 0 ) {
        die "Sample $sref->{Name} did not have any reads for primer $sref->{Primer}! Additionally, no reads with barcode $sref->{Barcode} were detected in the FASTQ file $tmpseqfile as well";
      } else {
        die "Sample $sref->{Name} did not have any reads for primer $sref->{Primer}! But the barcode $sref->{Barcode} was found on some reads in the FASTQ file $tmpseqfile";
      }
    }

    open(O,">",$sref->{StatsFile}) or die "Cannot write to $sref->{StatsFile}";
    print O "Sample ID: $sref->{Name}\n";
    print O "Number of sequenced reads: $tmpfilesize{$tmpseqfile}\n";
    print O "Number of reads with barcode $sref->{Barcode}: $sref->{NBarcode}\n";
    print O "Number of reads with barcode and primer $sref->{Barcode}$sref->{Primer}: $sref->{NPrimer}\n";
    print O "Bowtie command: $sref->{BowtieCmd}\n";
    close(O);
  
    if( $geo ) {
      $fh = ${$sref->{GeoFileHandle}};
      close($fh);
      `gzip $sref->{Name}.geo.fq`;
      `md5sum $sref->{Name}.geo.fq.gz > $sref->{Name}.geo.fq.gz.md5sum`;
    }
  }
}

### cleanup tmp files
for my $file (@tmpfiles) {
  unlink($file);
}

