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

sub converttobigwig {
  my ($origfile,$bwfile) = @_;
}

sub tomegabaserpm {
  my ($in,$out,$chr,$ostart,$oend) = @_;
  
  my ($start,$end) = ($ostart-1e6,$oend+1e6);
  
  open(I,"<",$in) or die "Cannot read $in: $!";
  open(O,">",$out) or die "Cannot write to $out: $!";
  
  my @entries;
  my $onchrom = 0;
  my $megabasecounts = 0;
  
  while(<I>) {
    if(/^track/) {
      print O;
      next;
    }
    
    if(/^variableStep.*chrom=(chr\w+)/) {
      if($1 eq $chr) {
        print O;
        $onchrom = 1;
      } else {
        $onchrom = 0;
      }
      next;
    }
    
    chomp;
    
    my ($pos,$val) = split /\t/;
    
    next unless $onchrom;
    next unless ($start <= $pos && $pos <= $end);
    
    $megabasecounts += $val;
    
    #print O "$pos\t" . sprintf("%.3f",$val/$mapped*1e6) . "\n";
    push @entries, [$pos,$val];
  }
  
  close(I);
  
  print O "$_->[0]\t" . sprintf("%.3f",$_->[1]/$megabasecounts*1e6) . "\n" for(@entries);
  
  close(O);
}

sub executebootstrap {
  my ($basedir,$infile,$outfile,$num) = @_;
  
  my $output = `Rscript $basedir/bootstrap.r $infile $outfile $num`;
  die "Failed to run bootstrap script: $output" unless $? == 0;
}

if(scalar(@ARGV)<4) {
  die "map-to-fragments.pl <sample table> <basedir> <scriptdir> <organism database>";
}


my ($sampletable,$basedir,$scriptdir,$organismdatabase) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $organismdatabase!" unless -e $organismdatabase;

my %organisms;

open(D,"<","$organismdatabase") or die "Cannot read $organismdatabase: $!";

while(<D>) {
  chomp;
  
  my ($id,$bowtie,$fasta,$chromsizes) = split /\t/;
  
  next if $id =~ /^$/;
  
  die "In organism database, cannot find bowtie index for $id!" unless -e "$bowtie.1.ebwt";
  die "In organism database, cannot find FASTA file for $id!" unless -e $fasta;
  die "In organism database, cannot find chromosome sizes file for $id!" unless -e $chromsizes;
  
  for(split(/,/,$id)) {
    $organisms{lc $_} = [$bowtie,$fasta,$chromsizes];
  }
}
close(D);

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my %stats;
my %sampletypes;

for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$sampletype,$condition,$rep,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;

  next if $name =~ /^#/;
  
  my $chromsizes = $organisms{lc $organism}->[2];

  $viewpointchrom = "chr$viewpointchrom" unless $viewpointchrom =~ /^chr/;
  
  my $fragmentfile = "$e1-$e2-$organism/fragments.txt";
  die "Cannot find $fragmentfile! Make sure to run re-fragment-identification.pl first!" unless( -e $fragmentfile);
  
  my $output = `$basedir/mapping-from-bam-file $fragmentfile bamfiles/$name.sorted.bam wigfiles/$name.raw.wig wigfiles/$name.filtered.wig $name bootstrap/$name.raw.counts.txt bootstrap/$name.filtered.counts.txt bootstrap/$name.raw.rpm.txt bootstrap/$name.filtered.rpm.txt stats/$name.txt $viewpointchrom $readstart $readend`;
  
  my $sampleconditionkey = "$sampletype-$condition";
  
  $sampletypes{$sampleconditionkey} = [[],[],[],[],$chromsizes] unless defined($sampletypes{$sampleconditionkey});
  
  push @{$sampletypes{$sampleconditionkey}->[0]}, "bootstrap/$name.raw.counts.txt";
  push @{$sampletypes{$sampleconditionkey}->[1]}, "bootstrap/$name.filtered.counts.txt";
  push @{$sampletypes{$sampleconditionkey}->[2]}, "bootstrap/$name.raw.rpm.txt";
  push @{$sampletypes{$sampleconditionkey}->[3]}, "bootstrap/$name.filtered.rpm.txt";
   
  die "Error in mapping fragments, output is: $output" unless $? == 0;
  
  open(O,"<","stats/$name.txt") or die "Cannot read stats file! $!";
  <O>; <O>; <O>; <O>; <O>; <O>; <O>; <O>; <O>; # skip the 9 lines
  my $l = <O>; ## only need 10th line (# of mapped reads)
  my (undef,$num) = split /: /, $l;
  close(O);
  
  torpm("wigfiles/$name.raw.wig","wigfiles/$name.raw.rpm.wig",$num);
  torpm("wigfiles/$name.filtered.wig","wigfiles/$name.filtered.rpm.wig",$num);
  
#  tomegabaserpm("wigfiles/$name.raw.wig","wigfiles/$name.raw.MB.rpm.wig",$viewpointchrom,$viewpointstart,$viewpointend);
#  tomegabaserpm("wigfiles/$name.filtered.wig","wigfiles/$name.filtered.MB.rpm.wig",$viewpointchrom,$viewpointstart,$viewpointend);
  
  executebootstrap($scriptdir,"bootstrap/$name.raw.counts.txt","bootstrap/$name.raw.counts.bootstrap.txt",1);
  executebootstrap($scriptdir,"bootstrap/$name.filtered.counts.txt","bootstrap/$name.filtered.counts.bootstrap.txt",1);
  executebootstrap($scriptdir,"bootstrap/$name.raw.counts.txt","bootstrap/$name.raw.rpm.bootstrap.txt",$num);
  executebootstrap($scriptdir,"bootstrap/$name.filtered.counts.txt", "bootstrap/$name.filtered.rpm.bootstrap.txt", $num);
  
  `wigToBigWig wigfiles/$name.raw.wig $chromsizes wigfiles/$name.raw.bw`;
  `wigToBigWig wigfiles/$name.filtered.wig $chromsizes wigfiles/$name.filtered.bw`;
  `wigToBigWig wigfiles/$name.raw.rpm.wig $chromsizes wigfiles/$name.raw.rpm.bw`;
  `wigToBigWig wigfiles/$name.filtered.rpm.wig $chromsizes wigfiles/$name.filtered.rpm.bw`;


  ### compress
  `gzip wigfiles/$name.raw.wig;`;
  `gzip wigfiles/$name.filtered.wig;`;
  `gzip wigfiles/$name.raw.rpm.wig;`;
  `gzip wigfiles/$name.filtered.rpm.wig; `;
}

for my $sck (keys(%sampletypes)) {
  next if scalar($sampletypes{$sck}->[0]) < 2; ## pass on samples with only one replicate
  
  my @outputfiles = ("bootstrap/$sck.raw.counts.txt","bootstrap/$sck.filtered.counts.txt","bootstrap/$sck.raw.rpm.txt","bootstrap/$sck.filtered.rpm.txt");
  my @wigoutputfiles = ("wigfiles/$sck.raw.counts.wig","wigfiles/$sck.filtered.counts.wig","wigfiles/$sck.raw.rpm.wig","wigfiles/$sck.filtered.rpm.wig");
  my @bwoutputfiles = ("wigfiles/$sck.raw.counts.bw","wigfiles/$sck.filtered.counts.bw","wigfiles/$sck.raw.rpm.bw","wigfiles/$sck.filtered.rpm.bw");


  my $chromsizes = $sampletypes{$sck}->[4];
  
  for( my $i = 0; $i <= $#outputfiles; $i++ ) {
    my $ofile = $outputfiles[$i];
    my $sfiles = join(" ",@{$sampletypes{$sck}->[$i]});
    my $wfile = $wigoutputfiles[$i];
    my $bwfile = $bwoutputfiles[$i];
    
    my $output = `Rscript $scriptdir/combine-profiles.r $ofile $sfiles`;
    die "Failed to combine files in $sck" unless $? == 0;
    
    open(T,"<",$ofile) or die "Cannot read $ofile: $!";
    open(W,">",$wfile) or die "Cannot write to $wfile: $!";
    
    print W "track type=wiggle_0 name=\"$sck\" description=\"$sck\"\n";
    my $lastchr = "";
    
    while(<T>) {
      chomp;
      my ($chr,$pos,$s) = split /\t/;
      
      unless($lastchr eq $chr) {
        print W "variableStep chrom=$chr\n";
        $lastchr = $chr;
      }
      
      print W "$pos\t$s\n";
      
      `wigToBigwig $wfile $chromsizes $bwfile`;
    }
    
    close(W);
    close(T);
  }
}

