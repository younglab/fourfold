#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use FourCOpts::OrganismDatabase qw(loadorgdatabase);
use FourCOpts::WigFile qw(writewigfile);

sub in {
  my ($test,$ref) = @_;
  
  my $ret = 0;
  
  if(scalar(@{$ref})==1) {
    $ret = 1 if ($ref->[0] eq "all" || $test =~ m/$ref->[0]/);
  } else {
  
    for(@{$ref}) {
      if($test eq $_) {
        $ret = 1;
        last;
      }
    }
  }
  
  return $ret;
}

if(scalar(@ARGV)<7) {
  die "Arguments <template> <organism database> <basedir> <pseudocount> <input dir> <output dir> <file 1> [file2...]";
}

my ($sampletable,$organismdatabase,$basedir,$pseudocount,$inputdir,$outputdir,@files) = @ARGV;
my $binsize = 5000;
my $stepsize = 50;

die "Cannot find file $sampletable!" unless -e $sampletable;
die "Cannot find organism database $organismdatabase!" unless -e $organismdatabase;

my %organisms = %{loadorgdatabase($organismdatabase)};

#open(D,"<","$organismdatabase") or die "Cannot read $organismdatabase: $!";

#while(<D>) {
#  chomp;
  
#  my ($id,$bowtie,$fasta) = split /\t/;
  
#  next if $id =~ /^$/;
  
#  die "In organism database, cannot find bowtie index for $id!" unless -e "$bowtie.1.ebwt";
#  die "In organism database, cannot fine FASTA file for $id!" unless -e $fasta;
  
#  for(split(/,/,$id)) {
#    $organisms{lc $_} = [$bowtie,$fasta];
#  }
#}
#close(D);

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

if( !defined($sheet) || !defined($nr) ) {
  print "ERROR: Cannot seem to read $sampletable as an Excel file\n";
  exit 1;
}

my %samplegroups;
my %sampleorg;

for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$sampletype,$condition,$rep,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;

  next if $name =~ /^#/;
  
  next unless in($name,\@files);
  
  my $samplekey = "$sampletype-$condition";
  
  push @{$samplegroups{$samplekey}}, "$inputdir/$name.filtered.rpm.txt";
  $sampleorg{$samplekey} = $organism; ### expand this later
}

my @sgroups = keys(%samplegroups);

for( my $i = 0; $i <= $#sgroups; $i++ ) {
  for( my $j = $i+1; $j <= $#sgroups; $j++) {
    my ($group1,$group2) = @sgroups[$i,$j];
    
    my $outputfile1 = "$outputdir/$group1-$group2-comparison.txt";
    my $outputsmoothedfile1 = "$outputdir/$group1-$group2-comparison.smoothed.txt";

    my $outputfile2 = "$outputdir/$group2-$group1-comparison.txt";
    my $outputsmoothedfile2 = "$outputdir/$group2-$group1-comparison.smoothed.txt";


    my $f1 = join(" ",@{$samplegroups{$group1}});
    my $f2 = join(" ",@{$samplegroups{$group2}});
    
    my $org = $sampleorg{$group1}; ## expand this later
    
    my $output = `Rscript $basedir/../scripts/comparison/comparison.r $outputfile1 $outputsmoothedfile1 $organisms{$org}->[2] $binsize $stepsize $pseudocount $f1 SEP $f2`;
    
    die "Failed to compare samples: $output" unless $? == 0;
    
    $output = `Rscript $basedir/../scripts/comparison/comparison.r $outputfile2 $outputsmoothedfile2 $organisms{$org}->[2] $binsize $stepsize $pseudocount $f2 SEP $f1`;
    
    die "Failed to compare samples: $output" unless $? == 0;
    
    writewigfile($outputsmoothedfile1,"$outputdir/$group1-$group2-comparison-smoothed.wig","$group1-$group2-comparison smoothed");
    writewigfile($outputsmoothedfile2,"$outputdir/$group2-$group1-comparison-smoothed.wig","$group2-$group1-comparison smoothed");
    
    my $wigfile = "$outputdir/$group1-$group2-comparison.wig";
    open(O,"<",$outputfile1) or die "Cannot read $outputfile1: $!";
    open(W,">",$wigfile) or die "Cannot write to $wigfile: $!";
    
    print W "track type=wiggle_0 name=\"$group1-$group2-comparison\" description=\"$group1-$group2-comparison\"\n";
    
    my $lastchr = '';
    while(<O>) {
      chomp;
      
      my ($chr,$pos,$signal) = split /\t/;
      
      unless($lastchr eq $chr) {
        print W "variableStep chrom=$chr\n";
        $lastchr =$chr;
      }
      print W "$pos\t$signal\n";
    }
    close(O);
    close(W);
    
    $wigfile = "$outputdir/$group2-$group1-comparison.wig";
    open(O,"<",$outputfile2) or die "Cannot read $outputfile2: $!";
    open(W,">",$wigfile) or die "Cannot write to $wigfile: $!";
    
    print W "track type=wiggle_0 name=\"$group2-$group1-comparison\" description=\"$group2-$group1-comparison\"\n";
    
    $lastchr = '';
    while(<O>) {
      chomp;
      
      my ($chr,$pos,$signal) = split /\t/;
      
      unless($lastchr eq $chr) {
        print W "variableStep chrom=$chr\n";
        $lastchr =$chr;
      }
      print W "$pos\t$signal\n";
    }
    close(O);
    close(W);
  }
}

