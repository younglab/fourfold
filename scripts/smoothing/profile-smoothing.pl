#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use FourCOpts::OrganismDatabase qw(loadorgdatabase);
use FourCOpts::WigFile qw(writewigfile);

sub in {
  my ($test,$ref) = @_;
  
  my $ret = 0;
  
  if(scalar(@{$ref})==1) {
    $ret = 1 if $test =~ m/$ref->[0]/;
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

if(scalar(@ARGV)<8) {
  die "profile-smoothing.pl <sample table> <basedir> <organism databse> <input dir> <output dir> <bin size> <step size> <samples...>";
}

my ($sampletable,$basedir,$organismdatabase,$inputdir,$outputdir,$binsize,$stepsize,@samples) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find directory $inputdir" unless -e $inputdir;
die "Cannot find $outputdir!" unless -e $outputdir;
die "Cannot find $organismdatabase" unless -e $organismdatabase;

my %organisms = %{loadorgdatabase($organismdatabase)};

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my %samplegroups;
my %sampleorganism;

for( my $i = 0; $i < $nr; $i++ ) { 
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$celltype,$condition,$rep,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  next if $name =~ /^#/;
  
  unless( $samples[0] eq "all") {
    next unless in($name,\@samples);
  }
  
  my $key ="$celltype:$condition";
  
  push @{$samplegroups{$key}}, ["$inputdir/$name.filtered.rpm.txt", "$inputdir/$name.filtered.rpm.bootstrap.txt"];
  
  my $chromsizes = $organisms{$organism}->[2];
  $sampleorganism{$key} = $chromsizes;

  my $output = `Rscript $basedir/smooth-single-profile.r $binsize $stepsize $chromsizes $inputdir/$name.filtered.rpm.txt $inputdir/$name.filtered.rpm.bootstrap.txt $outputdir/$name.filtered.smoothed.rpm.txt 2> /dev/null`;
  
  die "Smoothing failed with an error: $output" unless( $? == 0 );
  
  ### generate a WIG file
  
  my $outwig = "$outputdir/$name.filtered.smoothed.rpm.wig";
  my $outbw = "$outputdir/$name.filtered.smoothed.rpm.bw";

  
  open(T,"<","$outputdir/$name.filtered.smoothed.rpm.txt") or die "Cannot read $outputdir/$name.filtered.smoothed.rpm.txt: $!";
  open(W,">","$outwig") or die "Cannot write $outwig: $!";
  
  print W "track type=wiggle_0 name=\"$name\" description=\"$name\"\n";
  my $lastchr = "";
  while(<T>) {
    chomp;
    
    my ($chr,$start,$end,$s,$lb,$ub) = split /\t/;
    
    unless($lastchr eq $chr) {
      print W "variableStep chrom=$chr span=$stepsize\n";
      $lastchr=$chr;
    }
    
    print W "$start\t$s\n";
  }
  
  close(W);
  close(T);
  
  `wigToBigWig $outwig $chromsizes $outbw`; 
  `gzip $outwig`;
}

for my $group (keys(%samplegroups)) {
  my ($celltype,$condition) = split /:/, $group;
  
  ### code
  
  my @samples = @{$samplegroups{$group}};
  my $chromsizes = $sampleorganism{$group};
  
  my @tmp;
  
  push @tmp, join(" ",@{$_}) for (@samples);
  my $sampleexe = join(" ",@tmp);
  
  my $outfile = "$outputdir/$celltype-$condition.filtered.rpm.txt";
  my $outfilewig = "$outputdir/$celltype-$condition.filtered.rpm.wig";
  my $outfilebw = "$outputdir/$celltype-$condition.filtered.rpm.bw";

  my $sid = "$celltype-$condition";

  my $output = `Rscript $basedir/multi-sample-profile.r $outfile $binsize $stepsize $chromsizes $sampleexe 2>&1`;
  
  die "Smoothing of profile failed: $output" unless $? == 0;
  
  open(O,"<",$outfile) or die "Cannot read $outfile: $!";
  open(W,">",$outfilewig) or die "Cannot write to $outfilewig: $!";
  
  my $lastchr = "";
  
  print W "track type=wiggle_0 name=\"$sid\" description=\"$sid\"\n";
  
  while(<O>) {
    my ($chr,$s,$m) = split /\t/;
    
    
    unless($lastchr eq $chr) {
      print W "variableStep chrom=$chr span=$stepsize\n";
      $lastchr = $chr;
    }
    
    print W "$s\t$m\n";
  }
  
  close(O);
  close(W);
  
  `wigToBigWig $outfilewig $chromsizes $outfilebw`; 
  `gzip $outfilewig`;
}

