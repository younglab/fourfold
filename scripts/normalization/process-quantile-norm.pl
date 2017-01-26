#!/usr/bin/perl

use strict;
use Spreadsheet::Read;
use FourCOpts::OrganismDatabase qw(loadorgdatabase);


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

#print "DEBUG: @ARGV\n";
if(scalar(@ARGV) < 5) {
  die "process-quantile-norm.pl <sample table> <basedir> <organism database> <output dir> <only cis positions?> <samples...>";
}

my ($sampletable,$basedir,$organismdatabase,$outputdir,$cisonly,@samples) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;
die "Cannot find $outputdir!" unless -e $outputdir;

my %organisms = %{loadorgdatabase($organismdatabase)};

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my @norm;
my %samplegroups;
my $chromfile = "";
my $vchrom = "";

for( my $i = 0; $i < $nr; $i++ ) { 
  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my($name,$celltype,$condition,$rep,$organism,$viewpointid,$viewpointchrom,$viewpointstart,$viewpointend,$readstart,$readend,$fastq,$barcode,$primer,$revprimer,$e1,$e2,$e1s,$e2s) = @arr;
  
  next if $name =~ /^#/;
  
  if(in($name,\@samples)) {
    my $cfile = "bootstrap/$name.filtered.counts.txt";
    my $cbfile = "bootstrap/$name.filtered.counts.bootstrap.txt";
    push @norm, [$name,$cfile,$cbfile];
    
    if($cisonly) {
      if( $vchrom eq "" ) {
        $vchrom = $viewpointchrom;
      } elsif( $vchrom ne $viewpointchrom ) {
        die "The viewpoint chromsome is different across all samples when cis-only option is on!";
      }
    }
    
    my $samplekey = "$celltype-$condition";
    $samplegroups{$samplekey} = [] unless defined($samplegroups{$samplekey});
    
    push @{$samplegroups{$samplekey}}, ["$outputdir/$name.filtered.counts.txt","$outputdir/$name.filtered.rpm.txt"];
    
    $chromfile = $organisms{$organism}->[2];
  }
}

die "not enough samples" if scalar(@norm) < 2;

my $exe = "";
$exe .= join(" ",@{$_}) . " " for(@norm);

my $qnormtmp = "$outputdir/quantile-normalized-samples.txt";
my $qnormbtmp = "$outputdir/quantile-normalized-samples-bootstrap.txt";

$vchrom = "NA" unless $cisonly;
my $output = `Rscript $basedir/quantile-normalization.r $outputdir $vchrom $exe 2>&1`;

die "Failed to normalize (cmd: Rscript $basedir/quantile-normalization.r $outputdir $vchrom $exe): $output" unless $? == 0;

open(N,"<","$qnormtmp") or die "Cannot read normalized samples: $!";

<N>; ## skip header

my @output;

for(@norm) {
  my $ref = [$_->[0],[]]; ### place sample name and start array to build output wig file
  push @{$ref->[1]}, "track type=wiggle_0 name=\"$_->[1]QuantNorm\" description=\"$_->[1] Quantile Normalized\"";
  push @output, $ref;
}

my $lastchr = "";

while(<N>) {
  chomp;
  
  my ($chr,$pos,@vals) = split /\t/;
  
  if($chr ne $lastchr) {
    my $s = "variableStep chrom=$chr";
    
    push @{$_->[1]},$s for(@output);
    $lastchr = $chr;
  }
  
  for( my $i = 0; $i <= $#vals; $i++ ) {
    push @{$output[$i]->[1]}, "$pos\t$vals[$i]";
  }
}

close(N);

for(@output) {
  my $outputfile ="$outputdir/$$_[0].filtered.rpm.wig";
  my $bigwigfile = "$outputdir/$$_[0].filtered.rpm.bw";
  open(O,">",$outputfile) or die "Cannot write to $outputfile: $!";
  print O "$_\n" for(@{$_->[1]});
  close(O);
  
  `wigToBigWig $outputfile $chromfile $bigwigfile`;
}

unlink("$qnormtmp");
unlink("$qnormbtmp");

for my $skey (keys(%samplegroups)) { 
  #print "DEBUG: $skey\n";
  my @files = @{$samplegroups{$skey}};
  my @countfiles;
  my @rpmfiles;
  
  next unless scalar(@files)>=2;
  
  for(@files) {
    push @countfiles, $_->[0];
    push @rpmfiles, $_->[1];
  }
  
  
  my $ofile = "$outputdir/$skey.filtered.rpm.txt";
  my $wfile = "$outputdir/$skey.filtered.rpm.wig";
  my $bwfile = "$outputdir/$skey.filtered.rpm.bw";


  my $f = join(" ",@rpmfiles);

  
  my $output = `Rscript $basedir/../core/combine-profiles.r $ofile $f 2>&1`;
  
  die "Failed to combine samples: $output" unless $? == 0;
  
  open(T,"<","$ofile") or die "Cannot read $ofile: $!";
  open(W,">",$wfile) or die "Cannot write to $wfile: $!";
  
  print W "track type=wiggle_0 name=\"$skey-QuantNorm\" description=\"$skey-QuantNorm\"\n";
  
  my $lastchr = "";
  
  while(<T>) {
    chomp;
    
    my ($chr,$pos,$s) = split /\t/;
    
    unless($lastchr eq $chr) {
      print W "variableStep chrom=$chr\n";
      $lastchr =$chr;
    }
    
    print W "$pos\t$s\n";
  }
  
  close(T);
  close(W);
  
  `wigToBigWig $wfile $chromfile $bwfile`;
}
  
