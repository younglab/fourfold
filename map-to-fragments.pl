#!/usr/bin/perl

use strict;
use Spreadsheet::Read;

if(scalar(@ARGV)<1) {
  die "map-to-fragments.pl <sample table>";
}

my ($sampletable) = @ARGV;

die "Cannot find $sampletable!" unless -e $sampletable;

my $database = ReadData($sampletable);
my $sheet = $database->[1]; ## get first spreadsheet
my $nr = ${$sheet}{"maxrow"};

my %stats;

for( my $i = 0; $i < $nr; $i++ ) { 

  my @arr = Spreadsheet::Read::cellrow($sheet,$i+1);
  
  my ($samplename,undef,undef,undef,$vchr,$vstart,$vend,undef,undef,$primer,$re1,$re2,$re1seq,$re2seq) = @arr;

  next if $samplename =~ /^#/;
  
  my $fragmentfile = "$re1-$re2/fragments.txt";
  die "Cannot find $fragmentfile! Make sure to run re-fragment-identification.pl first!" unless( -e $fragmentfile);
  
  my %fragments;
  
  print "$samplename: reading fragments...\n";
  open(F,"<",$fragmentfile) or die "Cannot read $fragmentfile: $!";
  while(<F>) {
    chomp;
    
    my ($chr,$s,$e,$l,$r) = split /\t/;
    
    unless(defined($fragments{$chr})) {
      $fragments{$chr} = [];
    }
    
    #my $pos = "$s-$e";
    
    push @{$fragments{$chr}}, [$s,$e,$l,$r,0,0];
  }
  close(F);
  
  my $lastchr = "NA";
  my $posref;
  my $posidx;
  my $previdx;
  my $possize = 0;
  
  print "$samplename: mapping reads on the positive strand...\n";
  open(S,"samtools view -F 0x14 $samplename.sorted.bam |") or die "Cannot read $samplename.sorted.bam! $!";
  while(<S>) { ### all mapped in the positive direction
    chomp;
    
    my (undef,undef,$chr,$mpos,undef,undef,undef,undef,undef,$sequence) = split /\t/;
    
    unless(defined($fragments{$chr})) {
      warn "$chr appeared in the mapping but not in the fragments directory";
      next; ## ignore for now
    }
    
    if($lastchr ne $chr) {
      $lastchr = $chr;
      $posref = $fragments{$chr};
      $posidx = 0;
      $previdx= 0;
      $possize = scalar(@{$posref});
    }
    
    my $l = $mpos-1; ### !!! SAM has a 1-offset hence the subtraction ADJUST IF CHANGING

    $previdx = $posidx;
    my $found = 0;
    while( $posidx < $possize) { ## assumes the positions are sorted
      #print "DEBUG: $posref->[$posidx]->[0] $posref->[$posidx]->[1] $l\n";
      if($posref->[$posidx]->[0] == $l || $posref->[$posidx]->[1] == $l) {
        $posref->[$posidx]->[5]++;
        $found = 1;
        last;
      }
      
      $posidx++;
    }
    $posidx = $previdx if $found == 0; ### sanity check, in case there is a bad read in the file, don't throw everything away
  }
  close(S);
  
  $lastchr = "NA";
  
  print "$samplename: mapping reads on the negative strand...\n";
  open(S,"samtools view -F 0x4 -f 0x10 $samplename.sorted.bam |") or die "Cannot read $samplename.sorted.bam! $!";
  while(<S>) { ### all mapped in the negative direction
    chomp;
    
    my (undef,undef,$chr,$mpos,undef,undef,undef,undef,undef,$sequence) = split /\t/;
    
    unless(defined($fragments{$chr})) {
      warn "$chr appeared in the mapping but not in the fragments directory";
      next; ## ignore for now
    }
    
    if($lastchr ne $chr) {
      $lastchr = $chr;
      $posref = $fragments{$chr};
      $posidx = 0;
      $previdx= 0;
      $possize = scalar(@{$posref});
    }
    
    my $l = $mpos+length($sequence)-1; ### !!! SAM has a 1-offset hence the subtraction ADJUST IF CHANGING

    $previdx = $posidx;
    my $found = 0;
    while( $posidx < $possize) { ## assumes the positions are sorted
      #print "DEBUG: $posref->[$posidx]->[0] $posref->[$posidx]->[1] $l\n";
      if($posref->[$posidx]->[0] == $l || $posref->[$posidx]->[1] == $l) {
        $posref->[$posidx]->[4]++;
        $found = 1;
        last;
      }
      
      $posidx++;
    }
    $posidx = $previdx if $found == 0; ### sanity check, in case there is a bad read in the file, don't throw everything away
  }
  close(S);
  
  print "$samplename: writing output...\n";
  
  open(R,">","wigfiles/$samplename.raw.wig") or die "Cannot write raw: $!";
  open(F,">","wigfiles/$samplename.filtered.wig") or die "Cannot write filetered: $!";

  for my $chr (keys(%fragments)) {
    my @rawstrings;
    my @filteredstrings;
    
    
    push @rawstrings, "variableStep chrom=$chr\n";
    push @filteredstrings, "variableStep chrom=$chr\n";

    for my $arrref (@{$fragments{$chr}}) {
      my @arr = @{$arrref};
      push @rawstrings, "$arr[0]\t$arr[4]\n" unless $arr[4] == 0;
      push @rawstrings, "$arr[1]\t$arr[5]\n" unless $arr[5] == 0;
      
      push @filteredstrings, "$arr[0]\t$arr[4]\n" unless ( $arr[4] == 0 && $arr[2] eq "NA" );
      push @filteredstrings, "$arr[1]\t$arr[5]\n" unless ( $arr[5] == 0 && $arr[3] eq "NA" );
    }
    
    if(scalar(@rawstrings)>1) {
      print R for(@rawstrings);
    }
    
    if(scalar(@filteredstrings)>1) {
      print F for(@filteredstrings);
    }
  }
  
  close(R);
  close(F);
  
  ### find information about the viewpoint

  $vchr = "chr$vchr" unless $vchr =~ /^chr/;
  
  my @positions = @{$fragments{$vchr}};
  my @left;
  my @viewpoint;
  my @right;
  my $idx;
  
  for( $idx = 0; $idx < $#positions; $idx++ ) {
    my @d = @{$positions[$idx]};
    
    last if($d[0] == $vstart || $d[1] == $vend);
  }
  
  if($idx == $#positions) {
    warn "The provided viewpoint was not found for $samplename!";
    $stats{$samplename} = ["NA","NA"];
  } else {
    #my @left = @{$positions[$idx-1]};
    #my @right = @{$positions[$idx+1]};
    my @d = @{$positions[$idx]};
    
    if($vstart == $viewpoint[0]) {
      $stats{$samplename} = [$d[4],$d[5]];
    } else {
      $stats{$samplename} = [$d[5],$d[4]];

    }
  }
}

open(S,">","stats/self-and-noncut-freq.txt") or die "Cannot write to stats file: $!";
print S "Sample\t$Non-cut\tSelf-ligation\n";
for(keys(%stats)) {
  my @a = @{$stats{$_}};
  
  print S "$_\t$a[0]\t$a[1]\n";
}
close(S);

