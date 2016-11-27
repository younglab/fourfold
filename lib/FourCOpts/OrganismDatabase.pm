package FourCOpts::OrganismDatabase;

use strict;

use Exporter qw(import);

our @EXPORT_OK = qw(loadorgdatabase);

sub loadorgdatabase {
  die "Not enough arguments" if(scalar(@_)<1);
  
  my ($organismdatabase) = @_;
  
  my %organisms;
  
  open(D,"<","$organismdatabase") or die "Cannot read $organismdatabase: $!";

  while(<D>) {
    chomp;
  
    my ($id,$bowtie,$fasta,$chromsizes) = split /\t/;
  
    next if $id =~ /^$/;
  
    die "In organism database, cannot find bowtie index for $id!" unless -e "$bowtie.1.ebwt";
    die "In organism database, cannot find FASTA file for $id!" unless -e $fasta;
    die "In organism database, cannot find chromosome sizes for $id!" unless -e $chromsizes;
  
    for(split(/,/,$id)) {
      $organisms{lc $_} = [$bowtie,$fasta,$chromsizes];
    }
  }
  close(D);
  
  return \%organisms;
}

1;
