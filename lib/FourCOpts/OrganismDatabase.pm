package FourCOpts::OrganismDatabase;

use strict;

use Exporter qw(import);

our @EXPORT_OK = qw(loadorgdatabase lookup_organism);

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

sub lookup_organism {

  die "Not enough arguments for lookup_organism (@_)" if scalar(@_) < 2;
  
  my ($dbref,$organism) = @_;
  
  my $o = lc $organism;
  
  die "Cannot find $organism within organism database" if !defined($dbref->{$o});
  
  return $dbref->{$o};
}

1;
