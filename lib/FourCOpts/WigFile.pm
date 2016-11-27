package FourCOpts::WigFile;

use strict;

use Exporter qw(import);

our @EXPORT_OK = qw(writewigfile);

sub writewigfile {
  die "Not enough arguments" if(scalar(@_)<3);
  my ($tablefile,$wigfile,$desc) = @_;

  open(O,"<",$tablefile) or die "Cannot read $tablefile: $!";
  open(W,">",$wigfile) or die "Cannot write to $wigfile: $!";
    
  print W "track type=wiggle_0 name=\"$desc\" description=\"$desc\"\n";
    
  my $lastchr = '';
  while(<O>) {
    chomp;
      
    my ($chr,$pos,$signal) = split /\t/;
      
    unless($lastchr eq $chr) {
      print W "variableStep chrom=$chr\n";
      $lastchr = $chr;
    }
    print W "$pos\t$signal\n";
  }
  close(O);
  close(W);
}
    
1;
