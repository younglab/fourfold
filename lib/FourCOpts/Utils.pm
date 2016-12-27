package FourCOpts::Utils;

use strict;

use Exporter qw(import);

our @EXPORT_OK = qw(issamplein);

sub issamplein {
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

1;
