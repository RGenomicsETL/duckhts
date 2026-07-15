package DuckVEPNMDState;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  return {
    DuckVEP_NMD_CDS =>
      'TranscriptVariation CDS start,end observed by the NMD plugin'
  };
}

sub run {
  my ($self, $tva) = @_;
  my $tv = $tva->transcript_variation;
  my $start = $tv->cds_start;
  my $end = $tv->cds_end;

  return {
    DuckVEP_NMD_CDS =>
      defined($start) && defined($end) ? "$start,$end" : 'undefined'
  };
}

1;
