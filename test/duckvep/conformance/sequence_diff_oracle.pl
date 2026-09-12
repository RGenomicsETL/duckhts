#!/usr/bin/env perl
use strict;
use warnings;
use JSON;
use Bio::EnsEMBL::Variation::TranscriptHaplotype;

# Observe the pinned pure-Perl alignment and raw-difference implementation.
# This lane does not override the existing end-to-end Haplosaurus observer.
die "This contract requires the pure-Perl NW path, not Bio::Ext::Align\n"
    if eval { require Bio::Ext::Align; 1 };
{
    package DuckHTS::SequencePair;
    use parent 'Bio::EnsEMBL::Variation::TranscriptHaplotype';
    sub reference_seq { return $_[0]->{_reference}; }
}
@ARGV == 1 or die "usage: sequence_diff_oracle.pl pairs.jsonl\n";
open(my $input, '<', $ARGV[0]) or die "Cannot read pairs: $!\n";
my $json = JSON->new->canonical;
while (my $line = <$input>) {
    my $pair = $json->decode($line);
    my $haplotype = DuckHTS::SequencePair->new(
        -type => 'protein', -seq => $pair->{alternate}, -indel => $pair->{align});
    $haplotype->{_reference} = $pair->{reference};
    my ($aligned_ref, $aligned_alt) = @{$haplotype->get_aligned_sequences};
    my @differences;
    foreach my $raw (@{$haplotype->_get_raw_diffs}) {
        my $column = $raw->{p} - length($raw->{a1}) + 1;
        my $prefix_ref = substr($aligned_ref, 0, $column);
        my $prefix_alt = substr($aligned_alt, 0, $column);
        my ($ref, $alt) = @{$raw}{qw(a1 a2)};
        $ref =~ s/-//g; $alt =~ s/-//g;
        push @differences, {ref_start0 => ($prefix_ref =~ tr/-//c),
            alt_start0 => ($prefix_alt =~ tr/-//c), reference => $ref, alternate => $alt,
            alignment_start0 => $column};
    }
    print $json->encode({id => $pair->{id}, differences => \@differences}), "\n";
}
close($input) or die "Cannot close pairs: $!\n";
