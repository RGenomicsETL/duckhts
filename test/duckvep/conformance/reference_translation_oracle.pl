#!/usr/bin/env perl
use strict;
use warnings;
use JSON;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::TranscriptVariation;

# Actual single-exon source objects, not cached or overridden peptide methods.
# Slice attributes normally require a database adaptor; supply only the declared
# input table attribute. No translation, sequence or container method is replaced.
{
    package DuckHTS::TranslationSlice;
    use parent 'Bio::EnsEMBL::Slice';
    sub get_all_Attributes {
        my ($self, $code) = @_;
        return [] if defined($code) && lc($code) ne 'codon_table';
        return [$self->{_input_codon_table}];
    }
}
# Observe both core reference translation and the Haplosaurus container's
# reference/alternate paths. The original end-to-end observer is unchanged.
@ARGV == 1 or die "usage: reference_translation_oracle.pl cases.jsonl\n";
open(my $input, '<', $ARGV[0]) or die "Cannot read cases: $!\n";
my $json = JSON->new->canonical;
while (my $line = <$input>) {
    my $case = $json->decode($line);
    my $cds = $case->{cds};
    die "Require at least one complete ACGTN codon\n"
        unless length($cds) >= 3 && $cds =~ /^[ACGTNacgtn]+$/;
    my $slice = DuckHTS::TranslationSlice->new(-SEQ => $cds, -START => 1,
        -END => length($cds), -STRAND => 1, -SEQ_REGION_NAME => 'translation_witness',
        -COORD_SYSTEM => Bio::EnsEMBL::CoordSystem->new(-NAME => 'chromosome', -RANK => 1));
    $slice->{_input_codon_table} = Bio::EnsEMBL::Attribute->new(
        -CODE => 'codon_table', -VALUE => $case->{table});
    my $exon = Bio::EnsEMBL::Exon->new(-START => 1, -END => length($cds),
        -STRAND => 1, -PHASE => 0, -END_PHASE => length($cds) % 3, -SLICE => $slice);
    my $transcript = Bio::EnsEMBL::Transcript->new(-STABLE_ID => $case->{id},
        -SLICE => $slice, -STRAND => 1);
    $transcript->add_Exon($exon);
    my $translation = Bio::EnsEMBL::Translation->new(-STABLE_ID => $case->{id} . '_protein',
        -SEQ_START => 1, -SEQ_END => length($cds), -START_EXON => $exon, -END_EXON => $exon);
    $transcript->translation($translation);
    $transcript->edits_enabled(1);
    for my $edit (@{$case->{edits}}) {
        $translation->add_Attributes(Bio::EnsEMBL::Attribute->new(-CODE => $edit->{code},
            -VALUE => join(' ', $edit->{position1}, $edit->{position1}, $edit->{alternate})));
    }
    my $sample = Bio::EnsEMBL::Variation::Sample->new(-NAME => 'translation_sample');
    my $container = Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
        -TRANSCRIPT => $transcript, -GENOTYPES => [], -SAMPLES => [$sample]);
    my $mutated = $container->_mutate_sequences([], $sample->name);
    die "Unexpected no-edit lane count\n" unless @$mutated == 2;
    my $prepared = $transcript->translateable_seq;
    die "Source CDS changed for $case->{id}: '$cds' became '$prepared'\n"
        unless uc($prepared) eq uc($cds);
    my @independent;
    for my $variant (@{$case->{variants} || []}) {
        die "Independent witness requires a matching literal SNV\n" unless
            length($variant->{reference}) == 1 && length($variant->{alternate}) == 1 &&
            uc(substr($cds, $variant->{position1} - 1, 1)) eq uc($variant->{reference});
        my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
            -start => $variant->{position1}, -end => $variant->{position1}, -strand => 1,
            -slice => $slice, -allele_string => $variant->{reference} . '/' . $variant->{alternate},
            -variation_name => $variant->{id});
        my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
            -variation_feature => $vf, -transcript => $transcript);
        my $alleles = $tv->get_all_alternate_TranscriptVariationAlleles;
        die "No independent ALT\n" unless @$alleles;
        for my $allele (@$alleles) {
            my @terms = map { $_->SO_term } @{$allele->get_all_OverlapConsequences};
            push @independent, {id => $variant->{id}, allele => $allele->variation_feature_seq,
                consequences => \@terms, hgvsp => $allele->hgvs_protein};
        }
    }
    print $json->encode({id => $case->{id}, prepared_cds => $prepared,
        core_reference => $translation->seq, reference => $transcript->{protein},
        alternate_full => $container->_get_translation($prepared, $case->{table}),
        alternate => $mutated->[0]->{protein},
        (exists($case->{variants}) ? (independent_hgvs => \@independent) : ())}), "\n";
}
close($input) or die "Cannot close cases: $!\n";
