#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::VEP::Haplo::Runner;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;
use JSON;
use Scalar::Util qw(refaddr);
use File::Basename qw(dirname);
my $phase_output;
my %phase_calls;
my %lane_calls;
my $container_json = @ARGV && $ARGV[0] eq '--container-json' ? !!shift(@ARGV) : 0;

# Run the release-116 parser, mapper and container. Construction delegates to
# upstream unchanged; the observer copies shared genotype/mapping fields before
# another transcript can overwrite them. Complete sequence output and owned
# observations are serialized without changing biological operations.
{
    package DuckHTS::HaploObserver;
    use parent 'Bio::EnsEMBL::VEP::Haplo::Runner';

    sub dump_TranscriptHaplotypeContainer {
        my ($self, $container) = @_;
        if ($container_json) {
            $self->{_json} ||= JSON->new->canonical->convert_blessed;
            $self->SUPER::dump_TranscriptHaplotypeContainer($container);
        } else {
            my @haplotypes;
            foreach my $cds (@{$container->get_all_CDSHaplotypes}) {
                my $protein = $cds->get_ProteinHaplotype;
                push @haplotypes, {
                    cds => $cds->seq,
                    protein => $protein->seq,
                    flags => [sort @{$cds->get_all_flags}],
                    contributors => [sort map {$_->variation_name}
                                          @{$cds->get_all_VariationFeatures}],
                    samples => $cds->get_all_sample_counts,
                    count => $cds->count,
                };
            }
            print JSON->new->canonical->encode({
                transcript => $container->transcript->stable_id,
                total_haplotype_count => $container->total_haplotype_count,
                haplotypes => [sort {$a->{cds} cmp $b->{cds}} @haplotypes],
            }), "\n";
            $self->{_output_lines_count}++;
        }
        if ($phase_output) {
            my $calls = delete $phase_calls{Scalar::Util::refaddr($container)};
            die "missing construction-time genotype observations\n" unless defined $calls;
            my $lanes = delete $lane_calls{Scalar::Util::refaddr($container)} || [];
            print {$phase_output} JSON->new->canonical->encode({
                transcript => $container->transcript->stable_id,
                reference_cds => $container->transcript->{cds},
                ($container_json ? (reference_protein => $container->transcript->{protein}) : ()),
                default_ploidy => $container->_default_ploidy,
                sample_ploidy => $container->_sample_ploidy,
                source_buffer => [map {{
                    ids => $_->{ids}, chrom => $_->{chr},
                    start => $_->{start}, end => $_->{end}, alleles => $_->{alleles},
                }} @{$self->get_InputBuffer->buffer}],
                calls => $calls,
                replay_lanes => [sort {$a->{sample} cmp $b->{sample} || $a->{lane1} <=> $b->{lane1}} @$lanes],
            }), "\n";
        }
    }
}

(@ARGV == 3 || @ARGV == 4) or die
    "usage: haplotype_oracle.pl [--container-json] input.vcf reference.fa model.gff3.gz [phase-observations.jsonl]\n";
my ($vcf, $fasta, $gff, $phase_path) = @ARGV;
if (defined $phase_path) {
    open($phase_output, '>', $phase_path) or die "cannot write $phase_path: $!";
    # Copy each file lane before sequence grouping unions contributing sources.
    # Samples handled only by _add_reference_haplotypes do not enter this mutator.
    my $mutator = \&Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer::_mutate_sequences;
    {
        no warnings 'redefine';
        *Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer::_mutate_sequences = sub {
            my ($container, $genotypes, $sample) = @_;
            my $mutated = $mutator->(@_);
            for my $lane (0..$#$mutated) {
                my $value = $mutated->[$lane];
                push @{$lane_calls{refaddr($container)}}, {
                    sample => $sample, lane1 => $lane + 1,
                    cds => $value->{cds}, protein => $value->{protein},
                    applied_sources => [map {
                        my $vf = $value->{vfs}->{$_};
                        {
                            allele_key => $_, source_id => $vf->variation_name,
                            source_key => $vf->{_th_identifier},
                        }
                    } sort keys %{$value->{vfs}}],
                };
            }
            return $mutated;
        };
    }
    my $constructor = \&Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer::new;
    no warnings 'redefine';
    *Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer::new = sub {
        my $container = $constructor->(@_);
        if ($container) {
            $phase_calls{refaddr($container)} = [map {
                my $vf = $_->variation_feature;
                my $mapping = $vf->{_cds_mapping};
                {
                    source_id => $vf->variation_name,
                    source_key => $vf->{_th_identifier},
                    mapping_start => $mapping ? $mapping->start : undef,
                    mapping_end => $mapping ? $mapping->end : undef,
                    sample => $_->sample->name,
                    genotype => [@{$_->genotype}],
                }
            } @{$container->get_all_SampleGenotypeFeatures}];
        }
        return $container;
    };
}
my $runner = DuckHTS::HaploObserver->new({
    input_file => $vcf,
    fasta => $fasta,
    gff => $gff,
    dir => dirname($vcf),
    output_file => 'STDOUT',
    warning_file => $vcf . '.warnings',
    database => 0,
    no_stats => 1,
    json => $container_json,
});
$runner->run;
die "unemitted construction-time genotype observations\n" if keys %phase_calls;
die "unemitted mutation observations\n" if keys %lane_calls;
close($phase_output) or die "cannot close phase observations: $!" if $phase_output;
