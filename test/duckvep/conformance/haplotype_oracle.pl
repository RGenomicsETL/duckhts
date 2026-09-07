#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::VEP::Haplo::Runner;
use JSON;
use File::Basename qw(dirname);
my $phase_output;

# Observe the real release-116 parser, mapper and container. Only output is
# replaced: upstream's tabular output omits the complete sequences, and its
# JSON prefetch does not include CDS frameshift flags. No genotype, projection,
# mutation, translation or phase behavior is overridden here.
{
    package DuckHTS::HaploObserver;
    use parent 'Bio::EnsEMBL::VEP::Haplo::Runner';

    sub dump_TranscriptHaplotypeContainer {
        my ($self, $container) = @_;
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
        if ($phase_output) {
            print {$phase_output} JSON->new->canonical->encode({
                transcript => $container->transcript->stable_id,
                default_ploidy => $container->_default_ploidy,
                sample_ploidy => $container->_sample_ploidy,
                calls => [map {{
                    source_id => $_->variation_feature->variation_name,
                    sample => $_->sample->name,
                    genotype => $_->genotype,
                }} @{$container->get_all_SampleGenotypeFeatures}],
            }), "\n";
        }
        $self->{_output_lines_count}++;
    }
}

(@ARGV == 3 || @ARGV == 4) or die
    "usage: haplotype_oracle.pl input.vcf reference.fa model.gff3.gz [phase-observations.jsonl]\n";
my ($vcf, $fasta, $gff, $phase_path) = @ARGV;
open($phase_output, '>', $phase_path) or die "cannot write $phase_path: $!" if defined($phase_path);
my $runner = DuckHTS::HaploObserver->new({
    input_file => $vcf,
    fasta => $fasta,
    gff => $gff,
    dir => dirname($vcf),
    output_file => 'STDOUT',
    warning_file => $vcf . '.warnings',
    database => 0,
    no_stats => 1,
});
$runner->run;
close($phase_output) or die "cannot close phase observations: $!" if $phase_output;
