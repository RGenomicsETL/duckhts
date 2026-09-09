use strict;
use warnings;
use Bio::EnsEMBL::VEP::Haplo::Runner;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;
use File::Basename qw(dirname);
use JSON;

(@ARGV == 3 || @ARGV == 4) or die "usage: grouped_flags.pl calls.vcf reference.fa model.gff3.gz [lanes.jsonl]\n";
my ($vcf, $fasta, $gff, $lane_path) = @ARGV;
my $output;
my $groups_output;
if (defined $lane_path) {
    open($output, '>', $lane_path) or die "cannot write $lane_path: $!";
    open($groups_output, '>', "$lane_path.groups.jsonl") or die "cannot write group observations: $!";
    my $ordinal = 0;
    my $mutator = \&Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer::_mutate_sequences;
    my $dump = \&Bio::EnsEMBL::VEP::Haplo::Runner::dump_TranscriptHaplotypeContainer;
    no warnings 'redefine';
    *Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer::_mutate_sequences = sub {
        my ($container, $genotypes, $sample) = @_;
        my $mutated = $mutator->(@_);
        for my $lane (0..$#$mutated) {
            my $value = $mutated->[$lane];
            print {$output} JSON->new->canonical->encode({
                sample => $sample, lane1 => $lane + 1,
                traversal_ordinal => ++$ordinal,
                cds => $value->{cds}, protein => $value->{protein},
                flags => {%{$value->{flags}}},
            }), "\n";
        }
        return $mutated;
    };
    *Bio::EnsEMBL::VEP::Haplo::Runner::dump_TranscriptHaplotypeContainer = sub {
        my ($runner, $container) = @_;
        $dump->(@_);
        print {$groups_output} JSON->new->canonical->encode({
            transcript => $container->transcript->stable_id,
            reference_cds => $container->transcript->{cds},
            groups => [map {{cds => $_->seq,
                flags => {indel => $_->has_indel || 0,
                    frameshift => $_->{_frameshift} || 0, length_diff => $_->length_diff},
                categories => [sort @{$_->get_all_flags}],
            }} @{$container->get_all_CDSHaplotypes}],
        }), "\n";
        $ordinal = 0;
    };
}
my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({
    input_file => $vcf, fasta => $fasta, gff => $gff, dir => dirname($vcf),
    output_file => 'STDOUT', warning_file => $vcf . '.warnings',
    database => 0, no_stats => 1, json => 1,
});
$runner->run;
close($output) or die "cannot close lane observations: $!" if $output;
close($groups_output) or die "cannot close group observations: $!" if $groups_output;
