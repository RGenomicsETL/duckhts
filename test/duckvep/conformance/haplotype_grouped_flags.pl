use strict;
use warnings;
use Bio::EnsEMBL::VEP::Haplo::Runner;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;
use File::Basename qw(dirname);
use JSON;

(@ARGV == 3 || @ARGV == 4) or die "usage: grouped_flags.pl calls.vcf reference.fa model.gff3.gz [lanes.jsonl]\n";
my ($vcf, $fasta, $gff, $lane_path) = @ARGV;
my $output;
if (defined $lane_path) {
    open($output, '>', $lane_path) or die "cannot write $lane_path: $!";
    my $mutator = \&Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer::_mutate_sequences;
    no warnings 'redefine';
    *Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer::_mutate_sequences = sub {
        my ($container, $genotypes, $sample) = @_;
        my $mutated = $mutator->(@_);
        for my $lane (0..$#$mutated) {
            my $value = $mutated->[$lane];
            print {$output} JSON->new->canonical->encode({
                sample => $sample, lane1 => $lane + 1,
                cds => $value->{cds}, protein => $value->{protein},
                flags => {%{$value->{flags}}},
            }), "\n";
        }
        return $mutated;
    };
}
my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({
    input_file => $vcf, fasta => $fasta, gff => $gff, dir => dirname($vcf),
    output_file => 'STDOUT', warning_file => $vcf . '.warnings',
    database => 0, no_stats => 1, json => 1,
});
$runner->run;
close($output) or die "cannot close lane observations: $!" if $output;
