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
use Bio::EnsEMBL::Variation::SampleGenotypeFeature;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::TranscriptVariation;

# Actual single-exon source objects, not cached or overridden peptide methods.
# Slice attributes normally require a database adaptor; supply only the declared
# input table attribute. This shim changes no sequence or translation method;
# VCF mode uses the standard upstream FASTA transport below.
{
    package DuckHTS::TranslationSlice;
    use parent 'Bio::EnsEMBL::Slice';
    sub get_all_Attributes {
        my ($self, $code) = @_;
        return [] if defined($code) && lc($code) ne 'codon_table';
        return [$self->{_input_codon_table}];
    }
}

# Indel HGVS expands genomic slices. Use VEP's own FASTA transport over the
# explicitly supplied contigs, not an attached Slice sequence or invented flank.
sub prepare_vcf_fasta {
    my ($input, $json) = @_;
    require File::Temp;
    require File::Spec;
    require Bio::EnsEMBL::Variation::Utils::FastaSequence;
    my $directory = File::Temp::tempdir('duckhts-translation-XXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $path = File::Spec->catfile($directory, 'reference.fa');
    open(my $fasta, '>', $path) or die "Cannot create oracle FASTA: $!\n";
    my %contigs;
    while (my $line = <$input>) {
        my $case = $json->decode($line);
        die "VCF mode requires homogeneous variant_format='vcf' cases\n" unless
            defined($case->{variant_format}) && $case->{variant_format} eq 'vcf';
        my ($id, $genome, $cds, $start) =
            @{$case}{qw(id genomic_sequence cds cds_start1)};
        die "VCF mode requires unique literal contig IDs\n" unless
            defined($id) && $id =~ /^[^\s,<>]+$/ && !$contigs{$id}++;
        die "VCF mode requires declared literal genomic and CDS sequences\n" unless
            defined($genome) && $genome =~ /^[ACGTN]+$/i &&
            defined($cds) && length($cds) >= 3 && $cds =~ /^[ACGTN]+$/i;
        die "VCF mode requires CDS matching its declared genomic span\n" unless
            defined($start) && $start =~ /^[1-9][0-9]*$/ && $start <= length($genome) &&
            length($cds) <= length($genome) - $start + 1 &&
            uc(substr($genome, $start - 1, length($cds))) eq uc($cds);
        print {$fasta} ">$id\n$genome\n" or die "Cannot write oracle FASTA: $!\n";
    }
    close($fasta) or die "Cannot close oracle FASTA: $!\n";
    seek($input, 0, 0) or die "Cannot rewind oracle cases: $!\n";
    Bio::EnsEMBL::Variation::Utils::FastaSequence::setup_fasta(
        -FASTA => $path, -OFFLINE => 1);
}

# Keep original VCF records intact until the pinned parser removes their anchors.
# In particular, a matching N anchor is legal here even though a retained N in
# the parsed allele can make the subsequent TVA peptide unavailable.
sub parse_vcf_variants {
    my ($case, $slice, $config) = @_;
    my $cds = $case->{cds};
    my $chromosome = $slice->seq_region_name;
    my $vcf = "##fileformat=VCFv4.4\n" .
        "##contig=<ID=$chromosome,length=" . $slice->seq_region_length . ">\n" .
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    my %ids;
    for my $variant (@{$case->{variants} || []}) {
        my ($id, $position, $reference, $alternate) =
            @{$variant}{qw(id position1 reference alternate)};
        die "VCF witness requires a unique literal source ID\n" unless
            defined($id) && $id ne '.' && $id =~ /^[^\s;]+$/ && !$ids{$id}++;
        die "VCF witness requires one literal ACGTN REF and ALT\n" unless
            defined($reference) && $reference =~ /^[ACGTN]+$/i &&
            defined($alternate) && $alternate =~ /^[ACGTN]+$/i;
        die "VCF witness requires a complete matching uploaded REF\n" unless
            defined($position) && $position =~ /^[1-9][0-9]*$/ &&
            $position <= length($cds) && length($reference) <= length($cds) - $position + 1 &&
            uc(substr($cds, $position - 1, length($reference))) eq uc($reference);
        my $genomic_position = $case->{cds_start1} - 1 + $position;
        $vcf .= join("\t", $chromosome, $genomic_position, $id, $reference, $alternate,
            '.', 'PASS', '.') . "\n";
    }
    open(my $records, '<', \$vcf) or die "Cannot open in-memory VCF: $!\n";
    my $parser = Bio::EnsEMBL::VEP::Parser::VCF->new({
        config => $config, file => $records, valid_chromosomes => [$chromosome]});
    my @features;
    for my $variant (@{$case->{variants} || []}) {
        my $vf = $parser->next;
        die "VCF parser lost or reassigned source $variant->{id}\n" unless
            $vf && $vf->isa('Bio::EnsEMBL::Variation::VariationFeature') &&
            $vf->variation_name eq $variant->{id} &&
            $vf->{nontrimmed_allele_string} eq
                $variant->{reference} . '/' . $variant->{alternate};
        $vf->slice($slice);
        push @features, $vf;
    }
    die "VCF parser produced an extra source\n" if $parser->next;
    close($records) or die "Cannot close in-memory VCF: $!\n";
    return \@features;
}

# Observe core reference translation, Haplosaurus reference/alternate paths,
# and independent TVA consequences/HGVS from the declared source records.
my $raw_records = @ARGV && $ARGV[0] eq '--raw-records' ? !!shift(@ARGV) : 0;
@ARGV == 1 or die "usage: reference_translation_oracle.pl [--raw-records] cases.jsonl\n";
open(my $input, '<', $ARGV[0]) or die "Cannot read cases: $!\n";
my $json = JSON->new->canonical;
my $vcf_config;
my $first_line = <$input>;
die "Oracle requires at least one case\n" unless defined($first_line);
my $vcf_mode = exists($json->decode($first_line)->{variant_format});
seek($input, 0, 0) or die "Cannot rewind oracle cases: $!\n";
prepare_vcf_fasta($input, $json) if $vcf_mode;
while (my $line = <$input>) {
    my $case = $json->decode($line);
    die "Unknown variant_format\n" if exists($case->{variant_format}) &&
        (!defined($case->{variant_format}) || $case->{variant_format} ne 'vcf');
    my $is_vcf = exists($case->{variant_format});
    die "Cannot mix VCF and default SNV cases\n" if $is_vcf != $vcf_mode;
    die "Raw-record observations require original VCF cases\n" if $raw_records && !$is_vcf;
    if ($is_vcf && !$vcf_config) {
        require Bio::EnsEMBL::VEP::Config;
        require Bio::EnsEMBL::VEP::Parser::VCF;
        $vcf_config = Bio::EnsEMBL::VEP::Config->new({
            offline => 1, database => 0, minimal => 0, check_ref => 0, lookup_ref => 0,
            dir => File::Spec->devnull, warning_file => 'STDERR', quiet => 1});
    }
    my $cds = $case->{cds};
    die "Require at least one complete ACGTN codon\n"
        unless length($cds) >= 3 && $cds =~ /^[ACGTNacgtn]+$/;
    my $cds_start = $is_vcf ? $case->{cds_start1} : 1;
    my $slice = DuckHTS::TranslationSlice->new(
        ($is_vcf ? () : (-SEQ => $cds)), -START => 1,
        -END => $is_vcf ? length($case->{genomic_sequence}) : length($cds), -STRAND => 1,
        -SEQ_REGION_NAME => $is_vcf ? $case->{id} : 'translation_witness',
        -COORD_SYSTEM => Bio::EnsEMBL::CoordSystem->new(-NAME => 'chromosome', -RANK => 1));
    $slice->{_input_codon_table} = Bio::EnsEMBL::Attribute->new(
        -CODE => 'codon_table', -VALUE => $case->{table});
    my $exon = Bio::EnsEMBL::Exon->new(-START => $cds_start, -END => $cds_start + length($cds) - 1,
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
    if ($raw_records) {
        # Fixed raw GT 1|1 gives two literal ALT slots. Construct the same VF
        # and SGF fields as Haplo::AnnotationType::Transcript, without ordinary
        # VEP minimisation. This observes container mechanics, not raw-GT parsing.
        my @observations;
        my %ids;
        for my $variant (@{$case->{variants} || []}) {
            my ($id, $position, $reference, $alternate) =
                @{$variant}{qw(id position1 reference alternate)};
            die "Raw witness requires unique source IDs and matching literal alleles\n" unless
                defined($id) && $id =~ /^[^\s;]+$/ && !$ids{$id}++ &&
                defined($position) && $position =~ /^[1-9][0-9]*$/ &&
                defined($reference) && $reference =~ /^[ACGTN]+$/i &&
                defined($alternate) && $alternate =~ /^[ACGTN]+$/i &&
                $position <= length($cds) && length($reference) <= length($cds) - $position + 1 &&
                uc(substr($cds, $position - 1, length($reference))) eq uc($reference);
            my $start = $cds_start - 1 + $position;
            my $end = $start + length($reference) - 1;
            my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
                start => $start, end => $end, strand => 1, map_weight => 1,
                allele_string => "$reference/$alternate", variation_name => $id,
                chr => $slice->seq_region_name, slice => $slice,
            });
            my $sample = Bio::EnsEMBL::Variation::Sample->new(-NAME => 'translation_sample');
            my $gt = Bio::EnsEMBL::Variation::SampleGenotypeFeature->new_fast({
                variation_feature => $vf, sample => $sample,
                genotype => [$alternate, $alternate], phased => 1,
                start => $start, end => $end, strand => 1, slice => $slice,
            });
            my $container = Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
                -TRANSCRIPT => $transcript, -GENOTYPES => [$gt], -SAMPLES => [$sample]);
            die "No raw container for source $id\n" unless $container;
            my $genotypes = $container->_filter_and_sort_genotypes(
                $container->get_all_SampleGenotypeFeatures);
            my $mutated = $container->_mutate_sequences($genotypes, $sample->name);
            my $mapping = $vf->{_cds_mapping};
            my @lanes;
            for my $lane (0..$#$mutated) {
                my $value = $mutated->[$lane];
                push @lanes, {
                    lane1 => $lane + 1, sample => $sample->name,
                    cds => $value->{cds}, protein => $value->{protein},
                    flags => {%{$value->{flags}}},
                    applied_sources => [sort map {$_->variation_name} values %{$value->{vfs}}],
                };
            }
            # Serialize before the next event reuses this transcript's caches.
            my $container_json = JSON->new->canonical->convert_blessed->encode($container);
            push @observations, {
                id => $id, source_reference => $reference, source_alternate => $alternate,
                source_position1 => $start, raw_gt => '1|1',
                mapping_start => $mapping ? $mapping->start : undef,
                mapping_end => $mapping ? $mapping->end : undef,
                reference_cds => $transcript->{cds}, reference_protein => $transcript->{protein},
                lanes => \@lanes, container => $json->decode($container_json),
            };
        }
        print $json->encode({id => $case->{id}, raw_haplotypes => \@observations}), "\n";
        next;
    }
    my $sample = Bio::EnsEMBL::Variation::Sample->new(-NAME => 'translation_sample');
    my $container = Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
        -TRANSCRIPT => $transcript, -GENOTYPES => [], -SAMPLES => [$sample]);
    my $mutated = $container->_mutate_sequences([], $sample->name);
    die "Unexpected no-edit lane count\n" unless @$mutated == 2;
    my $prepared = $transcript->translateable_seq;
    die "Source CDS changed for $case->{id}: '$cds' became '$prepared'\n"
        unless uc($prepared) eq uc($cds);
    my @parsed = $is_vcf ? @{parse_vcf_variants($case, $slice, $vcf_config)} : ();
    my @independent;
    for my $variant (@{$case->{variants} || []}) {
        my $vf;
        my %source;
        if ($is_vcf) {
            $vf = shift @parsed;
            # Parser coordinates are genomic, one-based inclusive. Insertions
            # have start=end+1; the source position1 remains CDS-local.
            %source = (source_reference => $variant->{reference},
                source_alternate => $variant->{alternate}, parser_start => 0 + $vf->start,
                parser_end => 0 + $vf->end, parser_allele_string => $vf->allele_string);
        } else {
            die "Independent witness requires a matching literal SNV\n" unless
                length($variant->{reference}) == 1 && length($variant->{alternate}) == 1 &&
                uc(substr($cds, $variant->{position1} - 1, 1)) eq uc($variant->{reference});
            $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
                -start => $variant->{position1}, -end => $variant->{position1}, -strand => 1,
                -slice => $slice, -allele_string => $variant->{reference} . '/' . $variant->{alternate},
                -variation_name => $variant->{id});
        }
        my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
            -variation_feature => $vf, -transcript => $transcript);
        my $alleles = $tv->get_all_alternate_TranscriptVariationAlleles;
        die "No independent ALT\n" unless @$alleles;
        die "VCF witness did not retain exactly one ALT\n" if $is_vcf && @$alleles != 1;
        for my $allele (@$alleles) {
            my @terms = map { $_->SO_term } @{$allele->get_all_OverlapConsequences};
            push @independent, {id => $variant->{id}, allele => $allele->variation_feature_seq,
                consequences => \@terms, hgvsp => $allele->hgvs_protein, %source};
        }
    }
    print $json->encode({id => $case->{id}, prepared_cds => $prepared,
        core_reference => $translation->seq, reference => $transcript->{protein},
        alternate_full => $container->_get_translation($prepared, $case->{table}),
        alternate => $mutated->[0]->{protein},
        (exists($case->{variants}) ? (independent_hgvs => \@independent) : ())}), "\n";
}
close($input) or die "Cannot close cases: $!\n";
