#!/usr/bin/env perl
use strict;
use warnings;

use Text::ParseWords qw(parse_line);

my $require_current = 0;
if (@ARGV && $ARGV[0] eq '--require-current') {
    shift @ARGV;
    $require_current = 1;
}

my $path = shift //
    'test/duckvep/conformance/data/state_machine_transitions.tsv';
@ARGV == 0
    or die "usage: $0 [--require-current] [state_machine_transitions.tsv]\n";

open my $fh, '<', $path or die "cannot open $path: $!\n";
my $header = <$fh>;
defined $header or die "$path is empty\n";
chomp $header;
my @expected = qw(
    transition_id input_states to_state paths observed_dimensions implementation
    vep_authority fixed_test randomized_property executable_campaign terminal_outcomes
    proof_status
);
my @columns = split /\t/, $header, -1;
join("\t", @columns) eq join("\t", @expected)
    or die "$path has an unexpected header\n";

my %valid_status = map { $_ => 1 } qw(
    finite_basis_proved property_equivalent executable_sampled mechanics_only
    not_implemented
);
my %path_spec = (
    small                 => ['annotation_rows', 'executable'],
    small_nmd             => ['annotation_rows', 'executable'],
    small_hgvs            => ['annotation_rows', 'executable'],
    small_nmd_hgvs        => ['annotation_rows', 'executable'],
    structural            => ['annotation_rows', 'executable'],
    interval              => ['annotation_rows', 'executable'],
    haplotype_mechanics   => ['haplotype_coding_context', 'executable'],
    haplotype_future      => ['haplotype_consequence_set', 'planned'],
);
my (
    %seen_transition, %produced_state, %consumed_state, %status_count, %path_seen,
    %used_campaign, %used_outcome
);
my @rows;
my $line_number = 1;
while (my $line = <$fh>) {
    $line_number++;
    chomp $line;
    next if $line eq '';
    my @values = split /\t/, $line, -1;
    @values == @columns
        or die "$path:$line_number has " . scalar(@values) .
            " fields; expected " . scalar(@columns) . "\n";
    my %row;
    @row{@columns} = @values;
    for my $column (@columns) {
        $row{$column} ne ''
            or die "$path:$line_number has empty $column\n";
    }
    $row{transition_id} =~ /^[a-z][a-z0-9_]*$/
        or die "$path:$line_number has invalid transition_id\n";
    $row{to_state} =~ /^[a-z][a-z0-9_]*$/
        or die "$path:$line_number has invalid to_state\n";
    my @input_states = split /\+/, $row{input_states}, -1;
    @input_states && !grep { $_ !~ /^[a-z][a-z0-9_]*$/ } @input_states
        or die "$path:$line_number has invalid input_states\n";
    my %input_seen;
    !grep { $input_seen{$_}++ } @input_states
        or die "$path:$line_number duplicates an input state\n";
    my @paths = split /;/, $row{paths}, -1;
    @paths && !grep { $_ !~ /^[a-z][a-z0-9_]*$/ || !exists $path_spec{$_} } @paths
        or die "$path:$line_number has invalid paths\n";
    my %path_in_row;
    !grep { $path_in_row{$_}++ } @paths
        or die "$path:$line_number duplicates a path\n";
    !$seen_transition{$row{transition_id}}++
        or die "$path:$line_number duplicates $row{transition_id}\n";
    $valid_status{$row{proof_status}}
        or die "$path:$line_number has invalid proof_status\n";
    if ($row{proof_status} eq 'not_implemented') {
        $row{implementation} eq '-' && $row{fixed_test} eq '-' &&
            $row{randomized_property} eq '-' &&
            $row{executable_campaign} eq '-'
            or die "$path:$line_number claims evidence for an unimplemented transition\n";
    } else {
        for my $implementation (split /;/, $row{implementation}) {
            -f $implementation
                or die "$path:$line_number missing implementation $implementation\n";
        }
    }
    if ($row{proof_status} eq 'finite_basis_proved') {
        $row{fixed_test} ne '-'
            or die "$path:$line_number finite proof lacks an executable test\n";
    } elsif ($row{proof_status} eq 'property_equivalent') {
        $row{fixed_test} ne '-' && $row{randomized_property} ne '-'
            or die "$path:$line_number property proof lacks fixed/randomized evidence\n";
    } elsif ($row{proof_status} eq 'executable_sampled') {
        $row{fixed_test} ne '-' && $row{executable_campaign} ne '-'
            or die "$path:$line_number sampled proof lacks fixed/executable evidence\n";
    } elsif ($row{proof_status} eq 'mechanics_only') {
        $row{fixed_test} ne '-' || $row{randomized_property} ne '-'
            or die "$path:$line_number mechanics claim lacks executable evidence\n";
    }
    if ($row{executable_campaign} ne '-') {
        for my $campaign (split /;/, $row{executable_campaign}) {
            $campaign =~ /^[a-z][a-z0-9_]*$/
                or die "$path:$line_number has invalid campaign id $campaign\n";
            $used_campaign{$campaign}++;
        }
    }
    if ($row{terminal_outcomes} ne '-') {
        for my $outcome (split /;/, $row{terminal_outcomes}) {
            $outcome =~ /^[a-z][a-z0-9_]*$/
                or die "$path:$line_number has invalid terminal outcome $outcome\n";
            $used_outcome{$outcome}++;
        }
    }
    $produced_state{$row{to_state}} = 1;
    $consumed_state{$_}++ for @input_states;
    $path_seen{$_}++ for @paths;
    $status_count{$row{proof_status}}++;
    push @rows, \%row;
}
close $fh or die "cannot close $path: $!\n";

@rows or die "$path has no transitions\n";
for my $row (@rows) {
    for my $input_state (split /\+/, $row->{input_states}) {
        next if $input_state eq 'raw_allele';
        $produced_state{$input_state}
            or die "$row->{transition_id} consumes unproduced state $input_state\n";
    }
}
my %terminal_state = map { $_ => 1 } qw(
    annotation_rows haplotype_consequence_set
);
for my $state (sort keys %produced_state) {
    next if $terminal_state{$state};
    $consumed_state{$state}
        or die "$path produces disconnected state $state\n";
}

my %path_class_count;
for my $execution_path (sort keys %path_spec) {
    $path_seen{$execution_path}
        or die "$path has no transitions for execution path $execution_path\n";
    my ($terminal, $path_class) = @{$path_spec{$execution_path}};
    my %reachable = (raw_allele => 1);
    my $changed = 1;
    while ($changed) {
        $changed = 0;
        for my $row (@rows) {
            my %row_paths = map { $_ => 1 } split /;/, $row->{paths};
            next unless $row_paths{$execution_path};
            if ($path_class eq 'executable' &&
                $row->{proof_status} eq 'not_implemented') {
                die "$execution_path includes unimplemented transition " .
                    "$row->{transition_id}\n";
            }
            my @inputs = split /\+/, $row->{input_states};
            next if grep { !$reachable{$_} } @inputs;
            if (!$reachable{$row->{to_state}}) {
                $reachable{$row->{to_state}} = 1;
                $changed = 1;
            }
        }
    }
    for my $row (@rows) {
        my %row_paths = map { $_ => 1 } split /;/, $row->{paths};
        next unless $row_paths{$execution_path};
        my @missing = grep { !$reachable{$_} } split /\+/, $row->{input_states};
        @missing == 0
            or die "$row->{transition_id} is unsatisfiable on $execution_path: " .
                join(', ', @missing) . " cannot be reached\n";
    }
    $reachable{$terminal}
        or die "$execution_path cannot reach terminal state " .
            "$terminal\n";
    $path_class_count{$path_class}++;
}

my $test_source = 'test/duckvep/property/duckvep_kernel_prop.c';
open my $tfh, '<', $test_source or die "cannot open $test_source: $!\n";
my $tests;
{
    local $/;
    $tests = <$tfh>;
}
close $tfh or die "cannot close $test_source: $!\n";
for my $row (@rows) {
    next if $row->{fixed_test} eq '-';
    my $name = quotemeta($row->{fixed_test});
    $tests =~ /TEST\s+$name\s*\(/
        or die "$row->{transition_id} names absent fixed test $row->{fixed_test}\n";
}

my %property_name;
while ($tests =~ /(?:cfg|config)\.name\s*=\s*"([^"]+)"/sg) {
    $property_name{$1} = 1;
}
for my $row (@rows) {
    next if $row->{randomized_property} eq '-';
    for my $property (split /;/, $row->{randomized_property}) {
        $property_name{$property}
            or die "$row->{transition_id} names absent randomized property " .
                "$property\n";
    }
}

my $campaign_path =
    'test/duckvep/conformance/data/state_machine_campaigns.tsv';
open my $cfh, '<', $campaign_path
    or die "cannot open $campaign_path: $!\n";
my $campaign_header = <$cfh>;
defined $campaign_header or die "$campaign_path is empty\n";
chomp $campaign_header;
$campaign_header eq join("\t", qw(
    campaign_id evidence_path source_revision artifact_md5 corpus model
    oracle_version oracle_build_token authority
))
    or die "$campaign_path has an unexpected header\n";
my %campaign;
my $campaign_line = 1;
my %evidence_cache;
my $head;
if ($require_current) {
    open my $git, '-|', 'git', 'rev-parse', 'HEAD'
        or die "cannot start git rev-parse: $!\n";
    $head = <$git>;
    close $git or die "cannot determine current Git revision\n";
    chomp $head;
}
while (my $line = <$cfh>) {
    $campaign_line++;
    chomp $line;
    next if $line eq '';
    my ($id, $evidence_path, $source_revision, $artifact_md5, $corpus,
        $model, $oracle_version, $oracle_build_token, $authority, @extra) =
        split /\t/, $line, -1;
    !@extra && defined $authority &&
        !grep { !defined $_ || $_ eq '' } (
            $id, $evidence_path, $source_revision, $artifact_md5, $corpus,
            $model, $oracle_version, $oracle_build_token, $authority
        )
        or die "$campaign_path:$campaign_line has an invalid row\n";
    $id =~ /^[a-z][a-z0-9_]*$/
        or die "$campaign_path:$campaign_line has invalid campaign id\n";
    !$campaign{$id}
        or die "$campaign_path:$campaign_line duplicates campaign $id\n";
    -f $evidence_path
        or die "$campaign_path:$campaign_line missing evidence $evidence_path\n";
    $source_revision =~ /^[0-9a-f]{40}$/
        or die "$campaign_path:$campaign_line has invalid source revision\n";
    $artifact_md5 =~ /^[0-9a-f]{32}$/
        or die "$campaign_path:$campaign_line has invalid artifact MD5\n";
    open my $commit, '-|', 'git', 'cat-file', '-e', "$source_revision^{commit}"
        or die "cannot start git cat-file for $source_revision: $!\n";
    close $commit
        or die "$campaign_path:$campaign_line source revision is absent from Git\n";
    if ($require_current && $source_revision ne $head) {
        die "$campaign_path:$campaign_line campaign $id was measured at " .
            "$source_revision, not current HEAD $head\n";
    }

    if (!exists $evidence_cache{$evidence_path}) {
        open my $efh, '<', $evidence_path
            or die "cannot open $evidence_path: $!\n";
        my $evidence_header = <$efh>;
        defined $evidence_header or die "$evidence_path is empty\n";
        chomp $evidence_header;
        my @evidence_names = parse_line(',', 0, $evidence_header);
        my %column = map { $evidence_names[$_] => $_ } 0 .. $#evidence_names;
        my @rows;
        while (my $evidence_line = <$efh>) {
            chomp $evidence_line;
            next if $evidence_line eq '';
            my @value = parse_line(',', 0, $evidence_line);
            @value == @evidence_names
                or die "$evidence_path contains a malformed CSV row\n";
            push @rows, \@value;
        }
        close $efh or die "cannot close $evidence_path: $!\n";
        $evidence_cache{$evidence_path} = {
            names => \@evidence_names,
            column => \%column,
            rows => \@rows,
        };
    }
    my @evidence_names = @{ $evidence_cache{$evidence_path}->{names} };
    my %column = %{ $evidence_cache{$evidence_path}->{column} };
    for my $required (qw(source_revision artifact_md5 corpus model
                         oracle_version oracle_build)) {
        exists $column{$required}
            or die "$evidence_path lacks required column $required\n";
    }
    my @selected;
    for my $value_ref (@{ $evidence_cache{$evidence_path}->{rows} }) {
        my @value = @$value_ref;
        next unless
            $value[$column{source_revision}] eq $source_revision &&
            $value[$column{artifact_md5}] eq $artifact_md5 &&
            $value[$column{corpus}] eq $corpus &&
            $value[$column{model}] eq $model &&
            $value[$column{oracle_version}] eq $oracle_version;
        index($value[$column{oracle_build}], $oracle_build_token) >= 0
            or die "$campaign_path:$campaign_line selected row lacks " .
                "oracle-build token $oracle_build_token\n";
        push @selected, \@value;
    }
    @selected
        or die "$campaign_path:$campaign_line exact evidence key selects no rows\n";

    if (exists $column{stratum_kind}) {
        for my $required (qw(n exact_agree exact_discordant unresolved
                             resolved_discordant term_mismatch engine_extra
                             engine_missing)) {
            exists $column{$required}
                or die "$evidence_path lacks consequence counter $required\n";
        }
        my @all = grep {
            $_->[$column{stratum_kind}] eq 'all'
        } @selected;
        @all == 1
            or die "$campaign_path:$campaign_line must select one all-row\n";
        my $all = $all[0];
        $all->[$column{n}] > 0 &&
            $all->[$column{exact_agree}] == $all->[$column{n}] &&
            $all->[$column{exact_discordant}] == 0 &&
            $all->[$column{unresolved}] == 0 &&
            $all->[$column{resolved_discordant}] == 0 &&
            $all->[$column{term_mismatch}] == 0 &&
            $all->[$column{engine_extra}] == 0 &&
            $all->[$column{engine_missing}] == 0
            or die "$campaign_path:$campaign_line all-row is not exact\n";
        for my $row (@selected) {
            for my $counter (qw(exact_discordant resolved_discordant
                                term_mismatch engine_extra engine_missing)) {
                $row->[$column{$counter}] == 0
                    or die "$campaign_path:$campaign_line contains nonzero " .
                        "$counter\n";
            }
        }
    } elsif (exists $column{metric} && exists $column{comparison} &&
             exists $column{n}) {
        my %seen_metric;
        for my $row (@selected) {
            my $comparison = $row->[$column{comparison}];
            ($comparison eq 'match' || $comparison eq 'both_absent')
                or die "$campaign_path:$campaign_line contains HGVS " .
                    "comparison $comparison\n";
            $row->[$column{n}] > 0
                or die "$campaign_path:$campaign_line contains empty HGVS row\n";
            $seen_metric{$row->[$column{metric}]} = 1;
        }
        $seen_metric{hgvsc} && $seen_metric{hgvsp}
            or die "$campaign_path:$campaign_line lacks HGVSc or HGVSp evidence\n";
    } else {
        die "$campaign_path:$campaign_line has unsupported evidence schema\n";
    }
    $campaign{$id} = 1;
}
close $cfh or die "cannot close $campaign_path: $!\n";
for my $id (sort keys %used_campaign) {
    $campaign{$id} or die "$path names undefined executable campaign $id\n";
}
for my $id (sort keys %campaign) {
    $used_campaign{$id} or die "$campaign_path contains unused campaign $id\n";
}

my $outcome_path =
    'test/duckvep/conformance/data/state_machine_outcomes.tsv';
open my $ofh, '<', $outcome_path
    or die "cannot open $outcome_path: $!\n";
my $outcome_header = <$ofh>;
defined $outcome_header or die "$outcome_path is empty\n";
chomp $outcome_header;
$outcome_header eq
    "outcome_id\toutcome_class\tsemantics\timplementation\tfixed_test"
    or die "$outcome_path has an unexpected header\n";
my %outcome;
my $outcome_line = 1;
while (my $line = <$ofh>) {
    $outcome_line++;
    chomp $line;
    next if $line eq '';
    my ($id, $class, $semantics, $implementation, $fixed_test, @extra) =
        split /\t/, $line, -1;
    !@extra && defined $fixed_test &&
        $id ne '' && $class ne '' && $semantics ne '' &&
        $implementation ne '' && $fixed_test ne ''
        or die "$outcome_path:$outcome_line has an invalid row\n";
    $id =~ /^[a-z][a-z0-9_]*$/
        or die "$outcome_path:$outcome_line has invalid outcome id\n";
    $class eq 'error' || $class eq 'status'
        or die "$outcome_path:$outcome_line has invalid outcome class\n";
    !$outcome{$id}
        or die "$outcome_path:$outcome_line duplicates outcome $id\n";
    for my $source (split /;/, $implementation) {
        -f $source
            or die "$outcome_path:$outcome_line missing implementation $source\n";
    }
    my $name = quotemeta($fixed_test);
    $tests =~ /TEST\s+$name\s*\(/
        or die "$outcome_path:$outcome_line names absent fixed test $fixed_test\n";
    $outcome{$id} = 1;
}
close $ofh or die "cannot close $outcome_path: $!\n";
for my $id (sort keys %used_outcome) {
    $outcome{$id} or die "$path names undefined terminal outcome $id\n";
}
for my $id (sort keys %outcome) {
    $used_outcome{$id} or die "$outcome_path contains unused outcome $id\n";
}

exists $seen_transition{select_so_terms} or die "missing select_so_terms transition\n";
my ($selection) = grep { $_->{transition_id} eq 'select_so_terms' } @rows;
$selection->{proof_status} eq 'finite_basis_proved'
    or die "select_so_terms must retain the finite-basis proof\n";

print scalar(@rows), " state transitions: ",
    join(", ", map { "$_=$status_count{$_}" } sort keys %status_count),
    "; $path_class_count{executable} executable paths; ",
    "$path_class_count{planned} structurally connected planned path",
    $path_class_count{planned} == 1 ? "\n" : "s\n";
