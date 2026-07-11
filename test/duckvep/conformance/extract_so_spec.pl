#!/usr/bin/env perl
# Extract VEP's %OVERLAP_CONSEQUENCES table without loading the Ensembl Perl
# dependency graph. Empty fields use the explicit TSV sentinel \N so generated
# files have no ambiguous trailing tabs.
use strict;
use warnings;
use Digest::SHA qw(sha256_hex);
use Getopt::Long qw(GetOptions);

my ($output, $check, $expected_sha256);
GetOptions(
    'output=s' => \$output,
    'check=s' => \$check,
    'sha256=s' => \$expected_sha256,
) or die "usage: $0 [--output file | --check file] [--sha256 hex] Constants.pm\n";
die "--output and --check are mutually exclusive\n"
    if defined $output && defined $check;
my $path = shift @ARGV
    // die "usage: $0 [--output file | --check file] [--sha256 hex] Constants.pm\n";
die "unexpected argument: $ARGV[0]\n" if @ARGV;

open my $fh, '<', $path or die "cannot open $path: $!";
binmode $fh;
local $/;
my $src = <$fh>;
close $fh or die "cannot close $path: $!";
if (defined $expected_sha256) {
    my $actual = sha256_hex($src);
    die "$path sha256 mismatch: expected $expected_sha256, got $actual\n"
        if lc($actual) ne lc($expected_sha256);
}

my @cols = qw(SO_term SO_accession impact rank tier feature_SO_term predicate);
my @rows;
while ($src =~ /->new_fast\(\{(.*?)\}\s*\)/sg) {
    my $block = $1;
    my %field;
    for my $column (@cols) {
        $field{$column} =
            ($block =~ /'$column'\s*=>\s*'([^']*)'/) ? $1 : '';
    }
    next unless $field{SO_term};
    $field{predicate} =~ s/.*:://;
    push @rows, \%field;
}
die "no VEP consequences extracted from $path\n" unless @rows;

my $text = join("\t", @cols) . "\n";
for my $row (sort { ($a->{rank} || 999) <=> ($b->{rank} || 999) } @rows) {
    $text .= join("\t", map {
        length($row->{$_}) ? $row->{$_} : '\\N'
    } @cols) . "\n";
}

if (defined $check) {
    open my $actual_fh, '<', $check or die "cannot open $check: $!";
    local $/;
    my $actual = <$actual_fh>;
    close $actual_fh or die "cannot close $check: $!";
    die "$check is stale; regenerate with make duckvep-so-spec\n"
        if $actual ne $text;
} elsif (defined $output) {
    open my $output_fh, '>', $output or die "cannot open $output: $!";
    print {$output_fh} $text;
    close $output_fh or die "cannot close $output: $!";
} else {
    print $text;
}
