#!/usr/bin/env perl
use strict;
use warnings;

# Emit only compiled property sources after checking complete test registration.
my $directory = shift // 'test/duckvep/property';
@ARGV == 0 or die "usage: $0 [property-directory]\n";
sub read_file {
    my ($path) = @_;
    open my $fh, '<', $path or die "cannot open $path: $!\n";
    local $/;
    my $text = <$fh>;
    close $fh or die "cannot close $path: $!\n";
    defined $text && length $text or die "$path is empty\n";
    return $text;
}

my (%sources, %definitions, %registered);
my @texts;
for my $source (split /\n/, read_file("$directory/sources.tsv")) {
    $source =~ /^[a-z][a-z0-9_]*\.c$/
        or die "invalid property source: $source\n";
    !$sources{$source}++ or die "duplicate property source: $source\n";
    my $text = read_file("$directory/$source");
    while ($text =~ /^TEST\s+(\w+)\s*\(void\)\s*\{/mg) {
        !$definitions{$1}++ or die "duplicate test definition: $1\n";
    }
    push @texts, $text;
}
opendir my $dir, $directory or die "cannot open $directory: $!\n";
for my $source (sort readdir $dir) {
    next unless $source =~ /\.c$/;
    $sources{$source} or die "uncompiled property source: $source\n";
}
closedir $dir or die "cannot close $directory: $!\n";
for my $entry (split /\n/, read_file("$directory/duckvep_property_tests.def")) {
    $entry =~ /^DUCKVEP_PROPERTY_TEST\(([A-Za-z_][A-Za-z0-9_]*)\)$/
        or die "invalid property test registration: $entry\n";
    my $name = $1;
    !$registered{$name}++ or die "duplicate test registration: $name\n";
    $definitions{$name} or die "registered test lacks definition: $name\n";
}
for my $name (sort keys %definitions) {
    $registered{$name} or die "unregistered test definition: $name\n";
}
print join("\n", @texts);
