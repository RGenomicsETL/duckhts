#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

sub write_file {
    my ($path, $text) = @_;
    open my $fh, '>', $path or die "cannot write $path: $!\n";
    print {$fh} $text;
    close $fh or die "cannot close $path: $!\n";
}
sub check_inventory {
    my ($files) = @_;
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/$_", $files->{$_}) for keys %$files;
    open my $saved_stderr, '>&', \*STDERR or die "cannot save stderr: $!\n";
    open STDERR, '>', "$dir/error" or die "cannot capture stderr: $!\n";
    open my $pipe, '-|', $^X, "$Bin/read_sources.pl", $dir
        or die "cannot start source audit: $!\n";
    my $output = do { local $/; <$pipe> };
    my $ok = close $pipe;
    open STDERR, '>&', $saved_stderr or die "cannot restore stderr: $!\n";
    open my $error, '<', "$dir/error" or die "cannot read stderr: $!\n";
    my $message = do { local $/; <$error> };
    close $error or die "cannot close stderr: $!\n";
    return ($ok, $output // '', $message // '');
}

my %files = (
    'sources.tsv' => "alpha.c\nbeta.c\n",
    'alpha.c' => "TEST alpha(void) { PASS(); }\n",
    'beta.c' => "TEST beta(void) { PASS(); }\n",
    'duckvep_property_tests.def' =>
        "DUCKVEP_PROPERTY_TEST(beta)\nDUCKVEP_PROPERTY_TEST(alpha)\n",
);
my ($ok, $output, $error) = check_inventory(\%files);
ok($ok, 'complete inventory passes independently of test registration order');
is($output, "$files{'alpha.c'}\n$files{'beta.c'}", 'all compiled source bytes retained');
is($error, '', 'valid inventory has no diagnostic');
my @mutations = (
    ['sources.tsv', "alpha.c\n", qr/uncompiled property source: beta.c/],
    ['sources.tsv', "alpha.c\nbeta.c\nalpha.c\n", qr/duplicate property source/],
    ['sources.tsv', "alpha.c\nbeta.c\nabsent.c\n", qr/cannot open .*absent.c/],
    ['sources.tsv', "../alpha.c\n", qr/invalid property source/],
    ['alpha.c', $files{'beta.c'}, qr/duplicate test definition: beta/],
    ['duckvep_property_tests.def', "DUCKVEP_PROPERTY_TEST(alpha)\n",
        qr/unregistered test definition: beta/],
    ['duckvep_property_tests.def', $files{'duckvep_property_tests.def'} .
        "DUCKVEP_PROPERTY_TEST(alpha)\n", qr/duplicate test registration: alpha/],
    ['duckvep_property_tests.def', "DUCKVEP_PROPERTY_TEST(absent)\n",
        qr/registered test lacks definition: absent/],
    ['duckvep_property_tests.def', "RUN_TEST(alpha)\n", qr/invalid property test registration/],
);
for my $mutation (@mutations) {
    my ($path, $replacement, $diagnostic) = @$mutation;
    ($ok, $output, $error) = check_inventory({%files, $path => $replacement});
    ok(!$ok, "reject $diagnostic");
    is($output, '', 'invalid inventory publishes no source evidence');
    like($error, $diagnostic, 'named inventory failure');
}
done_testing();
