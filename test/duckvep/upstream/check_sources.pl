#!/usr/bin/env perl

use strict;
use warnings;

use Digest::SHA qw(sha256_hex);
use File::Basename qw(dirname);
use File::Find qw(find);
use File::Spec;
use FindBin qw($Bin);

my $verify_upstream_git = 0;
if (@ARGV) {
  @ARGV == 1 && $ARGV[0] eq '--verify-upstream-git'
    or die "usage: $0 [--verify-upstream-git]\n";
  $verify_upstream_git = 1;
}

my $manifest = File::Spec->catfile($Bin, 'sources.tsv');
open my $fh, '<', $manifest or die "cannot open $manifest: $!\n";

my $header = <$fh>;
defined $header or die "$manifest is empty\n";
chomp $header;
my @expected_header = qw(repository_url ref commit upstream_path local_path sha256 role);
my @header = split /\t/, $header, -1;
die "$manifest has an unexpected header\n"
  unless join("\t", @header) eq join("\t", @expected_header);

my %authority = (
  'https://github.com/Ensembl/ensembl-vep' => {
    ref => 'release/116',
    commit => '57ea5c52340acc1f156267f810ad162e26597082',
    local_prefix => 'ensembl-vep/release-116/',
    expected_count => 49,
    expected_role => 'oracle_self_test',
    git_env => 'DUCKVEP_ENSEMBL_VEP_GIT',
    complete_test_prefix => 't/',
  },
  'https://github.com/Ensembl/ensembl-variation' => {
    ref => 'release/116',
    commit => '2fb834b987ede3824e200197a838ce11e91aeb4b',
    local_prefix => 'ensembl-variation/release-116/',
    expected_count => 2,
    expected_role => 'semantic_fixture_source',
    expected_paths => {
      'modules/t/hgvs_parser.t' => 1,
      'modules/t/variation_effect.t' => 1,
    },
    git_env => 'DUCKVEP_ENSEMBL_VARIATION_GIT',
  },
  'https://github.com/Ensembl/ensembl-tools' => {
    ref => 'release/89',
    commit => '9d50991e514ea414981d3b5513256330daadce78',
    local_prefix => 'ensembl-tools/release-89/',
    expected_count => 5,
    expected_role => 'historical_lineage',
    git_env => 'DUCKVEP_ENSEMBL_TOOLS_GIT',
    complete_test_prefix => 'scripts/variant_effect_predictor/t/',
  },
);
my %roles = map { $_ => 1 } qw(oracle_self_test semantic_fixture_source historical_lineage);
my %receipted;
my %authority_count;
my %upstream_path;
my %receipt_sha256;
my $line_number = 1;

while (my $line = <$fh>) {
  $line_number++;
  chomp $line;
  next if $line eq '';
  my @field = split /\t/, $line, -1;
  die "$manifest:$line_number expected seven tab-separated columns\n" unless @field == 7;
  my ($repository_url, $ref, $commit, $upstream_path, $local_path, $sha256, $role) = @field;

  my $pin = $authority{$repository_url}
    or die "$manifest:$line_number unknown authority $repository_url\n";
  die "$manifest:$line_number ref $ref does not match $pin->{ref}\n"
    unless $ref eq $pin->{ref};
  die "$manifest:$line_number commit $commit does not match $pin->{commit}\n"
    unless $commit eq $pin->{commit};
  die "$manifest:$line_number unknown role $role\n" unless $roles{$role};
  die "$manifest:$line_number role $role does not match $pin->{expected_role}\n"
    unless $role eq $pin->{expected_role};
  die "$manifest:$line_number local path must end in .t\n" unless $local_path =~ /\.t\z/;
  die "$manifest:$line_number upstream path must end in .t\n" unless $upstream_path =~ /\.t\z/;
  die "$manifest:$line_number local path does not preserve the authority/ref prefix\n"
    unless index($local_path, $pin->{local_prefix}) == 0;
  my $expected_local = $pin->{local_prefix} . $upstream_path;
  die "$manifest:$line_number local path $local_path does not preserve $upstream_path\n"
    unless $local_path eq $expected_local;
  die "$manifest:$line_number invalid SHA-256 $sha256\n" unless $sha256 =~ /\A[0-9a-f]{64}\z/;
  die "$manifest:$line_number duplicate local path $local_path\n" if $receipted{$local_path};

  my $path = File::Spec->catfile($Bin, split m{/}, $local_path);
  open my $source, '<:raw', $path or die "cannot open mirrored source $path: $!\n";
  local $/;
  my $bytes = <$source>;
  close $source or die "cannot close mirrored source $path: $!\n";
  my $actual = sha256_hex($bytes);
  die "$path SHA-256 $actual does not match receipt $sha256\n" unless $actual eq $sha256;
  $receipted{$local_path} = 1;
  $authority_count{$repository_url}++;
  $upstream_path{$repository_url}->{$upstream_path} = 1;
  $receipt_sha256{$repository_url}->{$upstream_path} = $sha256;
}
close $fh or die "cannot close $manifest: $!\n";

for my $repository_url (sort keys %authority) {
  my $actual = $authority_count{$repository_url} || 0;
  my $expected = $authority{$repository_url}->{expected_count};
  die "$manifest receipts $actual files for $repository_url; expected $expected\n"
    unless $actual == $expected;
  if (defined $authority{$repository_url}->{expected_paths}) {
    my @missing = grep {
      !$upstream_path{$repository_url}->{$_}
    } sort keys %{ $authority{$repository_url}->{expected_paths} };
    my @extra = grep {
      !$authority{$repository_url}->{expected_paths}->{$_}
    } sort keys %{ $upstream_path{$repository_url} };
    die "$manifest selects the wrong semantic sources for $repository_url: " .
        "missing=[" . join(', ', @missing) . "] extra=[" .
        join(', ', @extra) . "]\n"
      if @missing || @extra;
  }
}

if ($verify_upstream_git) {
  for my $repository_url (sort keys %authority) {
    my $pin = $authority{$repository_url};
    my $repo = $ENV{$pin->{git_env}};
    defined $repo && $repo ne ''
      or die "$pin->{git_env} is required by --verify-upstream-git\n";
    -d $repo or die "$pin->{git_env} is not a directory: $repo\n";

    open my $commit, '-|', 'git', '-C', $repo, 'cat-file', '-e',
      "$pin->{commit}^{commit}"
      or die "cannot start git cat-file in $repo: $!\n";
    close $commit
      or die "$repo does not contain commit $pin->{commit}\n";

    for my $path (sort keys %{ $upstream_path{$repository_url} }) {
      open my $git, '-|', 'git', '-C', $repo, 'show',
        "$pin->{commit}:$path"
        or die "cannot start git show for $repository_url $path: $!\n";
      binmode $git;
      local $/;
      my $bytes = <$git>;
      close $git
        or die "cannot read $path from $repository_url commit $pin->{commit}\n";
      my $actual = sha256_hex($bytes);
      my $expected = $receipt_sha256{$repository_url}->{$path};
      die "$repository_url $pin->{commit}:$path SHA-256 $actual " .
          "does not match receipt $expected\n"
        unless $actual eq $expected;
    }

    next unless defined $pin->{complete_test_prefix};
    open my $tree, '-|', 'git', '-C', $repo, 'ls-tree', '-r', '--name-only',
      $pin->{commit}, '--', $pin->{complete_test_prefix}
      or die "cannot start git ls-tree in $repo: $!\n";
    my %expected_path;
    while (my $path = <$tree>) {
      chomp $path;
      $expected_path{$path} = 1 if $path =~ /\.t\z/;
    }
    close $tree
      or die "cannot enumerate $repository_url commit $pin->{commit}\n";
    my @missing = grep {
      !$upstream_path{$repository_url}->{$_}
    } sort keys %expected_path;
    my @extra = grep {
      !$expected_path{$_}
    } sort keys %{ $upstream_path{$repository_url} };
    die "$repository_url mirror inventory differs from $pin->{commit}: " .
        "missing=[" . join(', ', @missing) . "] extra=[" .
        join(', ', @extra) . "]\n"
      if @missing || @extra;
  }
}

my @unreceipted;
for my $entry (keys %authority) {
  my $prefix = $authority{$entry}->{local_prefix};
  my $root = File::Spec->catdir($Bin, split m{/}, $prefix);
  find(
    sub {
      return unless -f $_ && $_ =~ /\.t\z/;
      my $relative = File::Spec->abs2rel($File::Find::name, $Bin);
      $relative =~ s{\\}{/}g;
      push @unreceipted, $relative unless $receipted{$relative};
    },
    $root,
  );
}

die 'unreceipted upstream tests: ' . join(', ', sort @unreceipted) . "\n" if @unreceipted;
die "$manifest contains no upstream test receipts\n" unless keys %receipted;

my $self_test_receipts = File::Spec->catfile($Bin, 'self_test_receipts.tsv');
open my $self_test, '<', $self_test_receipts
  or die "cannot open $self_test_receipts: $!\n";
my $self_test_header = <$self_test>;
defined $self_test_header or die "$self_test_receipts is empty\n";
chomp $self_test_header;
my @expected_self_test_header = qw(
  run_date repository_url ref commit test_selector module_distribution
  module_distribution_sha256 perl_version environment_lock_path files assertions
  skipped_assertions skip_reasons result tap_path tap_sha256
);
die "$self_test_receipts has an unexpected header\n"
  unless $self_test_header eq join("\t", @expected_self_test_header);

my $self_test_count = 0;
my $vep_repository = 'https://github.com/Ensembl/ensembl-vep';
while (my $line = <$self_test>) {
  chomp $line;
  next if $line eq '';
  my @field = split /\t/, $line, -1;
  die "$self_test_receipts expected 16 tab-separated columns\n" unless @field == 16;
  my ($run_date, $repository_url, $ref, $commit, $test_selector,
      $module_distribution, $module_sha256, $perl_version,
      $environment_lock_path, $files, $assertions, $skipped, $skip_reasons,
      $result, $tap_path, $tap_sha256) = @field;
  my $pin = $authority{$repository_url}
    or die "$self_test_receipts names unknown authority $repository_url\n";
  die "$self_test_receipts self-test must name $vep_repository\n"
    unless $repository_url eq $vep_repository;
  die "$self_test_receipts receipt must use the pinned ref and commit\n"
    unless $ref eq $pin->{ref} && $commit eq $pin->{commit};
  die "$self_test_receipts has invalid run date $run_date\n"
    unless $run_date =~ /\A\d{4}-\d{2}-\d{2}\z/;
  die "$self_test_receipts must identify the complete VEP test selector\n"
    unless $test_selector eq 't/*.t';
  die "$self_test_receipts lacks a versioned module distribution\n"
    unless $module_distribution =~ /ensembl-vep=116\.0/;
  die "$self_test_receipts has invalid module SHA-256\n"
    unless $module_sha256 =~ /\A[0-9a-f]{64}\z/;
  die "$self_test_receipts has invalid Perl version\n"
    unless $perl_version =~ /\A\d+\.\d+\.\d+\z/;
  die "$self_test_receipts has invalid environment lock path\n"
    unless $environment_lock_path =~
      /\Areceipts\/[A-Za-z0-9_.-]+[.]conda-explicit[.]txt\z/;
  my $environment_lock = File::Spec->catfile(
    $Bin,
    split m{/},
    $environment_lock_path,
  );
  open my $lock, '<:raw', $environment_lock
    or die "cannot open environment lock $environment_lock: $!\n";
  my @lock_lines = <$lock>;
  close $lock or die "cannot close environment lock $environment_lock: $!\n";
  chomp @lock_lines;
  my @explicit = grep { $_ eq '@EXPLICIT' } @lock_lines;
  die $environment_lock . ' must contain exactly one @EXPLICIT marker' . "\n"
    unless @explicit == 1;
  my @package_url = grep { /\Ahttps:\/\// } @lock_lines;
  die "$environment_lock contains a non-URL package entry\n"
    if grep { $_ ne '' && $_ !~ /\A#/ && $_ ne '@EXPLICIT' &&
      $_ !~ /\Ahttps:\/\// } @lock_lines;
  my ($module_name, $module_version, $module_build) =
    $module_distribution =~ /\A[^:]+::([^=]+)=([^=]+)=([^=]+)\z/;
  defined $module_build
    or die "$self_test_receipts has malformed module distribution\n";
  my $module_package = join('-', $module_name, $module_version, $module_build);
  my @module_url = grep { m{/\Q$module_package\E[.](?:conda|tar[.]bz2)\z} }
    @package_url;
  die "$environment_lock must pin exactly one $module_package package\n"
    unless @module_url == 1;
  my @perl_url = grep {
    m{/perl-\Q$perl_version\E-[^/]+[.](?:conda|tar[.]bz2)\z}
  } @package_url;
  die "$environment_lock must pin exactly one Perl $perl_version package\n"
    unless @perl_url == 1;
  die "$self_test_receipts has invalid test counters\n"
    unless $files =~ /\A\d+\z/ && $assertions =~ /\A\d+\z/ &&
      $skipped =~ /\A\d+\z/ && $files == 49 && $assertions >= $skipped;
  my $skip_sum = 0;
  my %receipt_skip;
  for my $reason (split /;/, $skip_reasons, -1) {
    $reason =~ /\A(.+):(\d+)\z/
      or die "$self_test_receipts has malformed skip reason $reason\n";
    my ($label, $count) = ($1, $2);
    die "$self_test_receipts duplicates skip reason $label\n"
      if exists $receipt_skip{$label};
    $receipt_skip{$label} = $count;
    $skip_sum += $count;
  }
  die "$self_test_receipts skip reasons sum to $skip_sum, not $skipped\n"
    unless $skip_sum == $skipped;
  die "$self_test_receipts records non-passing result $result\n" unless $result eq 'PASS';
  die "$self_test_receipts has invalid TAP path $tap_path\n"
    unless $tap_path =~ /\Areceipts\/[A-Za-z0-9_.-]+\.tap\z/;
  die "$self_test_receipts has invalid TAP SHA-256\n"
    unless $tap_sha256 =~ /\A[0-9a-f]{64}\z/;
  my $tap_file = File::Spec->catfile($Bin, split m{/}, $tap_path);
  open my $tap, '<:raw', $tap_file or die "cannot open TAP receipt $tap_file: $!\n";
  local $/;
  my $tap_bytes = <$tap>;
  close $tap or die "cannot close TAP receipt $tap_file: $!\n";
  my $actual_tap_sha256 = sha256_hex($tap_bytes);
  die "$tap_file SHA-256 $actual_tap_sha256 does not match $tap_sha256\n"
    unless $actual_tap_sha256 eq $tap_sha256;
  my ($tap_files, $tap_assertions, $tap_result);
  my $tap_assertion_lines = 0;
  my $tap_plan_files = 0;
  my $tap_plan_assertions = 0;
  my $tap_file_passes = 0;
  my $tap_success_summaries = 0;
  my %tap_file;
  my %tap_skip;
  my $active_tap_file;
  my $active_assertions = 0;
  my $active_next_number = 1;
  my $active_plan;
  for my $tap_line (split /\n/, $tap_bytes, -1) {
    die "$tap_file contains a failed TAP status\n"
      if $tap_line =~ /\A\s*not ok(?:\s|\z)/;
    die "$tap_file contains a TAP bailout\n"
      if $tap_line =~ /\A\s*Bail out!/i;
    if ($tap_line =~ /\A(t\/[A-Za-z0-9_.\/-]+\.t)\s+\.*\s*\z/) {
      defined $active_tap_file
        and die "$tap_file starts $1 before $active_tap_file terminated\n";
      exists $upstream_path{$vep_repository}->{$1}
        or die "$tap_file executes unpinned test $1\n";
      die "$tap_file duplicates executed test $1\n" if $tap_file{$1}++;
      $active_tap_file = $1;
      $active_assertions = 0;
      $active_next_number = 1;
      $active_plan = undef;
      next;
    }
    if (defined $active_tap_file &&
        $tap_line =~ /\A(ok|not ok)\s+(\d+)(.*)\z/) {
      my ($status, $number, $rest) = ($1, $2, $3);
      $number == $active_next_number
        or die "$tap_file $active_tap_file expected assertion " .
          "$active_next_number, observed $number\n";
      $active_next_number++;
      $active_assertions++;
      $tap_assertion_lines++;
      die "$tap_file $active_tap_file contains failed assertion $number\n"
        if $status eq 'not ok';
      if ($rest =~ /#\s*skip\s*(.*?)\s*\z/i) {
        $tap_skip{$1}++;
      }
      next;
    }
    if (defined $active_tap_file &&
        $tap_line =~ /\A1\.\.(\d+)(?:\s|\z)/) {
      defined $active_plan
        and die "$tap_file $active_tap_file duplicates its TAP plan\n";
      $active_plan = $1;
      $tap_plan_files++;
      $tap_plan_assertions += $1;
      next;
    }
    if (defined $active_tap_file && $tap_line eq 'ok') {
      defined $active_plan && $active_plan == $active_assertions
        or die "$tap_file $active_tap_file plan does not match its assertions\n";
      $tap_file_passes++;
      $active_tap_file = undef;
      next;
    }
    if ($tap_line =~ /\AFiles=(\d+), Tests=(\d+),/) {
      defined $tap_files
        and die "$tap_file contains duplicate prove summaries\n";
      ($tap_files, $tap_assertions) = ($1, $2);
    }
    if ($tap_line =~ /\AResult:\s+(\S+)\s*\z/) {
      defined $tap_result
        and die "$tap_file contains duplicate result summaries\n";
      $tap_result = $1;
    }
    $tap_success_summaries++ if $tap_line eq 'All tests successful.';
  }
  defined $active_tap_file
    and die "$tap_file ends before $active_tap_file terminated\n";
  defined $tap_files && defined $tap_assertions && defined $tap_result
    or die "$tap_file lacks a complete prove summary\n";
  die "$tap_file summary does not match its receipt\n"
    unless $tap_files == $files && $tap_assertions == $assertions &&
      $tap_assertion_lines == $assertions && $tap_result eq $result;
  my @missing_tap_file = grep {
    !$tap_file{$_}
  } sort keys %{ $upstream_path{$vep_repository} };
  my @extra_tap_file = grep {
    !$upstream_path{$vep_repository}->{$_}
  } sort keys %tap_file;
  die "$tap_file executed-file inventory does not match the pinned VEP suite: " .
      "missing=[" . join(', ', @missing_tap_file) . "] extra=[" .
      join(', ', @extra_tap_file) . "]\n"
    if @missing_tap_file || @extra_tap_file;
  die "$tap_file contains failed assertions or incomplete test-file statuses\n"
    unless $tap_plan_files == $files &&
      $tap_plan_assertions == $assertions &&
      $tap_file_passes == $files &&
      $tap_success_summaries == 1;
  die "$tap_file skip reasons do not match its receipt\n"
    unless join("\t", map { $_ . '=' . $tap_skip{$_} } sort keys %tap_skip) eq
      join("\t", map { $_ . '=' . $receipt_skip{$_} } sort keys %receipt_skip);
  $self_test_count++;
}
close $self_test or die "cannot close $self_test_receipts: $!\n";
die "$self_test_receipts must contain exactly one pinned VEP receipt\n"
  unless $self_test_count == 1;

print 'verified ', scalar(keys %receipted),
  ' checksum-pinned upstream Ensembl test sources',
  ($verify_upstream_git ? ' against exact Git commits' : ''), "\n";
