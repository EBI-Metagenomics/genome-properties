use Test::More tests => 9;
use Cwd 'abs_path';
use File::Basename;
use File::Slurp;
use File::Temp qw/ tempdir /;

use strict;
use warnings;

#Regression test for https://github.com/ebi-pf-team/genome-properties/issues/30
#
#Perl's evaluate_step used to OR together every piece of evidence in a step
#regardless of whether it was marked "sufficient", meaning a matched
#non-sufficient evidence could incorrectly prove a step even when the
#sufficient evidence for that step was absent. This test locks in the fixed
#behaviour, which matches the pygenprop reimplementation:
#  - If a step has one or more evidences marked "sufficient", only those
#    evidences can prove the step (any match among them is enough).
#  - If a step has no "sufficient" evidences at all, every evidence on the
#    step must be matched for the step to succeed.

my $data_dir = abs_path($0);
$data_dir    =~ s|03SufficientEvidence\.t|data|;

my $tempdir = tempdir( CLEANUP => 0 );

use_ok('GenomeProperties');

my $gp = GenomeProperties->new;
isa_ok($gp, "GenomeProperties");

ok($gp->set_options({ matches  => $data_dir."/03_matches",
                     all      => 1,
                     name     => "test",
                     outfiles => ["summary"],
                     outdir   => $tempdir,
                     gpdir    => $data_dir,
                     gpff     => "03_gpff" },                 "cal"), "Set options");

ok($gp->open_outputfiles, "Open output files");
ok($gp->read_properties, "read_properties");
ok($gp->evaluate_properties, "Evaluate properties");
ok($gp->write_results, "Writing out results");
ok($gp->close_outputfiles, "Close file handle");

my $expected_results = read_file($data_dir."/03_summary_result");
my $observed_results = read_file($tempdir."/SUMMARY_FILE_test");

is($observed_results, $expected_results,
   "Sufficient/non-sufficient evidence steps resolve correctly (issue #30)");
