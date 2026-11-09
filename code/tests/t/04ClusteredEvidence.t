use Test::More tests => 17;
use Cwd 'abs_path';
use File::Temp qw/ tempdir /;

use strict;
use warnings;

#Work out where the data is
my $data_dir = abs_path($0);
$data_dir    =~ s|04ClusteredEvidence\.t|data|;

my $tempdir = tempdir( CLEANUP =>  1 );

use_ok( 'GenomeProperties');

#Step 3 of the test property only has clustered evidence (formerly HMM-CLUST).
#gntR_near is 600 bp from phnW_1; gntR_far is on a different contig.
sub run_gp {
  my (%extra) = @_;
  my $acc = delete($extra{acc}) // 'GenProp0713';
  my $gp = GenomeProperties->new;
  $gp->set_options({ matches  => $data_dir."/04_matches",
                     all      => 1,
                     name     => "clustered",
                     outfiles => ["summary"],
                     outdir   => $tempdir,
                     gpdir    => $data_dir,
                     gpff     => "04_gpff",
                     %extra }, "cal");
  $gp->open_outputfiles;
  $gp->read_properties;
  $gp->evaluate_properties;
  $gp->close_outputfiles;
  my $def = $gp->get_def($acc);
  my ($step) = grep { $_->order eq '3' } @{ $def->get_steps };
  return ($def, $step);
}

my ($def, $step) = run_gp();
ok($step->found, "Without gene positions the default treats clustered evidence as plain evidence");
is($def->result, 'YES', "Default without gene positions keeps the previous result");

($def, $step) = run_gp( clustered => 'strict' );
ok(!$step->found, "Strict mode without gene positions ignores clustered evidence");
is($def->result, 'PARTIAL', "Required clustered step missing gives PARTIAL");

($def, $step) = run_gp( gff => $data_dir."/04_gff" );
ok($step->found, "With gene positions the default checks proximity and the nearby gene counts");
is($def->result, 'YES', "Property complete with gene positions");
is_deeply($step->clustered_matches('PF00392'), ['gntR_near'],
  "Only the gene near the rest of the property is reported");

($def, $step) = run_gp( gff => $data_dir."/04_gff", cluster_distance => 500 );
ok(!$step->found, "Clustered evidence is rejected beyond cluster_distance");
is($def->result, 'PARTIAL', "Property partial when the gene is too far away");

($def, $step) = run_gp( clustered => 'permissive' );
ok($step->found, "Permissive mode treats clustered evidence as plain evidence");
is($def->result, 'YES', "Permissive mode reproduces the behaviour before clustered evidence");

#GenProp0736 step 3 has two clustered models, but only PF00392 is in the genome.
#Strict mode follows the JCVI HMM-CLUST rule: any one nearby clustered match is enough.
#Permissive mode keeps the plain-evidence rule, where every non-sufficient evidence must match.
($def, $step) = run_gp( acc => 'GenProp0736', gff => $data_dir."/04_gff" );
ok($step->found, "Strict mode: one nearby clustered match fills the step");

($def, $step) = run_gp( acc => 'GenProp0736', clustered => 'permissive' );
ok(!$step->found, "Permissive mode: all non-sufficient evidence must match, as before");

my $gp = GenomeProperties->new;
eval { $gp->set_options({ matches => $data_dir."/04_matches", all => 1, clustered => 'loose' }, "cal") };
like($@, qr/strict or permissive/, "Unknown clustered mode is rejected");

my $bad = GenomeProperties->new;
$bad->set_options({ matches => $data_dir."/04_matches", all => 1,
                    gpdir   => $data_dir, gpff => "04_gpff_bad_qualifier" }, "cal");
eval { $bad->read_properties };
like($@, qr/Unknown evidence qualifier \|nearby\|/, "Unknown EV qualifier is rejected");

my $ev = (grep { $_->order eq '3' } @{ $def->get_steps })[0]->get_evidence->[0];
ok($ev->clustered && !$ev->sufficient, "Evidence object reports the clustered qualifier");
