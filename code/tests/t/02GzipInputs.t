use strict;
use warnings;

use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

use_ok('GenomeProperties');

my $tempdir = tempdir(CLEANUP => 1);
my $fasta = ">seq1\nAAAA\n";
my $matches = join("\t", 'seq1', qw(md5 4 Pfam PF00001))."\n";

my $seq_file = "$tempdir/sequences.input";
my $match_file = "$tempdir/matches.input";
my $plain_seq_file = "$tempdir/sequences.fa";
my $plain_match_file = "$tempdir/matches.tsv";

open(my $seq_fh, '>', $plain_seq_file)
  or die "Could not create $plain_seq_file: $!";
print {$seq_fh} $fasta;
close($seq_fh);

open(my $match_fh, '>', $plain_match_file)
  or die "Could not create $plain_match_file: $!";
print {$match_fh} $matches;
close($match_fh);

ok(gzip(\$fasta => $seq_file), "Compress FASTA input")
  or diag($GzipError);
ok(gzip(\$matches => $match_file), "Compress match input")
  or diag($GzipError);

my $gp = GenomeProperties->new;
$gp->{seqs} = $seq_file;
$gp->{matches} = $match_file;

$gp->define_sequence_set;
ok(exists $gp->get_sequence_set->{seq1}, "Read sequences from gzip input");
ok($gp->signature_matches, "Read matches from gzip input");
is_deeply(
  $gp->get_sequence_set->{seq1},
  ['PF00001'],
  "Attach annotations read from gzip input",
);

my $plain_gp = GenomeProperties->new;
$plain_gp->{seqs} = $plain_seq_file;
$plain_gp->{matches} = $plain_match_file;

$plain_gp->define_sequence_set;
ok(
  exists $plain_gp->get_sequence_set->{seq1},
  "Read sequences from uncompressed input",
);
ok($plain_gp->signature_matches, "Read matches from uncompressed input");
is_deeply(
  $plain_gp->get_sequence_set->{seq1},
  ['PF00001'],
  "Attach annotations read from uncompressed input",
);

done_testing();
