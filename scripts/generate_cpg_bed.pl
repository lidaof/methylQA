#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $inputfile;
my $debug; my $verbose; my $simulate;
my $cmd; my $ret;
use File::Basename;
my $self = bless {};
# ($name,$path,$suffix) = fileparse($inputfile,@suffixlist);
# $name = fileparse($fullname,@suffixlist);
# $basename = basename($fullname,@suffixlist);
# $dirname  = dirname($fullname);
my $chunk; my $add_chr; my $chop_chr;
my $diff = 3; my $plus_sign;
my $samtools = 'samtools';
my $seqtk    = 'seqtk';
GetOptions(

           'i|input|ref|inputfile|reference:s'  => \$inputfile,  # input file is the genome reference fasta file
           'chunk|region:s'                     => \$chunk,      # region of the genome to calculate
           'samtools:s'                         => \$samtools,   # pointer to samtools binary
           'seqtk:s'                            => \$seqtk,      # pointer to seqtk binary
           'add_chr'                            => \$add_chr,    # add chr prefix
           'chop_chr'                           => \$chop_chr,   # chop chr prefix
           'plus_sign'                          => \$plus_sign,  # add plus sign to bed 4th column
           'diff:s'                             => \$diff,
           'debug'                              => \$debug,
           'verbose'                            => \$verbose,
           'simulate'                           => \$simulate,

          );

die "No inputfile [$inputfile] -- $!" unless (-s $inputfile);
die "No chunk [$chunk] -- $!" unless (defined $chunk);

my ($chunk_chrom,$chunk_start_end) = split(':',$chunk);
if ($add_chr) { my $tmp = $chunk_chrom; $chunk_chrom = 'chr' . $tmp; }
my ($chunk_start,$chunk_end) = split('-',$chunk_start_end);
my $region = "$chunk_chrom:$chunk_start-$chunk_end";
$cmd = "$samtools faidx $inputfile $region | seqtk seq -";
$ret = `$cmd`; chomp $ret;
my ($header,$seq) = split("\n",$ret);
$header =~ s/\>//; $header =~ s/\:/\_/;
print STDERR "$seq\n" if ($debug);
my @vals = $self->all_match_positions('CG',$seq);

my $offset = $chunk_start - $diff;
my $pair;
my $plus_sign_opt = ''; $plus_sign_opt = '+' if ($plus_sign);
while ($pair = shift @vals) {
  next if (!defined $pair);
  my ($this_start,$this_end) = @$pair;
  my $ofstart = $offset+$this_start;
  my $ofend   = $offset+$this_end;
  print STDERR "#$ofstart,$this_start,$ofend,$this_end\n" if ($verbose);
  my $sign_val = $ofstart+1;
  my $sign    = $plus_sign_opt; $sign .= $sign_val;
  my $score   = 0;
  my $strand  = '+';
  $chunk_chrom =~ s/^chr// if ($chop_chr);
  print $chunk_chrom . "\t" . $ofstart . "\t" . $ofend . "\t" . $sign  . "\t" . $score  . "\t" . $strand . "\n";
}

1;

########################################
# METHODS

sub all_match_positions {
    my ($self, $regex, $string) = @_;
    my @ret;
    while ($string =~ /($regex)/g) {
        push @ret, [pos($string), pos($string) + length $1];
    }
    return @ret;
}

# generate_cpg_bed.pl
#
# Cared for by Albert Vilella <avilella@gmail.com>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

generate_cpg_bed.pl - DESCRIPTION 

=head1 SYNOPSIS

Requires:

samtools
seqtk
bedtools

perl methylqa/scripts/chr_chunk.pl -i /data/Genomes/Human/hs38DH/hg38.chr.genome | while read chunk; do perl methylqa/scripts/generate_cpg_bed.pl -plus_sign -add_chr -i /data/Genomes/Human/hs38DH/hs38DH.fa -chunk $chunk ; done | bedtools sort -i stdin | bgzip -c > hs38DH.cpg.chr.bed.gz

=head1 DESCRIPTION

GetOptions(

           'i|input|ref|inputfile|reference:s'  => \$inputfile,  # input file is the genome reference fasta file
           'chunk|region:s'                     => \$chunk,      # region of the genome to calculate
           'samtools:s'                         => \$samtools,   # pointer to samtools binary
           'seqtk:s'                            => \$seqtk,      # pointer to seqtk binary
           'add_chr'                            => \$add_chr,    # add chr prefix
           'chop_chr'                           => \$chop_chr,   # chop chr prefix
           'plus_sign'                          => \$plus_sign,  # add plus sign to bed 4th column
           'diff:s'                             => \$diff,
           'debug'                              => \$debug,
           'verbose'                            => \$verbose,
           'simulate'                           => \$simulate,


=head1 AUTHOR - Albert Vilella

Email avilella@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut

