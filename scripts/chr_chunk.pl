#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $inputfile = "$ENV{GROUP}/CEGXP/hg38.chr.genome";
my $debug; my $verbose; my $simulate;
my $cmd; my $ret;
use File::Basename;
# ($name,$path,$suffix) = fileparse($inputfile,@suffixlist);
# $name = fileparse($fullname,@suffixlist);
# $basename = basename($fullname,@suffixlist);
# $dirname  = dirname($fullname);
my $size = 1000000;
my $keep_chr; my $overlap = 0; my $add_chr;
GetOptions(

           'i|input|inputfile:s'  => \$inputfile, # genome sizes file
           'add_chr'              => \$add_chr,   # add chr prefix
           'keep_chr'             => \$keep_chr,  # keep the chr prefix
           'size:s'               => \$size,      # size of each chunk emitted
           'overlap:s'            => \$overlap,   # overlap between one chunk and the next, makes end pos longer
           'debug'                => \$debug,
           'verbose'              => \$verbose,
           'simulate'             => \$simulate,

          );

open IN, "$inputfile" or die $!;
while (<IN>) {
  my $line = $_; chomp $line;
  $line =~ s/chr// unless ($keep_chr);
  my ($chr,$length) = split("\t",$line);
  if ($add_chr) { my $tmp = $chr; $chr = 'chr' . $tmp; }
  my $start = 1;
  my $end   = $start + $size - 1 + $overlap;
  my $outformat = "$chr:$start-$end";
  do {
    print "$outformat\n";
    $start += $size;
    $end   += $size;
    $outformat = "$chr:$start-$end";
  } while ($end < $length);
}

# chr_chunk.pl
#
# Cared for by Albert Vilella <avilella@gmail.com>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

chr_chunk.pl - DESCRIPTION 

=head1 SYNOPSIS

perl chr_chunk.pl -i /my/hg38.chr.genome

=head1 DESCRIPTION

GetOptions(

           'i|input|inputfile:s'  => \$inputfile, # genome sizes file
           'add_chr'              => \$add_chr,   # add chr prefix
           'keep_chr'             => \$keep_chr,  # keep the chr prefix
           'size:s'               => \$size,      # size of each chunk emitted
           'overlap:s'            => \$overlap,   # overlap between one chunk and the next, makes end pos longer
           'debug'                => \$debug,
           'verbose'              => \$verbose,
           'simulate'             => \$simulate,

=head1 AUTHOR - Albert Vilella

Email avilella@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut

