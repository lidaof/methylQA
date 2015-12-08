#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $inputfile;
my $debug; my $verbose; my $simulate;
my $cmd; my $ret;
my $self = bless {};
use File::Basename;
# ($name,$path,$suffix) = fileparse($inputfile,@suffixlist);
# $name = fileparse($fullname,@suffixlist);
# $basename = basename($fullname,@suffixlist);
# $dirname  = dirname($fullname);
my $outdir = "/bi/group/cegx/CEGXP";
my $outname;
GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
           'debug' => \$debug,
           'verbose' => \$verbose,
           'simulate' => \$simulate,
           'outname:s' => \$outname,
          );

$cmd = "samtools view -H $inputfile |";
open IN, "$cmd" or die $!;
while (<IN>) {
  my $line = $_;
  $self->{genome}{$1} = $2 if ($line =~ /^\@SQ\tSN\:(\S+)\tLN\:(\d+)/);
  $self->{reference} = $1 if ($line =~ /^\@PG.+reference\s+(\S+)\s+/);
}

my @suffixlist = ('.fa','.fasta','.fa.gz','.fasta.gz','.fas','.fas.gz');
my ($refname,$refpath,$refsuffix) = fileparse($self->{reference},@suffixlist);
my $outfile = $outdir . "/" . $refname . ".genome";
$outfile = $outdir . "/" . $outname if (defined $outname);

open OUT, ">$outfile" or die $!;
foreach my $contig (sort keys %{$self->{genome}}) {
  print OUT "$contig". "\t". $self->{genome}{$contig} . "\n";
}

close OUT;

print "outfile:\n";
print "$outfile\n";

# create_genome_sizes_from_bam.pl
#
# Cared for by Albert Vilella <avilella@gmail.com>
#
# Copyright Albert Vilella
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

create_genome_sizes_from_bam.pl - DESCRIPTION 

=head1 SYNOPSIS

perl create_genome_sizes_from_bam.pl -i /my/file.bam

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Albert Vilella

Email avilella@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=cut


