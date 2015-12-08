## methylation sequence data quality assessment tool

### install
simply run:

    make

will generate the binary ready for use. Please visit official site http://methylqa.sourceforge.net/ for documentation.

### update
After get the updated source (either by `git pull` or downloading the tarbal), run

    make clean 
    make

will generate the new version of binary.
(or `make cleanlocal` if no changes need for kent and samtools)

### Used Libraries
Libraries from UCSC Kent source (http://genome.ucsc.edu/admin/jk-install.html) and Samtools (http://samtools.sourceforge.net/) were included in this package. Files were distributed to their own licenses. methyQA currently use Kent source verson v284 and samtools version 0.1.18, the included library here were customized to reduce file size.

### Contact
Questions or comments, please contact our mailing list (methylQA@googlegroups.com), thanks!

### Addendum: generate cpg files for a reference genome

Given a bam file produced by aligning reads to a bisulfite converted reference genome, one can produce the required bed file describing the coordinates of CpG positions and the genome sizes files needed for methylqa bismark mode.

To produce the genome sizes file, use:

create_genome_sizes_from_bam.pl

To produce the cpg bed file, use:

chr_chunk.pl and generate_cpg_bed.pl

This will produce the files needed to run methylqa bismark mode together with the input bam file.
