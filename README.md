# methylQA (Version: 0.1.4 (r025))
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
