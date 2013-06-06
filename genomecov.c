#include "generic.h"

int genomecov_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Calculate genome coverage using sorted bedGraph file.\n\n");
    fprintf(stderr, "Usage:   methylQA genomecov [options] <chromosome size file> <sorted bedGraph file>\n\n");
    fprintf(stderr, "Options: -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main genomecov function */
int main_genomecov(int argc, char *argv[]){
    char *output, *outCovfile, *optoutput=NULL;
    int c;
    time_t start_time, end_time;
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "o:h?")) >= 0) {
        switch (c) {
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return genomecov_usage(); break;
            default: return 1;
        }
    }
    if (optind + 2 > argc)
        return genomecov_usage();
    char *chrsizefile = argv[optind];
    char *bedgraph = argv[optind+1];
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(bedgraph)));
    }
    
    if(asprintf(&outCovfile, "%s.genomeCoverage", output) < 0)
        errAbort("Mem Error.\n");

    fprintf(stderr, "* Calculating genome coverage\n");
    struct hash *hash = calGenomeCovBedGraph(chrsizefile, bedgraph);
    writeGenomeCov(hash, outCovfile);

    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

