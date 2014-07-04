#include "generic.h"

int bwscale_usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Apply scale for a bigWig file.\n\n");
    fprintf(stderr, "Usage:   methylQA bwscale [options] <bigWig> <scale>\n\n");
    fprintf(stderr, "Options: -o       output prefix [basename of input without extension_scaled]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main bwscale function */
int main_bwscale(int argc, char *argv[]){
    char *output, *outbwfile, *outbedgraph, *optoutput=NULL;
    int c;
    time_t start_time, end_time;
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "o:h?")) >= 0) {
        switch (c) {
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return bwscale_usage(); break;
            default: return 1;
        }
    }
    if (optind + 2 > argc)
        return bwscale_usage();
    char *bwfile = argv[optind];
    float scale = atof(argv[optind+1]);
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(bwfile)));
    }
    
    if(asprintf(&outbwfile, "%s_scaled.bw", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbedgraph, "%s.bedGraph.tmp", output) < 0)
        errAbort("Mem Error.\n");

    fprintf(stderr, "* Applying scale\n");
    bwScale(bwfile, outbwfile, outbedgraph, scale);

    unlink(outbedgraph);

    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

