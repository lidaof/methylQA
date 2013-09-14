#include "generic.h"

static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: methylQA (methylation sequence data quality assessment tool from Wang lab)\n");
    fprintf(stderr, "Version: %s\n\n", methylQA_VERSION);
    fprintf(stderr, "Usage:   methylQA <command> [options]\n\n");
    fprintf(stderr, "Command: medip       generate genome density for MeDIP-Seq data\n");
    fprintf(stderr, "         mre         generate CpG density for MRE-Seq data\n");
    fprintf(stderr, "         density     generate genome density for ChIP-Seq data\n");
    fprintf(stderr, "         genomecov   calculate genome coverage by reads using sorted bedGraph file\n");
    fprintf(stderr, "         bismark     analysis bismark output\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc < 2) return usage();
    if (strcmp(argv[1], "medip") == 0) return main_medip(argc-1, argv+1);
    else if (strcmp(argv[1], "mre") == 0) return main_cpg(argc-1, argv+1);
    else if (strcmp(argv[1], "density") == 0) return main_density(argc-1, argv+1);
    else if (strcmp(argv[1], "genomecov") == 0) return main_genomecov(argc-1, argv+1);
    else if (strcmp(argv[1], "bismark") == 0) return main_bismark(argc-1, argv+1);
    else {
        fprintf(stderr, "[methylQA] unrecognized command '%s'\n", argv[1]);
        return 1;
    }
    return 0;
}
