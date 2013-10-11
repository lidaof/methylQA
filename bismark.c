#include "generic.h"

int bismark_usage(){
    fprintf(stderr, "\nAnalyzing SAM/BAM generated by Bismark.\nInput multilple bam files need be separated by comma, support at most 100 bam/sam files.\n");
    fprintf(stderr, "Default output stats over each CpG, using -b to generate bismark like output, which based on each C (in CpG)\n");
    fprintf(stderr, "This tool was way faster than bismark_methylation_extractor, but have limit functions.\nCurrently only works with bismark output with bowtie 1.\nUsers are recommended to use the original bismark_methylation_extractor tool for a comprehensive result.\n\n");
    fprintf(stderr, "Usage:   methylQA bismark [options] <chromosome size file> <CpG bed file> <bam/sam alignment file1,file2,file3...>\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -C       Add 'chr' string as prefix of reference sequence [off]\n");
    fprintf(stderr, "         -s       only output report stats [off]\n");
    fprintf(stderr, "         -b       bismark-like output [off]\n");
    fprintf(stderr, "         -I       Insert length threshold [500]\n");
    fprintf(stderr, "         -o       output prefix [basename of first input filename without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main function */
int main_bismark (int argc, char *argv[]) {
    
    char *output, *outReportfile, *outCpGfile, *outbedGraphfile, *row[100], *samfilecopy;
    unsigned long long int *cnt;
    unsigned long long int *cnt2;
    int optSam = 0, c, optaddChr = 0, optStats = 0, optBis = 0;
    unsigned int optisize = 500;
    char *optoutput = NULL;
    struct hash *cpgHash = newHash(0);
    time_t start_time, end_time;
    start_time = time(NULL);
    
    while ((c = getopt(argc, argv, "SCsbo:I:h?")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'C': optaddChr = 1; break;
            case 's': optStats = 1; break;
            case 'b': optBis = 1; break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return bismark_usage(); break;
            default: return 1;
        }
    }
    if (optind + 3 > argc)
        return bismark_usage();

    char *chr_size_file = argv[optind];
    char *cpg_bed_file = argv[optind+1];
    char *sam_file = argv[optind+2];

    fprintf(stderr, "* CpG file provided: %s\n", cpg_bed_file);
    fprintf(stderr, "* Insert size cutoff: %u\n", optisize);
   
    struct hash *chrHash = hashNameIntFile(chr_size_file);
    
    samfilecopy = cloneString(sam_file);
    int numFields = chopByChar(samfilecopy, ',', row, ArraySize(row));
    fprintf(stderr, "* Provided %i BAM/SAM file(s)\n", numFields);
    
    if(optStats) {
        fprintf(stderr, "* Warning: will report stats only as -s specified\n");
    }
    // if use select bismark like output, read cpgHash at each C stats
    if(optBis) {
        fprintf(stderr, "* Warning: will output stats over each C (in CpG)\n");
        cpgHash = cpgBed2BinKeeperHashBismark(chrHash, cpg_bed_file);
    }else{
        fprintf(stderr, "* Warning: will output stats over each CpG\n");
        cpgHash = cpgBed2BinKeeperHash(chrHash, cpg_bed_file);
    }

    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(row[0])));
    }
    

    if(asprintf(&outCpGfile, "%s.CpG.bedGraph", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbedGraphfile, "%s.density.bedGraph", output) < 0)
        errAbort("Mem Error.\n");
    if (asprintf(&outReportfile, "%s.report", output) < 0)
        errAbort("Preparing output wrong");
    

    //sam file to bed file
    //fprintf(stderr, "* Parsing the SAM/BAM file\n");
    cnt = bismarkBamParse(sam_file, cpgHash, optSam, optaddChr);
    
    //write to file
    cnt2 = writecpgBismark(cpgHash, outbedGraphfile, outCpGfile, optStats);
    
    //sort output
    if(!optStats) {
        fprintf(stderr, "* Sorting output density\n");
        sortBedfile(outbedGraphfile);
    }
    
    //sort output
    if(!optStats) {
        fprintf(stderr, "* Sorting output CpG methylation call\n");
        sortBedfile(outCpGfile);
    }

    //generate bigWig
    //fprintf(stderr, "* Generating bigWig\n");
    //bigWigFileCreate(outbedGraphfile, chr_size_file, 256, 1024, 0, 1, outbigWigfile);
    //bedGraphToBigWig(outbedGraphfile, chr_size_file, outbigWigfile);
    
    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReportBismark(outReportfile, cnt, cnt2, numFields, row, optBis, hashIntSum(chrHash));
    
    //cleaning
    hashFree(&chrHash);
    hashFree(&cpgHash);
    free(outCpGfile);
    free(outbedGraphfile);
    //free(outbigWigfile);
    free(outReportfile);
    free(samfilecopy);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}
