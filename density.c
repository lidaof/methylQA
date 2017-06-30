#include "generic.h"

int density_usage(){
    fprintf(stderr, "\nAnalyzing ChIP-seq data, generating density and reports.\n");
    fprintf(stderr, "Please notice that if reads mapped to the chromosomes which are not in the size file, those reads will be discarded.\n\n");
    fprintf(stderr, "Usage:   methylQA density [options] <chromosome size file> <bam/sam alignment file>\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -Q       unique reads mapping Quality threshold [10]\n");
    fprintf(stderr, "         -r       do NOT remove redundant reads [off]\n");
    fprintf(stderr, "         -T       treat 1 paired-end read as 2 single-end reads [off]\n");
    fprintf(stderr, "         -D       do NOT discard if only one end mapped in a paired end reads [off]\n");
    fprintf(stderr, "         -C       Add 'chr' string as prefix of reference sequence [off]\n");
    fprintf(stderr, "         -E       extend reads to represent fragment [150], specify 0 if want no extension\n");
    fprintf(stderr, "         -I       Insert length threshold [500]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main function */
int main_density (int argc, char *argv[]) {
    
    char *output, *outReportfile, *outExtfile, *outbedGraphfile, *outbigWigfile, *outInsertfile, *outGenomeCov;
    unsigned long long int *cnt;
    int optSam = 0, c, optDup = 1, optaddChr = 0, optDis = 1, optTreat = 0;
    unsigned int optQual = 10, optExt = 150, optisize = 500;
    char *optoutput = NULL;
    time_t start_time, end_time;
    start_time = time(NULL);
    struct slInt *cpgCount = NULL;
    struct slInt *slPair = NULL;
    long long fragbase = 0;
    
    while ((c = getopt(argc, argv, "SQ:rTDCo:E:I:h?")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'Q': optQual = (unsigned int)strtol(optarg, 0, 0); break;
            case 'r': optDup = 0; break;
            case 'T': optTreat = 1; break;
            case 'D': optDis = 0; break;
            case 'C': optaddChr = 1; break;
            case 'E': optExt = (unsigned int)strtol(optarg, 0, 0); break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return density_usage(); break;
            default: return 1;
        }
    }
    if (optind + 2 > argc)
        return density_usage();

    char *chr_size_file = argv[optind];
    char *sam_file = argv[optind+1];
   
    //struct hash *genome = newHash(0);
    struct hash *hash = hashNameIntFile(chr_size_file);
    struct hash *cpgHash = newHash(0);

    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(sam_file)));
    }
    

    if(asprintf(&outExtfile, "%s.extended.bed", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbedGraphfile, "%s.extended.bedGraph", output) < 0)
        errAbort("Mem Error.\n");
    if(asprintf(&outbigWigfile, "%s.bigWig", output) < 0)
        errAbort("Mem Error.\n");
    if (asprintf(&outReportfile, "%s.report", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outInsertfile, "%s.insertdistro", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outGenomeCov, "%s.genomeCoverage", output) < 0)
        errAbort("Preparing output wrong");
    

    //sam file to bed file
    fprintf(stderr, "* Parsing the SAM/BAM file\n");
    cnt = sam2bedwithCpGstat(sam_file, outExtfile, hash, cpgHash, &cpgCount, &slPair, optSam, optQual, optDup, optaddChr, optDis, optisize, optExt, optTreat);
    //sort
    //fprintf(stderr, "\n* Sorting\n");
    //bedSortFile(outBedfile, outBedfile);

    //remove dup
    //fprintf(stderr, "* Removing duplication\n");
    //uniqueBed = removeBedDup(outBedfile, outFilterfile);

    //extend and write extend bed
    //fprintf(stderr, "* Extending to %d and writing extended bed\n", arguments.extend);
    //int extendWarn = extendBed(hash, arguments.extend, outFilterfile, outExtfile);
    //if (extendWarn == 1)
    //    outExtfile = cloneString(outFilterfile);
    //
    //if (extendWarn != 1){
        //sort extend bed
    //    fprintf(stderr, "* Sorting extended bed\n");
    //    bedSortFile(outExtfile, outExtfile);
    //}

    

    plotMappingStat(cnt, output);

    if(slPair != NULL){
        fprintf(stderr, "* Generating fragments size stats\n");
        //writeInsertsize(slPair, outInsertfile);
        fragbase = plotInsertsize(slPair, output); //quite time consuming -- fixed
    }
    fprintf(stderr, "* fragments total base: %lli\n", fragbase);
    
    //sort extend bed
    fprintf(stderr, "* Sorting extended bed\n");
    sortBedfile(outExtfile);
    
    //bedItemOverlap step
    fprintf(stderr, "* Generating bedGraph\n");
    bedItemOverlapCount(hash, outExtfile, outbedGraphfile);

    //generate bigWig
    fprintf(stderr, "* Generating bigWig\n");
    //bigWigFileCreate(outbedGraphfile, chr_size_file, 256, 1024, 0, 1, outbigWigfile);
    bedGraphToBigWig(outbedGraphfile, chr_size_file, outbigWigfile);

    fprintf(stderr, "* Calculating genome coverage\n");
    struct hash *covhash = calGenomeCovBedGraph(chr_size_file, outbedGraphfile);
    writeGenomeCov(covhash, outGenomeCov);
    //plotGenomeCov(covhash, output);
    
    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReportDensity(outReportfile, cnt, optQual);
    
    //cleaning
    hashFree(&hash);
    free(outExtfile);
    free(outbedGraphfile);
    free(outbigWigfile);
    free(outReportfile);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}
