#include "generic.h"

int medip_usage(){
    fprintf(stderr, "\nAnalyzing MeDIP-seq data, generating density and reports.\n");
    fprintf(stderr, "Please notice that if reads mapped to the chromosomes which are not in the size file, those reads will be discarded.\n\n");
    fprintf(stderr, "Usage:   methylQA medip [options] <chromosome size file> <bam/sam alignment file>\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -Q       unique reads mapping Quality threshold [10]\n");
    fprintf(stderr, "         -r       keep redundant reads [off]\n");
    fprintf(stderr, "         -T       treat 1 paired-end read as 2 single-end reads [off]\n");
    fprintf(stderr, "         -m       specify a CpG bed file for calculating CpG stats [null]\n");
    fprintf(stderr, "         -D       keep reads if only one end mapped in a pair [off]\n");
    fprintf(stderr, "         -C       add 'chr' string as prefix of reference sequence [off]\n");
    fprintf(stderr, "         -E       extend reads to represent fragment [150], specify 0 if want no extension\n");
    fprintf(stderr, "         -I       insert length threshold [500]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main function */
int main_medip (int argc, char *argv[]) {
    
    char *output, *outReportfile, *outExtfile, *outbedGraphfile, *outbigWigfile, *outCountfile, *outCovfile, *outInsertfile, *outGenomeCov;
    unsigned long long int *cnt;
    int optSam = 0, c, optDup = 1, optaddChr = 0, optDis = 1, optTreat = 0;
    unsigned int optQual = 10, optExt = 150, optisize = 500;
    char *optoutput = NULL, *optm = NULL;
    time_t start_time, end_time;
    start_time = time(NULL);
    struct slInt *cpgCount = NULL;
    struct slInt *slPair = NULL;
    long long fragbase = 0;
    int *covCnt=malloc(1024);
    long long *countCnt = malloc(1024);
    
    while ((c = getopt(argc, argv, "SQ:rTm:DCo:E:I:h?")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'Q': optQual = (unsigned int)strtol(optarg, 0, 0); break;
            case 'r': optDup = 0; break;
            case 'T': optTreat = 1; break;
            case 'm': optm = strdup(optarg); break;
            case 'D': optDis = 0; break;
            case 'C': optaddChr = 1; break;
            case 'E': optExt = (unsigned int)strtol(optarg, 0, 0); break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return medip_usage(); break;
            default: return 1;
        }
    }
    if (optind + 2 > argc)
        return medip_usage();

    char *chr_size_file = argv[optind];
    char *sam_file = argv[optind+1];
   
    //struct hash *genome = newHash(0);
    struct hash *hash = hashNameIntFile(chr_size_file);
    struct hash *cpgHash = newHash(0);
    if (optm != NULL){
        fprintf(stderr, "* CpG bed file %s provided, will calculate CpG stats\n", optm);
        fprintf(stderr, "* Reading the CpG bed file\n");
        cpgHash = cpgBed2BinKeeperHash(hash, optm);
    }

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
    if (asprintf(&outCountfile, "%s.cpgCount", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outCovfile, "%s.cpgCoverage", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outInsertfile, "%s.insertdistro", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outGenomeCov, "%s.genomeCoverage", output) < 0)
        errAbort("Preparing output wrong");
    
    
    //initilize genome coverage hash
    //genome = initGenomeCovHash(hash);
    //exit(1);

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

    if (optm != NULL){
        fprintf(stderr, "* Generating CpG stats\n");
        //writecpgCount(cpgCount, outCountfile);
        //writecpgCov(cpgHash, outCovfile);
        countCnt = plotcpgCount(cpgCount, output);
        covCnt = plotcpgCov(cpgHash, output);
        hashFree(&cpgHash);
        slFreeList(&cpgCount);
    }
    
    
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
    //writeGenomeCov(covhash, outGenomeCov);
    plotGenomeCov(covhash, output);
    
    //write report file
    fprintf(stderr, "* Preparing report file\n");
    writeReportDensity(outReportfile, cnt, optQual);

    // pdf report
    genMeDIPTex(output, optQual, cnt, fragbase, covCnt, countCnt, slPair, hash, covhash, optm);
    tex2pdf(output);

    
    
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
