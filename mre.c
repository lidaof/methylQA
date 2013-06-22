#include "generic.h"

int cpg_usage(){
    //fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   methylQA mre [options] <chromosome size file> <MRE fragment file> <bam/sam alignment file>\n\n");
    fprintf(stderr, "Generating MRE CpG density.\n\n");
    fprintf(stderr, "Options: -S       input is SAM [off]\n");
    fprintf(stderr, "         -Q       unique reads mapping Quality threshold [10]\n");
    fprintf(stderr, "         -c       base calling from which base [1]\n");
    fprintf(stderr, "         -n       MRE fragment minimal length cutoff [50]\n");
    fprintf(stderr, "         -x       MRE fragment maximal length cutoff [500]\n");
    fprintf(stderr, "         -R       remove redundant reads [off]\n");
    fprintf(stderr, "         -T       treat 1 paired-end read as 2 single-end reads [off]\n");
    fprintf(stderr, "         -m       specify a CpG bed file for calculating CpG stats [null]\n");
    fprintf(stderr, "         -D       discard if only one end mapped in a paired end reads [off]\n");
    fprintf(stderr, "         -C       Add 'chr' string as prefix of reference sequence [off]\n");
    //fprintf(stderr, "         -E       output each MRE enzyme reads [off]\n");
    fprintf(stderr, "         -I       Insert length threshold [500]\n");
    fprintf(stderr, "         -o       output prefix [basename of input without extension]\n");
    fprintf(stderr, "         -h       help message\n");
    fprintf(stderr, "         -?       help message\n");
    fprintf(stderr, "\n");
    return 1;
}

/* main stat function */
int main_cpg (int argc, char *argv[]) {
    
    char *output, *outbigWig, *outbedGraph, *outBed, *outFilterBed;
    unsigned long long int *cnt2, *cnt, cnt1;
    int optSam = 0, c, optDup = 0, optaddChr = 0, optDis = 0, optTreat = 0, optMin = 50, optMax = 500, optCall = 1; //optEach = 0;
    unsigned int optQual = 10, optisize = 500;
    char *optoutput = NULL, *optm = NULL;
    time_t start_time, end_time;
    struct hash *hash = newHash(0);
    long long *cnt3;
    start_time = time(NULL);
    while ((c = getopt(argc, argv, "SQ:c:n:x:RTm:DCo:I:h?")) >= 0) {
        switch (c) {
            case 'S': optSam = 1; break;
            case 'Q': optQual = (unsigned int)strtol(optarg, 0, 0); break;
            case 'c': optCall = (int)strtol(optarg, 0, 0); break;
            case 'n': optMin = (int)strtol(optarg, 0, 0); break;
            case 'x': optMax = (int)strtol(optarg, 0, 0); break;
            case 'R': optDup = 1; break;
            case 'T': optTreat = 1; break;
            case 'm': optm = strdup(optarg); break;
            case 'D': optDis = 1; break;
            case 'C': optaddChr = 1; break;
            //case 'E': optEach = 1; break;
            case 'I': optisize = (unsigned int)strtol(optarg, 0, 0); break;
            case 'o': optoutput = strdup(optarg); break;
            case 'h':
            case '?': return cpg_usage(); break;
            default: return 1;
        }
    }
    if (optind + 3 > argc)
        return cpg_usage();

    char *chr_size_file = argv[optind];
    char *mre_frag_file = argv[optind+1];
    char *sam_file = argv[optind+2];
    
    if(optoutput) {
        output = optoutput;
    } else {
        output = cloneString(get_filename_without_ext(basename(sam_file)));
    }
    
    if(asprintf(&outbigWig, "%s.CpG.bigWig", output) < 0)
        errAbort("Mem Error.\n");
    if (asprintf(&outFilterBed, "%s.filter.bed", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outBed, "%s.bed", output) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&outbedGraph, "%s.CpG.bedGraph", output) < 0)
        errAbort("Preparing output wrong");

    struct hash *chrHash = hashNameIntFile(chr_size_file);
    
    fprintf(stderr, "* Start to parse the MRE fragments file\n");
    hash = MREfrag2Hash(mre_frag_file, optMin, optMax);
    
    //sam file
    fprintf(stderr, "* Start to parse the SAM/BAM file\n");
    cnt2 = sam2bed(sam_file, outBed, chrHash, optSam, optQual, optDup, optaddChr, optDis, optisize, 0, optTreat);

    fprintf(stderr, "* Filtering reads by MRE site\n");
    //if (optEach){
        cnt = filterReadByMREsite(hash, outBed, outFilterBed, optCall - 1, output);
    //} else {
    //    cnt = filterReadByMREsite(hash, outBed, outFilterBed, optCall - 1, "NULL");
    //}

    fprintf(stderr, "* Generating CpG bedGraph\n");
    cnt1 = CpGscorebedGraph(hash, cnt, outbedGraph);

    //write report file
    fprintf(stderr, "* fragment stats and preparing report file\n");
    struct fragd *fragdistro = fragmentStats(hash, cnt2, optQual, cnt, cnt1, output, 40, 400, 20);
    
    fprintf(stderr, "* Generating bigWig files\n");
    sortBedfile(outbedGraph);
    //bigWigFileCreate(outbedGraph, chr_size_file, 256, 1024, 0, 1, outbigWig);
    bedGraphToBigWig(outbedGraph, chr_size_file, outbigWig);

    struct hash *cpgHash = newHash(0);
    if (optm != NULL){
        fprintf(stderr, "* CpG bed file %s provided, will calculate CpG stats\n", optm);
        fprintf(stderr, "* Reading the CpG bed file\n");
        cpgHash = cpgBed2BinKeeperHash(chrHash, optm);
    }
    fprintf(stderr, "* Calculating CpG stats over filtered bed\n");
    cnt3 = bedCpGstat(cpgHash, outFilterBed);
    genMRETex(output, cnt2, cnt, cnt1, chrHash, cpgHash, cnt3, fragdistro);
    tex2pdf(output);
    
    //cleaning
    hashFree(&chrHash);
    hashFree(&hash);
    free(outbigWig);
    free(outBed);
    free(outbedGraph);
    free(outFilterBed);
    end_time = time(NULL);
    fprintf(stderr, "* Done, time used %.0f seconds.\n", difftime(end_time, start_time));
    return 0;
}

