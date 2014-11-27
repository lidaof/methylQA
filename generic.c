#include "generic.h"

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

static int binOffsetsExtended[] =
	{4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};

#define _binFirstShift 17	/* How much to shift to get to finest bin. */
#define _binNextShift 3		/* How much to shift to get to next larger bin. */

char template_name[]="/tmp/methylQAXXXXXX";
char *cpglabel[19] = {"0","1","2","3","4","5","6","7","8","9","10","11-20","21-30","31-40","41-50","51-100","101-200","201-300","\\textgreater 300"};

/* definitions of functions */

char *strrev(char *str){
    int i = 0, j = strlen(str);
    char *rev = malloc(sizeof(char)*(j + 1));
    //while(str[i] != '\0')
    while(j > 0)
        rev[i++] = str[--j];
    //printf("%i\n", i);
    rev[i] = '\0';
    return rev;
}

char *get_filename_without_ext(char *filename) {
    char *s;
    s = malloc(strlen(filename) + 1);
    strcpy(s, filename);
    char *dot = strrchr(s, '.');
    if(!dot || dot == s) return s;
    *dot = '\0';
    return s;
}

char *get_filename_ext(char *filename) {
    char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

char * texTitleEscape(char *title){
    return replaceChars(title, "_", "\\_");
}

bool is_file(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

struct lineFile *lineFileOpen2(char *fileName, bool zTerm){
/* Open up a lineFile or die trying. */
if (is_dir(fileName))
    errAbort("Error: %s is a directory not a file", fileName);
struct lineFile *lf = lineFileMayOpen(fileName, zTerm);
if (lf == NULL)
    errAbort("Couldn't open %s , %s", fileName, strerror(errno));
return lf;
}

boolean binKeeperAnyInclude(struct binKeeper *bk, int start, int end){
/* Return TRUE if start/end includes any items in binKeeper. */
struct binElement *el;
int startBin, endBin;
int i,j, len;

if (start < bk->minPos) start = bk->minPos;
if (end > bk->maxPos) end = bk->maxPos;
if (start >= end) return FALSE;
startBin = (start>>_binFirstShift);
endBin = ((end-1)>>_binFirstShift);
for (i=0; i<ArraySize(binOffsetsExtended); ++i)
    {
    int offset = binOffsetsExtended[i];
    for (j=startBin+offset; j<=endBin+offset; ++j)
        {
	for (el=bk->binLists[j]; el != NULL; el = el->next)
	    {
            len = el->end - el->start;
	    if (rangeIntersection(el->start, el->end, start, end) >= len)
	        {
		return TRUE;
		}
	    }
	}
    startBin >>= _binNextShift;
    endBin >>= _binNextShift;
    }
return FALSE;
}

int binKeeperCpGstat(struct binKeeper *bk, int start, int end) {
    /* get stat for cpg */
    struct binElement *el;
    int startBin, endBin;
    int i,j, len, c=0;
    //struct slInt *cc;
    if (start < bk->minPos) start = bk->minPos;
    if (end > bk->maxPos) end = bk->maxPos;
    startBin = (start>>_binFirstShift);
    endBin = ((end-1)>>_binFirstShift);
    for (i=0; i<ArraySize(binOffsetsExtended); ++i){
        int offset = binOffsetsExtended[i];
        for (j=startBin+offset; j<=endBin+offset; ++j){
            for (el=bk->binLists[j]; el != NULL; el = el->next){
                len = el->end - el->start;
                if (rangeIntersection(el->start, el->end, start, end) >= len){
                    struct cpgC *oc = (struct cpgC *) el->val;
                    (oc->c)++;
                    c++;
        	}
           }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    //cc = slIntNew(c);
    //slAddHead(&cpgCount, cc);
    return c;
    //fprintf(stderr, "read have %d CpG\n", c);
}

void writeReportDensity(char *outfile, unsigned long long int *cnt, unsigned int mapQ){
    FILE *f = mustOpen(outfile, "w");
    fprintf(f, "total reads (pair): %llu\n", cnt[0]);
    fprintf(f, "  read ends 1: %llu\n", cnt[0]);
    fprintf(f, "  read ends 2: %llu\n", cnt[1]);
    fprintf(f, "  mapped read ends 1: %llu\n", cnt[2]);
    fprintf(f, "  mapped read ends 2: %llu\n", cnt[3]);
    fprintf(f, "  used read ends 1: %llu\n", cnt[4]);
    fprintf(f, "  used read ends 2: %llu\n", cnt[5]);
    fprintf(f, "mappable reads (pair): %llu\n", cnt[6]);
    //fprintf(f, "non-redundant mappable reads (pair): %llu\n", cnt[8]);
    fprintf(f, "uniquely mapped reads (pair) (mapQ >= %u): %llu\n", mapQ, cnt[7]);
    fprintf(f, "non-redundant uniquely mapped reads (pair): %llu\n", cnt[9]);
    carefulClose(&f);
}

void writeReportBismark(char *outfile, unsigned long long int *cnt, unsigned long long int *cnt2, int numFields, char *row[100], int bisMode, long long genomeBase){
    FILE *f = mustOpen(outfile, "w");
    int i;
    fprintf(f, "files provided: %i\n", numFields);
    for(i = 0; i < numFields; i++) fprintf(f, "  %s\n", row[i]);
    fprintf(f, "\n");
    fprintf(f, "uniquely mappable reads (pair): %llu\n", cnt[0]);
    fprintf(f, "quality failed mapped reads (pair) in the bismark bam: %llu\n", cnt[1]);
    fprintf(f, "oversized mapped reads (pair) in the bismark bam: %llu\n", cnt[4]);
    //fprintf(f, "duplicated mapped reads (pair) (per lane based): %llu\n", cnt[2]);
    fprintf(f, "total base of uniquely mapped reads (pair): %llu\n", cnt[3]);
    fprintf(f, "total base of uniquely mapped reads (pair) cover genome base (%lli): %.1fX\n", genomeBase, cnt[3]*1.0/genomeBase);
    fprintf(f, "\n");
    fprintf(f, "in all uniquely mapped reads (pair), found:\n");
    fprintf(f, "    number of methylated C in CHG context (was protected): %llu\n", cnt[5]);            // X for 
    fprintf(f, "    number of not methylated C in CHG context (was converted): %llu\n", cnt[6]);        // x for 
    fprintf(f, "        C->T convertion rate in CHG context: %.2f%%\n", cnt[6]*100.0/(cnt[6]+cnt[5]));
    fprintf(f, "    number of methylated C in CHH context (was protected): %llu\n", cnt[7]);            // H for 
    fprintf(f, "    number of not methylated C in CHH context (was converted): %llu\n", cnt[8]);        // h for 
    fprintf(f, "        C->T convertion rate in CHH context: %.2f%%\n", cnt[8]*100.0/(cnt[8]+cnt[7]));
    fprintf(f, "    number of methylated C in CpG context (was protected): %llu\n", cnt[9]);            // Z for 
    fprintf(f, "    number of not methylated C in CpG context (was converted): %llu\n", cnt[10]);        // z for 
    fprintf(f, "        C->T convertion rate in CpG context: %.2f%%\n", cnt[10]*100.0/(cnt[10]+cnt[9]));
    fprintf(f, "    number of methylated C in Unknown context (was protected): %llu\n", cnt[11]);        // U for 
    fprintf(f, "    number of not methylated C in Unknown context (was converted): %llu\n", cnt[12]);    // u for 
    fprintf(f, "        C->T convertion rate in Unknown context: %.2f%%\n", cnt[12]*100.0/(cnt[12]+cnt[11]));
    fprintf(f, "\n");
    if (cnt2 != NULL){
        if (bisMode){
            fprintf(f, "in the total %llu CpG Cytosine:\n", cnt2[19]);
        } else {
            fprintf(f, "in the total %llu CpG:\n", cnt2[19]);
        }
        fprintf(f, "%15s\t%15s\t%10s\t%c\n", "Times covered", "Count", "Percent", '|');
        for(i = 0; i < 18; i++){
            fprintf(f, "%15s\t%15llu\t%10.2f\t|%s\n", cpglabel[i], cnt2[i], cnt2[i]*100.0/cnt2[19], print_bar((int)(cnt2[i]*100.0/cnt2[19])));
        }
        fprintf(f, "%15s\t%15llu\t%10.2f\t|%s\n", ">300", cnt2[18], cnt2[18]*100.0/cnt2[19], print_bar((int)(cnt2[18]*100.0/cnt2[19])));
    }
    carefulClose(&f);
}

long long writeInsertsize(struct slInt *slPair, char *outfile){
    struct slInt *c;
    long long sum = 0;
    FILE *f = mustOpen(outfile, "w");
    for ( c = slPair; c != NULL; c = c->next){
        fprintf(f, "%d\n", c->val);
        sum += (long long) c->val;
    }
    carefulClose(&f);
    return sum;
}

long long plotInsertsize(struct slInt *slPair, char *prefix){
    long long sum = 0;
    char tmpRfile[50], tmpifile[50];
    strcpy(tmpRfile, template_name);
    strcpy(tmpifile, template_name);
    int fd = mkstemp(tmpRfile);
    if (fd == -1)
        errAbort("create temp file error.");
    int fd2 = mkstemp(tmpifile);
    if (fd2 == -1)
        errAbort("create temp file error.");
    sum = writeInsertsize(slPair, tmpifile);
    FILE *fout = fdopen(fd, "w");
    fprintf(fout, "dat <- read.table(\"%s\")\n", tmpifile);
    fprintf(fout, "pdf('%s.insertdistro.pdf')\n", prefix);
    fprintf(fout, "d <- density(dat$V1)\n");
    fprintf(fout, "plot(d, main=\"Fragments Size Distribution\")\n");
    fprintf(fout, "polygon(d, col=\"red\", border=\"blue\")\n");
    fprintf(fout, "dev.off()\n");
    fclose(fout);

    char *command;
    if (asprintf(&command, "Rscript %s", tmpRfile) < 0)
        errAbort("Preparing command wrong");
    if (system(command) == -1)
        fprintf(stderr, "failed to call R for plotting");
    unlink(tmpifile);
    unlink(tmpRfile);
    return sum;
}

void writecpgCount(struct slInt *cpgCount, char *outfile){
    struct slInt *c;
    //fprintf(stderr, "%d elements in cpgCount", slCount(cpgCount));
    FILE *f = mustOpen(outfile, "w");
    for ( c = cpgCount; c != NULL; c = c->next){
        fprintf(f, "%d\n", c->val);
    }
    carefulClose(&f);
}

long long * plotcpgCount(struct slInt *Count, char *prefix){
    long long *cnt = malloc(sizeof(long long)*20);
    struct slInt *c;
    long long sum = 0;
    int i, j;
    for(i=0;i<20;i++)
        cnt[i] = 0;
    for ( c = Count; c != NULL; c = c->next){
        j = c->val;
        sum += j;
        if (j == 0) cnt[0] ++;
        else if (j == 1) cnt[1]++;
        else if (j == 2) cnt[2]++;
        else if (j == 3) cnt[3]++;
        else if (j == 4) cnt[4]++;
        else if (j == 5) cnt[5]++;
        else if (j == 6) cnt[6]++;
        else if (j == 7) cnt[7]++;
        else if (j == 8) cnt[8]++;
        else if (j == 9) cnt[9]++;
        else if (j == 10) cnt[10]++;
        else if (j >=11 && j <= 20) cnt[11]++;
        else if (j >=21 && j <= 30) cnt[12]++;
        else if (j >=31 && j <= 40) cnt[13]++;
        else if (j >=41 && j <= 50) cnt[14]++;
        else if (j >=51 && j <= 100) cnt[15]++;
        else if (j >=101 && j <= 200) cnt[16]++;
        else if (j >=201 && j <= 300) cnt[17]++;
        else cnt[18]++; // j >=301
    }
    cnt[19] = sum;

    char tmpRfile[50];
    strcpy(tmpRfile, template_name);
    int fd = mkstemp(tmpRfile);
    if (fd == -1)
        errAbort("create temp file error.");
    FILE *fout = fdopen(fd, "w");
    fprintf(fout, "dat<-data.frame(label=c('0','1','2','3','4','5','6','7','8','9','10','11-20','21-30','31-40','41-50','51-100','101-200','201-300','>300'), \nvalue=c(");
    for(i=0;i<19;i++){
        fprintf(fout,"%lli", cnt[i]);
        if (i<18) fprintf(fout, ",");
    }
    fprintf(fout, "))\n");
    fprintf(fout, "dat <- subset(dat, value!=0)\n");
    fprintf(fout, "pdf('%s.cpgCount.pdf')\n", prefix);
    fprintf(fout, "op <- par(mar = c(5,7,4,2) + 0.1)\n");
    fprintf(fout, "barplot(rev(dat$value), names=rev(dat$label), main=\"CpG Count\", xlab=\"Number of Fragments\", ylab=\"\", col=4, las=1, horiz=TRUE)\n");
    fprintf(fout, "title(ylab = \"CpG Count in Fragments\", cex.lab = 1.5, line = 4.5)\n");
    fprintf(fout, "par(op)\n");
    fprintf(fout, "dev.off()\n");
    fclose(fout);

    char *command;
    if (asprintf(&command, "Rscript %s", tmpRfile) < 0)
        errAbort("Preparing command wrong");
    if (system(command) == -1)
        fprintf(stderr, "failed to call R for plotting");
    unlink(tmpRfile);
    return cnt;
}

int * plotcpgCov(struct hash *cpgHash, char *prefix){
    struct hashEl *hel;
    int *cnt = malloc(sizeof(int)*19);
    int i, j;
    for(i=0;i<19;i++)
        cnt[i] = 0;
    struct hashCookie cookie = hashFirst(cpgHash);
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) hel->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            struct cpgC *oc = (struct cpgC *) be->val;
            j = oc->c;
            if (j == 0) cnt[0] ++;
            else if (j == 1) cnt[1]++;
            else if (j == 2) cnt[2]++;
            else if (j == 3) cnt[3]++;
            else if (j == 4) cnt[4]++;
            else if (j == 5) cnt[5]++;
            else if (j == 6) cnt[6]++;
            else if (j == 7) cnt[7]++;
            else if (j == 8) cnt[8]++;
            else if (j == 9) cnt[9]++;
            else if (j == 10) cnt[10]++;
            else if (j >=11 && j <= 20) cnt[11]++;
            else if (j >=21 && j <= 30) cnt[12]++;
            else if (j >=31 && j <= 40) cnt[13]++;
            else if (j >=41 && j <= 50) cnt[14]++;
            else if (j >=51 && j <= 100) cnt[15]++;
            else if (j >=101 && j <= 200) cnt[16]++;
            else if (j >=201 && j <= 300) cnt[17]++;
            else cnt[18]++; // j >=301
        }
        binKeeperFree(&bk);
    }
    char tmpRfile[50];
    strcpy(tmpRfile, template_name);
    int fd = mkstemp(tmpRfile);
    if (fd == -1)
        errAbort("create temp file error.");
    FILE *fout = fdopen(fd, "w");
    fprintf(fout, "dat<-data.frame(label=c('0','1','2','3','4','5','6','7','8','9','10','11-20','21-30','31-40','41-50','51-100','101-200','201-300','>300'), \nvalue=c(");
    for(i=0;i<19;i++){
        fprintf(fout,"%i", cnt[i]);
        if (i<18) fprintf(fout, ",");
    }
    fprintf(fout, "))\n");
    fprintf(fout, "dat <- subset(dat, value!=0)\n");
    fprintf(fout, "pdf('%s.cpgCoverage.pdf')\n", prefix);
    fprintf(fout, "op <- par(mar = c(5,7,4,2) + 0.1)\n");
    fprintf(fout, "barplot(rev(dat$value), names=rev(dat$label), main=\"CpG Coverage\", xlab=\"Number of CpG\", ylab=\"\", col=3, las=1, horiz=TRUE)\n");
    fprintf(fout, "title(ylab = \"Times Covered\", cex.lab = 1.5, line = 4.5)\n");
    fprintf(fout, "par(op)\n");
    fprintf(fout, "dev.off()\n");
    fclose(fout);

    char *command;
    if (asprintf(&command, "Rscript %s", tmpRfile) < 0 )
        errAbort("Preparing command wrong");
    if (system(command) == -1)
        fprintf(stderr, "failed to call R for plotting");
    unlink(tmpRfile);
    return cnt;
}

void writecpgCov(struct hash *cpgHash, char *outfile){
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(cpgHash);
    FILE *f = mustOpen(outfile, "w");
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) hel->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            struct cpgC *oc = (struct cpgC *) be->val;
            fprintf(f, "%i\n", oc->c);
        }
        binKeeperFree(&bk);
    }
    carefulClose(&f);
}

unsigned long long int *writecpgBismark(struct hash *cpgHash, char *outfile, char *outcpg, int statsOnly, int covThres){
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int)*20);
    int i, j = 0;
    unsigned long long int tot=0; //tot will be the 20th element, which was not the coverage count, was the total C (in CpG) count
    for(i=0;i<20;i++)
        cnt[i] = 0;
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(cpgHash);
    FILE *f = NULL, *f2 = NULL;
    if(!statsOnly) {
        f = mustOpen(outfile, "w");
        f2 = mustOpen(outcpg, "w");
    }
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) hel->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            tot++;
            struct cpgC *oc = (struct cpgC *) be->val;
            if (oc->mc > 0 || oc->umc > 0){
                j = oc->mc + oc->umc;
                if(!statsOnly) {
                    if (j >= covThres){
                        fprintf(f, "%s\t%i\t%i\t%i\n", hel->name, be->start, be->end, j);
                        fprintf(f2, "%s\t%i\t%i\t%.4f\n", hel->name, be->start, be->end, (float)(oc->mc)/j);
                    }
                }
                //if (j == 0) cnt[0] ++;
                if (j == 1) cnt[1]++;
                else if (j == 2) cnt[2]++;
                else if (j == 3) cnt[3]++;
                else if (j == 4) cnt[4]++;
                else if (j == 5) cnt[5]++;
                else if (j == 6) cnt[6]++;
                else if (j == 7) cnt[7]++;
                else if (j == 8) cnt[8]++;
                else if (j == 9) cnt[9]++;
                else if (j == 10) cnt[10]++;
                else if (j >=11 && j <= 20) cnt[11]++;
                else if (j >=21 && j <= 30) cnt[12]++;
                else if (j >=31 && j <= 40) cnt[13]++;
                else if (j >=41 && j <= 50) cnt[14]++;
                else if (j >=51 && j <= 100) cnt[15]++;
                else if (j >=101 && j <= 200) cnt[16]++;
                else if (j >=201 && j <= 300) cnt[17]++;
                else cnt[18]++; // j >=301
            } else {
                cnt[0]++; 
            }
        }
        binKeeperFree(&bk);
    }
    cnt[19] = tot;
    if(!statsOnly) {
        carefulClose(&f2);
        carefulClose(&f);
    }
    return cnt;
}

void writecpgBismarkLite(struct hash *cpgHash, char *outfilefor, char *outfilerev, int covThres){
    int j = 0;
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(cpgHash);
    FILE *f = mustOpen(outfilefor, "w");
    FILE *f2 = mustOpen(outfilerev, "w");
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) hel->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            struct cpgC *oc = (struct cpgC *) be->val;
            if (oc->mc > 0 || oc->umc > 0){
                j = oc->mc + oc->umc;
                if (j >= covThres){
                    if (oc->strand == '+'){
                        fprintf(f, "%s\t%i\t%i\t%.4f\n", hel->name, be->start, be->end, (float)(oc->mc)/j);
                        //fprintf(f, "%s\t%i\t%i\t%i\t%i\t%.4f\n", hel->name, be->start, be->end, oc->mc, oc->umc, (float)(oc->mc)/j);
                    }else{
                        fprintf(f2, "%s\t%i\t%i\t%.4f\n", hel->name, be->start, be->end, (float)(oc->mc)/j);
                        //fprintf(f2, "%s\t%i\t%i\t%i\t%i\t%.4f\n", hel->name, be->start, be->end, oc->mc, oc->umc, (float)(oc->mc)/j);
                    }
                }
            }
        }
        binKeeperFree(&bk);
    }
    carefulClose(&f);
    carefulClose(&f2);
}

void writecpgBismarkLiteHash(struct hash *cpgHash, char *outfilefor, char *outfilerev, int covThres){
    int j = 0, k;
    struct hashEl *hel, *hel2;
    struct hashCookie cookie = hashFirst(cpgHash);
    FILE *f = mustOpen(outfilefor, "w");
    FILE *f2 = mustOpen(outfilerev, "w");
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct hash *hash2 = (struct hash *) hel->val;
        struct hashCookie cookie2 = hashFirst(hash2);
        while( (hel2 = hashNext(&cookie2)) != NULL ){
            struct cpgC *oc = (struct cpgC *) hel2->val;
            //if (oc->mc > 0 || oc->umc > 0){
            if (oc->mc > 0){
                j = oc->mc + oc->umc;
                if (j >= covThres) {
                    k = (int)strtol(hel2->name, 0, 0);
                    if (oc->strand == '+'){
                        fprintf(f, "%s\t%i\t%i\t%.4f\n", hel->name, k, k+1, (float)(oc->mc)/j);
                        //fprintf(f, "%s\t%i\t%i\t%i\t%i\t%.4f\n", hel->name, be->start, be->end, oc->mc, oc->umc, (float)(oc->mc)/j);
                    }else{
                        fprintf(f2, "%s\t%i\t%i\t%.4f\n", hel->name, k, k+1, (float)(oc->mc)/j);
                        //fprintf(f2, "%s\t%i\t%i\t%i\t%i\t%.4f\n", hel->name, be->start, be->end, oc->mc, oc->umc, (float)(oc->mc)/j);
                    }
                }
            }
        }
        hashFree(&hash2);
    }
    carefulClose(&f);
    carefulClose(&f2);
}

void writeGenomeCov(struct hash *cov, char *outfile){
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(cov);
    FILE *f = mustOpen(outfile, "w");
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct gcov *g = (struct gcov *) hel->val;
        fprintf(f, "%s\t%i\t%i\n", hel->name, g->total, g->cov);
    }
    carefulClose(&f);
}

void plotGenomeCov(struct hash *cov, char *prefix){
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(cov);
    int i = 1;
    int tot = hashNumEntries(cov);
    char tmpRfile[50];
    strcpy(tmpRfile, template_name);
    int fd = mkstemp(tmpRfile);
    if (fd == -1)
        errAbort("create temp file error.");
    FILE *fout = fdopen(fd, "w");
    fprintf(fout, "dat <- data.frame(label=rep(NA, %d), total=rep(0, %d), cov=rep(0, %d), stringsAsFactors=FALSE)\n", tot, tot, tot);
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct gcov *g = (struct gcov *) hel->val;
        fprintf(fout, "dat[%d, ] <- c('%s', %i, %i)\n", i, hel->name, g->total, g->cov);
        i++;
    }
    fprintf(fout, "dat$total <- as.numeric(dat$total)\n");
    fprintf(fout, "dat$cov <- as.numeric(dat$cov)\n");
    fprintf(fout, "pdf('%s.genomeCov.pdf')\n", prefix);
    fprintf(fout, "op <- par(mar = c(5,7,4,2) + 0.1)\n");
    fprintf(fout, "barplot(rev(dat$cov/dat$total), names=rev(dat$label), main=\"Genome Coverage\", xlab=\"Coverage Percentage\", ylab=\"\", col=5, las=1, horiz=TRUE, xlim=c(0,1))\n");
    fprintf(fout, "title(ylab = \"Chromosomes\", cex.lab = 1.5, line = 4.5)\n");
    fprintf(fout, "par(op)\n");
    fprintf(fout, "dev.off()\n");
    fclose(fout);

    char *command;
    if (asprintf(&command, "Rscript %s", tmpRfile) < 0 )
        errAbort("Preparing command wrong");
    if (system(command) == -1)
        fprintf(stderr, "failed to call R for plotting");
    unlink(tmpRfile);
}

void plotMappingStat(unsigned long long int *cnt, char *prefix){
    char tmpRfile[50];
    strcpy(tmpRfile, template_name);
    int fd = mkstemp(tmpRfile);
    if (fd == -1)
        errAbort("create temp file error.");
    FILE *fout = fdopen(fd, "w");
    fprintf(fout, "dat <- data.frame(count=c(%llu, %llu, %llu, %llu), label=c('total\\nfragments','mapped\\nfragments','uniquely\\nmapped\\nfragments', 'non-redundant\\nuniquely\\nmapped\\nfragments'))\n", cnt[0], cnt[6], cnt[7], cnt[9]);
    fprintf(fout, "dat$count <- as.numeric(dat$count)\n");
    fprintf(fout, "pdf('%s.mappingStat.pdf')\n", prefix);
    fprintf(fout, "op <- par(mar = c(7,8,4,2) + 0.1)\n");
    fprintf(fout, "bar <- barplot(dat$count, las=1, col=2:5, ylab='',xlab='', main='Mapping stats')\n");
    fprintf(fout, "title(ylab = 'Fragments count', cex.lab = 1.5, line = 4.5)\n");
    fprintf(fout, "axis(1, at=bar, labels=dat$label, padj=1, tick=FALSE)\n");
    fprintf(fout, "par(op)\n");
    fprintf(fout, "dev.off()\n");
    fclose(fout);

    char *command;
    if (asprintf(&command, "Rscript %s", tmpRfile) < 0)
        errAbort("Preparing command wrong");
    if (system(command) == -1)
        fprintf(stderr, "failed to call R for plotting");
    unlink(tmpRfile);
}

void assignCpGcount(struct hash *chrHash, struct hash *cpgHash, struct hash *chgHash, struct hash *chhHash, char *chrom, int start, char *methycall, char strand, int left, int right, unsigned long long int *methyCnt, int fullMode){
    int i, j;
    char key[20];
    struct binElement *hitList = NULL, *hit; // *hitList2 = NULL, *hit2;
    struct hashEl *hel = hashLookup(cpgHash, chrom);
    if (hel == NULL){
        return;
    }
    for(i = 0; i < strlen(methycall); i++){
        if (i < left) continue;
        if (i >= right) continue;
        j = methycall[i];
        if(j == 'Z' || j == 'z'){
            struct binKeeper *bk = (struct binKeeper *) hel->val;
            hitList = binKeeperFind(bk, start+i, start+i+1); //bismark have methyl info on each site, not each cpg -- FIXing
            if (hitList != NULL){
                //fprintf(stdout, "found C: %s %i %i with status %c\n", chrom, start+i, start+i+1, j);
                for (hit = hitList; hit !=NULL; hit = hit->next) {
                    //fprintf(stdout, "hit C: %s %i %i from hash\n", chrom, hit->start, hit->end);
                    struct cpgC *cg = (struct cpgC *) hit->val;
                    if (j == 'Z'){
                        (cg->mc)++;
                        methyCnt[4]++;
                    }else if(j == 'z'){
                        (cg->umc)++;
                        methyCnt[5]++;
                    }
                    //fprintf(stdout, "mC: %i umC: %i\n", cg->mc, cg->umc);
                    //break; // should be ok to comment out since just 1 CpG
                }
            }else{
                warn("not a CpG Cytosine found: %s %i %i", chrom, start+i, start+i+1);
                continue;
            }
        } else {
            if(sprintf(key, "%i", start+i) < 0)
                errAbort("Mem Error.\n");
            if(j == 'X'){
                methyCnt[0]++;
                if (fullMode){
                    /*
                    struct hashEl *hel2 = hashLookup(chgHash, chrom);
                    if (hel2 != NULL) {
                        struct binKeeper *bk2 = (struct binKeeper *) hel2->val;
                        hitList2 = binKeeperFind(bk2, start+i, start+i+1);
                        if (hitList2 == NULL){
                            struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                            c->c = 0;
                            c->mc = 1;
                            c->umc = 0;
                            c->strand = strand;
                            binKeeperAdd(bk2, start+i, start+i+1, c);
                        }else{
                            for (hit2 = hitList2; hit2 !=NULL; hit2 = hit2->next) {
                                struct cpgC *cg2 = (struct cpgC *) hit2->val;
                                (cg2->mc)++;
                            }
                        }
                    } else {
                        struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                        c->c = 0;
                        c->mc = 1;
                        c->umc = 0;
                        c->strand = strand;
                        int size = hashIntValDefault(chrHash, chrom, 0);
                        if (size == 0) {
                            continue;
                        }
                        struct binKeeper *bk2 = binKeeperNew(0, size);
                        binKeeperAdd(bk2, start+i, start+i+1, c);
                        hashAdd(chgHash, chrom, bk2);
                    }
                    */
                    struct hashEl *hel2 = hashLookup(chgHash, chrom);
                    if (hel2 != NULL) {
                        struct hash *hash2 = (struct hash *) hel2->val;
                        struct hashEl *hel3 = hashLookup(hash2, key);
                        if (hel3 == NULL){
                            struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                            c->c = 0;
                            c->mc = 1;
                            c->umc = 0;
                            c->strand = strand;
                            hashAdd(hash2, key, c);
                        }else{
                            struct cpgC *cg2 = (struct cpgC *) hel3->val;
                            (cg2->mc)++;
                        }
                    } else {
                        struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                        c->c = 0;
                        c->mc = 1;
                        c->umc = 0;
                        c->strand = strand;
                        struct hash *hash = newHash(0);
                        hashAdd(hash, key, c);
                        hashAdd(chgHash, chrom, hash);
                    }
                }
            }else if(j == 'x'){
                methyCnt[1]++;
                if (fullMode){
                    /*
                    struct hashEl *hel2 = hashLookup(chgHash, chrom);
                    if (hel2 != NULL) {
                        struct binKeeper *bk2 = (struct binKeeper *) hel2->val;
                        hitList2 = binKeeperFind(bk2, start+i, start+i+1);
                        if (hitList2 == NULL){
                            struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                            c->c = 0;
                            c->mc = 0;
                            c->umc = 1;
                            c->strand = strand;
                            binKeeperAdd(bk2, start+i, start+i+1, c);
                        }else{
                            for (hit2 = hitList2; hit2 !=NULL; hit2 = hit2->next) {
                                struct cpgC *cg2 = (struct cpgC *) hit2->val;
                                (cg2->umc)++;
                            }
                        }
                    } else {
                        struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                        c->c = 0;
                        c->mc = 0;
                        c->umc = 1;
                        c->strand = strand;
                        int size = hashIntValDefault(chrHash, chrom, 0);
                        if (size == 0) {
                            continue;
                        }
                        struct binKeeper *bk2 = binKeeperNew(0, size);
                        binKeeperAdd(bk2, start+i, start+i+1, c);
                        hashAdd(chgHash, chrom, bk2);
                    }
                    */
                    struct hashEl *hel2 = hashLookup(chgHash, chrom);
                    if (hel2 != NULL) {
                        struct hash *hash2 = (struct hash *) hel2->val;
                        struct hashEl *hel3 = hashLookup(hash2, key);
                        if (hel3 == NULL){
                            struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                            c->c = 0;
                            c->mc = 0;
                            c->umc = 1;
                            c->strand = strand;
                            hashAdd(hash2, key, c);
                        }else{
                            struct cpgC *cg2 = (struct cpgC *) hel3->val;
                            (cg2->umc)++;
                        }
                    } else {
                        struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                        c->c = 0;
                        c->mc = 0;
                        c->umc = 1;
                        c->strand = strand;
                        struct hash *hash = newHash(0);
                        hashAdd(hash, key, c);
                        hashAdd(chgHash, chrom, hash);
                    }
                }
            }else if(j == 'H'){
                methyCnt[2]++;
                if (fullMode){
                    /*
                    struct hashEl *hel2 = hashLookup(chhHash, chrom);
                    if (hel2 != NULL) {
                        struct binKeeper *bk2 = (struct binKeeper *) hel2->val;
                        hitList2 = binKeeperFind(bk2, start+i, start+i+1);
                        if (hitList2 == NULL){
                            struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                            c->c = 0;
                            c->mc = 1;
                            c->umc = 0;
                            c->strand = strand;
                            binKeeperAdd(bk2, start+i, start+i+1, c);
                        }else{
                            for (hit2 = hitList2; hit2 !=NULL; hit2 = hit2->next) {
                                struct cpgC *cg2 = (struct cpgC *) hit2->val;
                                (cg2->mc)++;
                            }
                        }
                    } else {
                        struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                        c->c = 0;
                        c->mc = 1;
                        c->umc = 0;
                        c->strand = strand;
                        int size = hashIntValDefault(chrHash, chrom, 0);
                        if (size == 0) {
                            continue;
                        }
                        struct binKeeper *bk2 = binKeeperNew(0, size);
                        binKeeperAdd(bk2, start+i, start+i+1, c);
                        hashAdd(chhHash, chrom, bk2);
                    }
                    */
                    struct hashEl *hel2 = hashLookup(chhHash, chrom);
                    if (hel2 != NULL) {
                        struct hash *hash2 = (struct hash *) hel2->val;
                        struct hashEl *hel3 = hashLookup(hash2, key);
                        if (hel3 == NULL){
                            struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                            c->c = 0;
                            c->mc = 1;
                            c->umc = 0;
                            c->strand = strand;
                            hashAdd(hash2, key, c);
                        }else{
                            struct cpgC *cg2 = (struct cpgC *) hel3->val;
                            (cg2->mc)++;
                        }
                    } else {
                        struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                        c->c = 0;
                        c->mc = 1;
                        c->umc = 0;
                        c->strand = strand;
                        struct hash *hash = newHash(0);
                        hashAdd(hash, key, c);
                        hashAdd(chhHash, chrom, hash);
                    }
                }
            }else if(j == 'h'){
                methyCnt[3]++;
                if (fullMode){
                    /*
                    struct hashEl *hel2 = hashLookup(chhHash, chrom);
                    if (hel2 != NULL) {
                        struct binKeeper *bk2 = (struct binKeeper *) hel2->val;
                        hitList2 = binKeeperFind(bk2, start+i, start+i+1);
                        if (hitList2 == NULL){
                            struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                            c->c = 0;
                            c->mc = 0;
                            c->umc = 1;
                            c->strand = strand;
                            binKeeperAdd(bk2, start+i, start+i+1, c);
                        }else{
                            for (hit2 = hitList2; hit2 !=NULL; hit2 = hit2->next) {
                                struct cpgC *cg2 = (struct cpgC *) hit2->val;
                                (cg2->umc)++;
                            }
                        }
                    } else {
                        struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                        c->c = 0;
                        c->mc = 0;
                        c->umc = 1;
                        c->strand = strand;
                        int size = hashIntValDefault(chrHash, chrom, 0);
                        if (size == 0) {
                            continue;
                        }
                        struct binKeeper *bk2 = binKeeperNew(0, size);
                        binKeeperAdd(bk2, start+i, start+i+1, c);
                        hashAdd(chhHash, chrom, bk2);
                    }
                    */
                    struct hashEl *hel2 = hashLookup(chhHash, chrom);
                    if (hel2 != NULL) {
                        struct hash *hash2 = (struct hash *) hel2->val;
                        struct hashEl *hel3 = hashLookup(hash2, key);
                        if (hel3 == NULL){
                            struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                            c->c = 0;
                            c->mc = 0;
                            c->umc = 1;
                            c->strand = strand;
                            hashAdd(hash2, key, c);
                        }else{
                            struct cpgC *cg2 = (struct cpgC *) hel3->val;
                            (cg2->umc)++;
                        }
                    } else {
                        struct cpgC *c = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
                        c->c = 0;
                        c->mc = 0;
                        c->umc = 1;
                        c->strand = strand;
                        struct hash *hash = newHash(0);
                        hashAdd(hash, key, c);
                        hashAdd(chhHash, chrom, hash);
                    }
                }
            }else if(j == 'U'){
                methyCnt[6]++;
            }else if(j == 'u'){
                methyCnt[7]++;
            }
        }
    }
}

unsigned long long int *bismarkBamParse(char *samfile, struct hash *chrHash, struct hash *cpgHash, struct hash *chgHash, struct hash *chhHash, char *forwardread, char *reverseread, int isSam, int addChr, int fullMode, unsigned int iSize) {
    /*
  ### . for bases not involving cytosines                       ###
  ### X for methylated C in CHG context (was protected)         ###
  ### x for not methylated C in CHG context (was converted)     ###
  ### H for methylated C in CHH context (was protected)         ###
  ### h for not methylated C in CHH context (was converted)     ###
  ### Z for methylated C in CpG context (was protected)         ###
  ### z for not methylated C in CpG context (was converted)     ###
  ### U for methylated C in Unknown context (was protected)     ###
  ### u for not methylated C in Unknown context (was converted) ###

                default                          old_flag
           ===================              ===================
           Read 1       Read 2              Read 1       Read 2
  OT:         99          147                  67          131
  OB:         83          163                 115          179
  CTOT:       99          147                  67          131
  CTOB:       83          163                 115          179

  TODO: currently works only for bismark with bowtie1
  FIXME: lack of process of CIGAR in order to support bowtie2
    */
    char chr[100], strand, read_cove[4], genome_cove[4], methycall[1000], *row[100], cstrand;
    //char key[100];
    int fi, start, left, right, distance=0; //cutoff used for remove PCR duplication, single end as 1, paired end as 2
    //int fstart, fend, fstrand, cutoff = 0;
    unsigned long long int linecnt = 0, dupCount = 0, failCount = 0, totalbase = 0, oversizeCount = 0;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 13);
    unsigned long long int *methyCnt = malloc(sizeof(unsigned long long int) * 8);
    FILE *forward_f = NULL, *reverse_f = NULL;
    int i, cend;
    struct hash *nochr = newHash(0);
    for(i = 0; i < 8; i++) cnt[i] = 0;
    for(i = 0; i < 8; i++) methyCnt[i] = 0;
    if (fullMode){
        forward_f = mustOpen(forwardread, "w");
        reverse_f = mustOpen(reverseread, "w");
    }
    //process sam/bam list
    int numFields = chopByChar(samfile, ',', row, ArraySize(row));
    for(fi = 0; fi < numFields; fi++){
        fprintf(stderr, "\n* Processing %s\n", row[fi]);
        samfile_t *samfp;
        bam1_t *b;
        bam_header_t *h;
        struct hash *dup = newHash(0);
        if (isSam) {
            if ( (samfp = samopen(row[fi], "r", 0)) == 0) {
                fprintf(stderr, "Fail to open SAM file %s\n", samfile);
                errAbort("Error\n");
            }
        } else {
            if ( (samfp = samopen(row[fi], "rb", 0)) == 0) {
                fprintf(stderr, "Fail to open BAM file %s\n", samfile);
                errAbort("Error\n");
            }
        }
        h = samfp->header;
        b = bam_init1();
        while ( samread(samfp, b) >= 0) {
            linecnt++;
            totalbase += b->core.l_qseq;
            if ((linecnt % 10000) == 0)
                fprintf(stderr, "\r* Processed lines: %llu", linecnt);
            //change chr name to chr1, chr2 ...
            strcpy(chr, h->target_name[b->core.tid]);
            if (addChr){
                if (startsWith("GL", h->target_name[b->core.tid])) {
                    continue;
                } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                    strcpy(chr,"chrM");
                } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                    strcpy(chr, "chr");
                    strcat(chr, h->target_name[b->core.tid]);
                }
            }
            //if (sameWord(chr, "chrM")){
            //    continue; //skip chrM
            //}
            //strand
            //check Ref reads mapped to existed in chromosome size file or not
            struct hashEl *he = hashLookup(nochr, chr);
            if (he != NULL)
                continue;
            cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
            if (cend == 1){
                hashAddInt(nochr, chr, 1);
                warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
                continue;
            }
            strcpy(read_cove, bam_aux2Z(bam_aux_get(b, "XR")));
            strcpy(genome_cove, bam_aux2Z(bam_aux_get(b, "XG")));
            if (sameWord( genome_cove, read_cove )){
                strand = '+';
            }else{
                strand = '-';
            }
            //read strand
            cstrand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if ( (b->core.flag == 0) || (b->core.flag == 16)){
                //single end
                //cutoff = 1;
                start = (int) b->core.pos;
                //fend = (int) b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                //fstart = start;
                //fstrand = strand;
                left = 0;
                right = b->core.l_qseq;
                if (fullMode){
                    if (strand == '+'){
                        //fprintf(forward_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, start+(b->core.l_qseq), bam1_qname(b), b->core.qual, strand);
                        fprintf(forward_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, start+(b->core.l_qseq), "N", b->core.qual, strand);
                    }else{
                        //fprintf(reverse_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, start+(b->core.l_qseq), bam1_qname(b), b->core.qual, strand);
                        fprintf(reverse_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, start+(b->core.l_qseq), "N", b->core.qual, strand);
                    }
                }
            }else if((b->core.flag ==99) || (b->core.flag == 147) || (b->core.flag == 83) || (b->core.flag == 163) || (b->core.flag == 67) || (b->core.flag == 131) || (b->core.flag == 115) || (b->core.flag == 179) ){
                if (abs(b->core.isize) > iSize || b->core.isize == 0){
                    oversizeCount++;
                    continue;
                }
                // paired end, both new and old flag
                //cutoff = 2;
                //if (strand == '+'){
                //    fstart = (int) b->core.pos;
                //}else{
                //    fstart = (int) b->core.mpos;
                //}
                start = (int) b->core.pos;
                if (fullMode){
                    if (strand == '+'){
                        //fprintf(forward_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, start+(b->core.l_qseq), bam1_qname(b), b->core.qual, strand);
                        fprintf(forward_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, start+(b->core.l_qseq), "N", b->core.qual, strand);
                    }else{
                        //fprintf(reverse_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, start+(b->core.l_qseq), bam1_qname(b), b->core.qual, strand);
                        fprintf(reverse_f, "%s\t%u\t%u\t%s\t%i\t%c\n", chr, start, start+(b->core.l_qseq), "N", b->core.qual, strand);
                    }
                }
                //fend = (int) fstart + abs(b->core.isize);
                if( (b->core.flag == 99) || (b->core.flag == 83) || (b->core.flag == 67) || (b->core.flag == 115)) {
                    //first pair
                    left = 0;
                    right = b->core.l_qseq;
                    //if (strand == '+'){
                        //fstrand = '+';
                    //}else{
                        //fstrand = '-';
                    //}
                }else{
                    //second pair
                    //need take care about the overlap of 1st pair
                    //the overlap region should only be called once
                    /*  
                     *  R1----------------->
                     *          R2<------------------
                     *
                     *  R2----------------->
                     *          R1<----------------
                     *
                     *
                    */
                    distance = abs(b->core.pos - b->core.mpos);
                    if (strand == '+'){
                        left = 0;
                        right = min(distance, b->core.l_qseq); //might not overlap
                        //fstrand = '-';
                    } else {
                        right = b->core.l_qseq;
                        left = max(0, b->core.l_qseq - distance); //might not overlap
                        //fstrand = '+';
                    }
                }
            } else {
                //other type of reads, might be fail of quality check
                failCount++;
                continue;
            }
            //remove dup
            //if (sprintf(key, "%s:%i:%i:%c:%s:%s", chr, fstart, fend, strand, read_cove, genome_cove) < 0)
            //if (sprintf(key, "%s:%i:%i:%c", chr, fstart, fend, fstrand) < 0)
            //    errAbort("Mem ERROR");
            //hashIncInt(dup, key);
            //fprintf(stderr, "key %s added\n", key);
            //int judge = hashIntVal(dup, key);
            //fprintf(stderr, "judge %i cutoff %i\n", judge, cutoff);
            //if (judge > cutoff){
            //    dupCount++;
            //    continue;
            //}
            //process methylation call
            //if(strand == '+'){
            //    strcpy(methycall, bam_aux2Z(bam_aux_get(b, "XM")));
            //}else{
            //    strcpy(methycall, strrev(bam_aux2Z(bam_aux_get(b, "XM"))));
            //}
            /* 
             * the bismark code have reverse the methylation call string when mapping to - strand
             * when I also did this, get many position not belong to CpG C
             * just used the methylation string no matter strand even works, why? FIXME
            */
            strcpy(methycall, bam_aux2Z(bam_aux_get(b, "XM")));
            assignCpGcount(chrHash, cpgHash, chgHash, chhHash, chr, start, methycall, cstrand, left, right, methyCnt, fullMode);
            //fprintf(stdout, "%s\t%i\t%i\t%i\t%c\t%i\t%i\t%s\n", chr, start, fstart, fend, strand, left, right, methycall);
        }
        samclose(samfp);
        bam_destroy1(b);
        //bam_header_destroy(h);
        freeHash(&dup);
    }
    fprintf(stderr, "\r* Processed lines: %llu\n", linecnt);
    fprintf(stderr, "* Oversized alignments: %llu\n", oversizeCount);
    fprintf(stderr, "* Quality Failed alignments: %llu\n", failCount);
    //fprintf(stderr, "* Duplicated alignments: %lu\n", dupCount);
    freeHash(&nochr);
    cnt[0] = linecnt;
    cnt[1] = failCount;
    cnt[2] = dupCount;
    cnt[3] = totalbase;
    cnt[4] = oversizeCount;
    for (i=5;i<13;i++) cnt[i] = methyCnt[i-5];
    if(fullMode){
        carefulClose(&forward_f);
        carefulClose(&reverse_f);
    }
    return cnt;
}

unsigned long long int *sam2bed(char *samfile, char *outbed, struct hash *chrHash, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat) {
    samfile_t *samfp;
    FILE *outbed_f = mustOpen(outbed, "w");
    char chr[100], key[100], strand;
    unsigned int start, end, cend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 10);
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    unsigned long long int map_supp = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    if (isSam) {
        if ( (samfp = samopen(samfile, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", samfile);
            errAbort("Error\n");
        }
    } else {
        if ( (samfp = samopen(samfile, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", samfile);
            errAbort("Error\n");
        }
    }
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    int8_t *buf;
    int max_buf;
    buf = 0;
    max_buf = 0;
    uint8_t *seq;
    while ( samread(samfp, b) >= 0) {
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1++;
            }else{
                if(treat)
                    read_end1++;
                else
                    read_end2++;
            }
        }else{
            read_end1++;
        }
        if (((read_end1 + read_end2) % 10000) == 0)
            fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
        if (b->core.flag & BAM_FUNMAP)
            continue;
        if (b->core.flag & BAM_SUPP){
            map_supp++;
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_mapped++;
            }else{
                if (treat)
                    read_end1_mapped++;
                else
                    read_end2_mapped++;
            }
        }else{
            read_end1_mapped++;
        }
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        //check Ref reads mapped to existed in chromosome size file or not
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_used++;
            }else{
                if (treat)
                    read_end1_used++;
                else
                    read_end2_used++;
            }
        }else{
            read_end1_used++;
        }
        //get mapping location for paired-end or single-end
        if (treat){
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                    //start = max(end - extension, 0);
                }
            }

        }else{
        if (b->core.flag & BAM_FPAIRED) {
            if (!(b->core.flag & BAM_FMUNMAP)){
                if (b->core.flag & BAM_FREAD1){
                    if (abs(b->core.isize) > iSize || b->core.isize == 0){
                        continue;
                    }else{
                        reads_mapped++;
                        if (b->core.qual >= mapQ)
                            reads_mapped_unique++;
                        if (b->core.isize > 0){
                            start = (unsigned int) b->core.pos;
                            strand = '+';
                            int tmpend = start + b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }else{
                            start = (unsigned int) b->core.mpos;
                            strand = '-';
                            int tmpend = start - b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }
                
                    }
                }else{
                    continue;
                }
            }else{
                if (discardWrongEnd){
                    continue;
                }else{
                    reads_mapped++;
                    if (b->core.qual >= mapQ)
                        reads_mapped_unique++;
                    start = (unsigned int) b->core.pos;
                    int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                    end = min(cend, (unsigned int)tmpend);
                    strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                    if (extension) {
                        if (strand == '+'){
                            end = min(start + extension, cend);
                        }else{
                            if (end < extension)
                                start = 0;
                            else
                                start = end - extension;
                            //start = max(end - extension, 0);
                        }
                    }
                }
            }
        }else{
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                }
            }
        }
    }
        //remove dup or not
        if (rmDup){
            if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                errAbort("Mem ERROR");
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        reads_nonredundant++;
        if (b->core.qual >= mapQ)
            reads_nonredundant_unique++;
        //output bed
        int i, qlen = b->core.l_qseq;
        if(b->core.qual >= mapQ){
            fprintf(outbed_f, "%s\t%u\t%u\t", chr, start, end);
            //print read sequence
            if (max_buf < qlen + 1 ) {
                max_buf = qlen + 1;
                kroundup32(max_buf);
                buf = realloc(buf, max_buf);
            }
            buf[qlen] = 0;
            seq = bam1_seq(b);
            for (i = 0; i < qlen; ++i)
                buf[i] = bam1_seqi(seq, i);
            if (b->core.flag & 16) {
                for (i = 0; i < qlen>>1; ++i){
                    int8_t t = seq_comp_table[buf[qlen - 1 - i]];
                    buf[qlen - 1 - i] = seq_comp_table[buf[i]];
                    buf[i] = t;
                }
                if (qlen&1) buf[i] = seq_comp_table[buf[i]];
            }
            for (i = 0; i < qlen; ++i)
                buf[i] = bam_nt16_rev_table[buf[i]];
            fprintf(outbed_f, "%s", (char*)buf);

            fprintf(outbed_f, "\t%i\t%c\n", b->core.qual, strand);
        }
    }
    fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
    fprintf(stderr, "* Skipped supplementary alignments: %llu\n", map_supp);
    samclose(samfp);
    free(buf);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    carefulClose(&outbed_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_nonredundant_unique;
    return cnt;
}

unsigned long long int *sam2bedwithCpGstat(char *samfile, char *outbed, struct hash *chrHash, struct hash *cpgHash, struct slInt **cpgCount, struct slInt **slPair, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat) {
    samfile_t *samfp;
    FILE *outbed_f = mustOpen(outbed, "w");
    struct slInt *countcpg = NULL;
    struct slInt *pairsl = NULL;
    char chr[100], key[100], strand;
    unsigned int start, end, cend;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 10);
    unsigned long long int read_end1 = 0, read_end2 = 0;
    unsigned long long int read_end1_mapped = 0, read_end2_mapped = 0;
    unsigned long long int read_end1_used = 0, read_end2_used = 0;
    unsigned long long int reads_nonredundant = 0;
    unsigned long long int reads_nonredundant_unique = 0;
    unsigned long long int reads_mapped = 0;
    unsigned long long int reads_mapped_unique = 0;
    unsigned long long int map_supp = 0;
    struct hash *nochr = newHash(0), *dup = newHash(0);
    if (isSam) {
        if ( (samfp = samopen(samfile, "r", 0)) == 0) {
            fprintf(stderr, "Fail to open SAM file %s\n", samfile);
            errAbort("Error\n");
        }
    } else {
        if ( (samfp = samopen(samfile, "rb", 0)) == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", samfile);
            errAbort("Error\n");
        }
    }
    bam1_t *b;
    bam_header_t *h;
    h = samfp->header;
    b = bam_init1();
    int8_t *buf;
    int max_buf;
    buf = 0;
    max_buf = 0;
    uint8_t *seq;
    while ( samread(samfp, b) >= 0) {
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1++;
            }else{
                if(treat)
                    read_end1++;
                else
                    read_end2++;
            }
        }else{
            read_end1++;
        }
        if (((read_end1 + read_end2) % 10000) == 0)
            fprintf(stderr, "\r* Processed read ends: %llu", (read_end1 + read_end2));
        if (b->core.flag & BAM_FUNMAP)
            continue;
        if (b->core.flag & BAM_SUPP){
            map_supp++;
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_mapped++;
            }else{
                if (treat)
                    read_end1_mapped++;
                else
                    read_end2_mapped++;
            }
        }else{
            read_end1_mapped++;
        }
        //change chr name to chr1, chr2 ...
        strcpy(chr, h->target_name[b->core.tid]);
        if (addChr){
            if (startsWith("GL", h->target_name[b->core.tid])) {
                continue;
            } else if (sameWord(h->target_name[b->core.tid], "MT")) {
                strcpy(chr,"chrM");
            } else if (!startsWith("chr", h->target_name[b->core.tid])) {
                strcpy(chr, "chr");
                strcat(chr, h->target_name[b->core.tid]);
            }
        }
        //check Ref reads mapped to existed in chromosome size file or not
        struct hashEl *he = hashLookup(nochr, chr);
        if (he != NULL)
            continue;
        cend = (unsigned int) (hashIntValDefault(chrHash, chr, 2) - 1);
        if (cend == 1){
            hashAddInt(nochr, chr, 1);
            warn("* Warning: read ends mapped to chromosome %s will be discarded as %s not existed in the chromosome size file", chr, chr);
            continue;
        }
        if (b->core.flag & BAM_FPAIRED) {
            if (b->core.flag & BAM_FREAD1){
                read_end1_used++;
            }else{
                if (treat)
                    read_end1_used++;
                else
                    read_end2_used++;
            }
        }else{
            read_end1_used++;
        }
        //get mapping location for paired-end or single-end
        if (treat){
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                    //start = max(end - extension, 0);
                }
            }

        }else{
        if (b->core.flag & BAM_FPAIRED) {
            if (!(b->core.flag & BAM_FMUNMAP)){
                if (b->core.flag & BAM_FREAD1){
                    if (abs(b->core.isize) > iSize || b->core.isize == 0){
                        continue;
                    }else{
                        reads_mapped++;
                        if (b->core.qual >= mapQ)
                            reads_mapped_unique++;
                        if (b->core.isize > 0){
                            start = (unsigned int) b->core.pos;
                            strand = '+';
                            int tmpend = start + b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }else{
                            start = (unsigned int) b->core.mpos;
                            strand = '-';
                            int tmpend = start - b->core.isize;
                            end = min(cend, (unsigned int)tmpend);
                        }
                
                    }
                }else{
                    continue;
                }
            }else{
                if (discardWrongEnd){
                    continue;
                }else{
                    reads_mapped++;
                    if (b->core.qual >= mapQ)
                        reads_mapped_unique++;
                    start = (unsigned int) b->core.pos;
                    int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
                    end = min(cend, (unsigned int)tmpend);
                    strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
                    if (extension) {
                        if (strand == '+'){
                            end = min(start + extension, cend);
                        }else{
                            if (end < extension)
                                start = 0;
                            else
                                start = end - extension;
                            //start = max(end - extension, 0);
                        }
                    }
                }
            }
        }else{
            reads_mapped++;
            if (b->core.qual >= mapQ)
                reads_mapped_unique++;
            start = (unsigned int) b->core.pos;
            int tmpend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + b->core.l_qseq;
            end = min(cend, (unsigned int)tmpend);
            strand = (b->core.flag&BAM_FREVERSE)? '-' : '+';
            if (extension) {
                if (strand == '+'){
                    end = min(start + extension, cend);
                }else{
                    if (end < extension)
                        start = 0;
                    else
                        start = end - extension;
                }
            }
        }
    }
        //remove dup or not
        if (rmDup){
            if (sprintf(key, "%s:%u:%u:%c", chr, start, end, strand) < 0)
                errAbort("Mem ERROR");
            struct hashEl *hel = hashLookup(dup, key);
            if (hel == NULL) {
                hashAddInt(dup, key, 1);
            } else {
                continue;
            }
        }
        reads_nonredundant++;
        if (b->core.qual >= mapQ)
            reads_nonredundant_unique++;
        //output bed
        int i, qlen = b->core.l_qseq;
        if(b->core.qual >= mapQ){
            fprintf(outbed_f, "%s\t%u\t%u\t", chr, start, end);
            //print read sequence
            if (max_buf < qlen + 1 ) {
                max_buf = qlen + 1;
                kroundup32(max_buf);
                buf = realloc(buf, max_buf);
            }
            buf[qlen] = 0;
            seq = bam1_seq(b);
            for (i = 0; i < qlen; ++i)
                buf[i] = bam1_seqi(seq, i);
            if (b->core.flag & 16) {
                for (i = 0; i < qlen>>1; ++i){
                    int8_t t = seq_comp_table[buf[qlen - 1 - i]];
                    buf[qlen - 1 - i] = seq_comp_table[buf[i]];
                    buf[i] = t;
                }
                if (qlen&1) buf[i] = seq_comp_table[buf[i]];
            }
            for (i = 0; i < qlen; ++i)
                buf[i] = bam_nt16_rev_table[buf[i]];
            fprintf(outbed_f, "%s", (char*)buf);

            fprintf(outbed_f, "\t%i\t%c\n", b->core.qual, strand);
            //output insert size
            slAddHead(&pairsl, slIntNew(end-start));
            //do cpg stat
            struct hashEl *hel5 = hashLookup(cpgHash, chr);
            if (hel5 == NULL)
                continue;
            struct binKeeper *bs5 = (struct binKeeper *) hel5->val;
            int c = binKeeperCpGstat(bs5, start, end);
            slAddHead(&countcpg, slIntNew(c));
        }
    }
    fprintf(stderr, "\r* Processed read ends: %llu\n", (read_end1 + read_end2));
    slReverse(&countcpg);
    *cpgCount = countcpg;
    if (read_end1 == read_end2) { //paired-end
        fprintf(stderr, "* Paired end data\n");
        slReverse(&pairsl);
        *slPair = pairsl;
        //fprintf(stderr, "* %d elements\n", slCount(slPair));
    }else{
        if (read_end2 == 0){
            fprintf(stderr, "* Single end data\n");
        } else {
            fprintf(stderr, "* Mixed of single & paired end data\n");
        }
    }
    fprintf(stderr, "* Skipped supplementary alignments: %llu\n", map_supp);
    samclose(samfp);
    free(buf);
    bam_destroy1(b);
    freeHash(&nochr);
    freeHash(&dup);
    carefulClose(&outbed_f);
    cnt[0] = read_end1;
    cnt[1] = read_end2;
    cnt[2] = read_end1_mapped;
    cnt[3] = read_end2_mapped;
    cnt[4] = read_end1_used;
    cnt[5] = read_end2_used;
    cnt[6] = reads_mapped;
    cnt[7] = reads_mapped_unique;
    cnt[8] = reads_nonredundant;
    cnt[9] = reads_nonredundant_unique;
    return cnt;
}

struct hash *initGenomeCovHash(struct hash *chrHash){
    struct hash *hash = newHash(0);
    struct hashEl *hel;
    int i;
    struct hashCookie cookie = hashFirst(chrHash);
    while((hel = hashNext(&cookie)) != NULL) {
        struct slInt *sl = NULL;
        int s = ptToInt(hel->val);
        for (i=0; i < s; i++){
            slAddHead(&sl, slIntNew(0));
        }
        slReverse(&sl);
        hashAdd(hash, hel->name, sl);
    }
    return hash;
}

struct hash *MREfrag2Hash (char *fragfile, int minlen, int maxlen){
    struct hash *hash = newHash(0);
    char *row[20], *line;
    int start, end, length, rstart, rend;
    char key1[100], key2[100];
    struct lineFile *frag_stream = lineFileOpen2(fragfile, TRUE);
    while ( lineFileNextReal(frag_stream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 4)
            errAbort("file %s doesn't appear to be in bed format. At least 3 fields required, got %d", fragfile, numFields);
        start = (int) strtol(row[1], NULL, 0);
        end = (int) strtol(row[2], NULL, 0);
        length  = end - start;
        if (length >= minlen && length <= maxlen){
            struct mreFrag *mre1 = malloc(sizeof(struct mreFrag));
            struct mreFrag *mre2 = malloc(sizeof(struct mreFrag));
            if (sameWord(row[3], "CGCG")) { //special case for CGCG
                rstart = start + 2;
                rend = end - 2;
                if (sprintf(key1, "%s:%i:%c", row[0], rstart, '+') < 0)
                    errAbort("Mem ERROR");
                if (sprintf(key2, "%s:%i:%c", row[0], rend, '-') < 0)
                    errAbort("Mem ERROR");
            } else { // case for CCGG, CCGC, GCGC, ACGT
                rstart = start + 1;
                rend = end - 1;
                if (sprintf(key1, "%s:%i:%c", row[0], rstart, '+') < 0)
                    errAbort("Mem ERROR");
                if (sprintf(key2, "%s:%i:%c", row[0], rend, '-') < 0)
                    errAbort("Mem ERROR");
            }
            mre1->head = 1; mre1->start = rstart; mre1->end = rend;
            mre2->head= 0; mre2->start = rstart; mre2->end = rend;
            mre1->reads_count = 0; mre2->reads_count = 0;
            strcpy(mre1->pair, key2); strcpy(mre2->pair, key1);
            strcpy(mre1->site, row[3]);
            strcpy(mre2->site, row[3]);
            strcpy(mre1->chr, row[0]);
            strcpy(mre2->chr, row[0]);
            hashAdd(hash, key1, mre1);
            hashAdd(hash, key2, mre2);
        } else {
            continue;
        }
    }
    lineFileClose(&frag_stream);
    return hash;
}

struct hash *cpgBed2BinKeeperHash (struct hash *chrHash, char *cpgbedfile){
    struct hash *hash = newHash(0);
    char *row[20], *line;
    int start, end;
    struct lineFile *stream = lineFileOpen2(cpgbedfile, TRUE);
    while ( lineFileNextReal(stream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 4)
            errAbort("file %s doesn't appear to be in bed format. At least 3 fields required, got %d", cpgbedfile, numFields);
        start = (int) strtol(row[1], NULL, 0);
        end = (int) strtol(row[2], NULL, 0);
        struct cpgC *c = malloc(sizeof(struct cpgC));
        c->c = 0;
        c->mc = 0;
        c->umc = 0;
        c->strand = '.';
        struct hashEl *hel = hashLookup(hash, row[0]);
        if (hel != NULL) {
            struct binKeeper *bk = (struct binKeeper *) hel->val;
            binKeeperAdd(bk, start, end, c); // c stores cpg coverage count
        } else {
            int size = hashIntValDefault(chrHash, row[0], 0);
            if (size == 0) {
                continue;
            }
            struct binKeeper *bk = binKeeperNew(0, size);
            binKeeperAdd(bk, start, end, c);
            hashAdd(hash, row[0], bk);
        }
    }
    lineFileClose(&stream);
    return hash;
}

struct hash *cpgBed2BinKeeperHashBismark (struct hash *chrHash, char *cpgbedfile){
    struct hash *hash = newHash(0);
    char *row[20], *line;
    int start, end;
    struct lineFile *stream = lineFileOpen2(cpgbedfile, TRUE);
    while ( lineFileNextReal(stream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 4)
            errAbort("file %s doesn't appear to be in bed format. At least 3 fields required, got %d", cpgbedfile, numFields);
        start = (int) strtol(row[1], NULL, 0);
        end = (int) strtol(row[2], NULL, 0);
        struct cpgC *c1 = malloc(sizeof(struct cpgC));
        struct cpgC *c2 = malloc(sizeof(struct cpgC)); //change from each CpG to 2 base as bismark does this
        c1->c = 0;
        c1->mc = 0;
        c1->umc = 0;
        c1->strand = '+';
        c2->c = 0;
        c2->mc = 0;
        c2->umc = 0;
        c2->strand = '-';
        struct hashEl *hel = hashLookup(hash, row[0]);
        if (hel != NULL) {
            struct binKeeper *bk = (struct binKeeper *) hel->val;
            binKeeperAdd(bk, start, start+1, c1); // c stores cpg coverage count
            binKeeperAdd(bk, start+1, end, c2); // c stores cpg coverage count
        } else {
            int size = hashIntValDefault(chrHash, row[0], 0);
            if (size == 0) {
                continue;
            }
            struct binKeeper *bk = binKeeperNew(0, size);
            binKeeperAdd(bk, start, start+1, c1);
            binKeeperAdd(bk, start+1, end, c2);
            hashAdd(hash, row[0], bk);
        }
    }
    lineFileClose(&stream);
    return hash;
}

long long * bedCpGstat(struct hash *cpgHash, char *bedfile){
    char *row[20], *line;
    int start, end;
    long long *cnt = malloc(sizeof(long long)*2);
    long long c1 = 0, c2 = 0;
    struct lineFile *stream = lineFileOpen2(bedfile, TRUE);
    while ( lineFileNextReal(stream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 4)
            errAbort("file %s doesn't appear to be in bed format. At least 3 fields required, got %d", bedfile, numFields);
        start = (int) strtol(row[1], NULL, 0);
        end = (int) strtol(row[2], NULL, 0);
        //output insert size
        c1 += (long long)(end - start);
        //do cpg stat
        struct hashEl *hel5 = hashLookup(cpgHash, row[0]);
        if (hel5 == NULL)
            continue;
        struct binKeeper *bs5 = (struct binKeeper *) hel5->val;
        int c = binKeeperCpGstat(bs5, start, end);
        c2 += (long long) c;
    }
    lineFileClose(&stream);
    //fprintf(stderr, "c1 : %lli", c1);
    //fprintf(stderr, "c2 : %lli", c2);
    cnt[0] = c1;
    cnt[1] = c2;
    return cnt;
}

unsigned long long int *filterReadByMREsite(struct hash *hash, char *inBed, char *outBed, int call, char *prefix, int guess){
    FILE *f  = mustOpen(outBed, "w");
    //if (strcmp(prefix, "NULL") != 0){
        char *out1 = malloc(200);
        char *out2 = malloc(200);
        char *out3 = malloc(200);
        char *out4 = malloc(200);
        char *out5 = malloc(200);
        char *out6 = malloc(200);
        if (sprintf(out1, "%s_CCGG_reads.txt", prefix) < 0)
            errAbort("Mem error");
        if (sprintf(out2, "%s_CCGC_reads.txt", prefix) < 0)
            errAbort("Mem error");
        if (sprintf(out3, "%s_GCGC_reads.txt", prefix) < 0)
            errAbort("Mem error");
        if (sprintf(out4, "%s_ACGT_reads.txt", prefix) < 0)
            errAbort("Mem error");
        if (sprintf(out5, "%s_CGCG_reads.txt", prefix) < 0)
            errAbort("Mem error");
        if (sprintf(out6, "%s_unknown_reads.txt", prefix) < 0)
            errAbort("Mem error");
        FILE *f1 = mustOpen(out1, "w");
        FILE *f2 = mustOpen(out2, "w");
        FILE *f3 = mustOpen(out3, "w");
        FILE *f4 = mustOpen(out4, "w");
        FILE *f5 = mustOpen(out5, "w");
        FILE *f6 = mustOpen(out6, "w");
    //}
    char strand, *row[20], *line;
    struct lineFile *inBedStream = lineFileOpen2(inBed, TRUE);
    char key1[100]; //key2[100];
    char key11[100]; //guess MRE site with read position 1 or 4
    int start, end, rstart, rend;
    struct hashEl *hel;
    unsigned long long int CCGG = 0, CCGC = 0, GCGC = 0, ACGT = 0, CGCG = 0, unknown = 0;
    unsigned long long int *cnt = malloc(sizeof(unsigned long long int) * 6);
    while( lineFileNextReal(inBedStream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 6)
            errAbort("file %s doesn't appear to be in bed format. At least 6 fields required, got %d", inBed, numFields);
        strand = row[5][0];
        start = (int)strtol(row[1], NULL, 0);
        end = (int)strtol(row[2], NULL, 0);
        if (strand == '+'){
            rstart = start - call;
            rend = end;
            if (sprintf(key1, "%s:%i:%c", row[0], rstart, '+') < 0)
                errAbort("Mem ERROR");
            if (guess){
                if (sprintf(key11, "%s:%i:%c", row[0], rstart - 3, '+') < 0)
                    errAbort("Mem ERROR");
            }
            //if (sprintf(key2, "%s:%i:%c", bed->chrom, bed->chromStart - 3, '+') < 0)
            //    errAbort("Mem ERROR");
        } else {
            rstart = start;
            rend = end + call;
            if (sprintf(key1, "%s:%i:%c", row[0], rend, '-') < 0)
                errAbort("Mem ERROR");
            if (guess){
                if (sprintf(key11, "%s:%i:%c", row[0], rend + 3, '-') < 0)
                    errAbort("Mem ERROR");
            }
            //if (sprintf(key2, "%s:%i:%c", bed->chrom, bed->chromEnd + 3, '-') < 0)
            //    errAbort("Mem ERROR");
        }
        hel = hashLookup(hash, key1);
        if (hel == NULL){
            if (guess){
                hel = hashLookup(hash, key11);
            }
        }
        //struct hashEl *hel2 = hashLookup(hash, key2);
        if (hel != NULL) {
            //bedOutputN(bed, 6, f, '\t', '\n');
            fprintf(f, "%s\t%i\t%i\t%s\t%s\t%c\n", row[0], rstart, rend, row[3], row[4], strand);
            struct mreFrag *mre = (struct mreFrag *) hel->val;
            mre->reads_count++;
            if (sameWord(mre->site, "CCGG")) {
                CCGG++;
                //if (strcmp(prefix, "NULL") != 0)
                    fprintf(f1, "%s\t%i\t%i\t%s\t%s\t%c\n", row[0], rstart, rend, row[3], row[4], strand);
            } else if (sameWord(mre->site, "CCGC") || sameWord(mre->site, "GCGG")) {
                CCGC++;
                //if (strcmp(prefix, "NULL") != 0)
                    fprintf(f2, "%s\t%i\t%i\t%s\t%s\t%c\n", row[0], rstart, rend, row[3], row[4], strand);
            } else if (sameWord(mre->site, "GCGC")){
                GCGC++;
                //if (strcmp(prefix, "NULL") != 0)
                    fprintf(f3, "%s\t%i\t%i\t%s\t%s\t%c\n", row[0], rstart, rend, row[3], row[4], strand);
            } else if (sameWord(mre->site, "ACGT")){
                ACGT++;
                //if (strcmp(prefix, "NULL") != 0)
                    fprintf(f4, "%s\t%i\t%i\t%s\t%s\t%c\n", row[0], rstart, rend, row[3], row[4], strand);
            } else if (sameWord(mre->site, "CGCG")){
                CGCG++;
                //if (strcmp(prefix, "NULL") != 0)
                    fprintf(f5, "%s\t%i\t%i\t%s\t%s\t%c\n", row[0], rstart, rend, row[3], row[4], strand);
            }
        /*} else if(hel2 != NULL){
            bedOutputN(bed, 6, f, '\t', '\n');
            struct mreFrag *mre = (struct mreFrag *) hel2->val;
            mre->reads_count++;
            if (sameWord(mre->site, "CCGG")) {
                CCGG++;
            } else if (sameWord(mre->site, "CCGC") || sameWord(mre->site, "GCGG")) {
                CCGC++;
            } else if (sameWord(mre->site, "GCGC")){
                GCGC++;
            } else if (sameWord(mre->site, "ACGT")){
                ACGT++;
            } else if (sameWord(mre->site, "CGCG")){
                CGCG++;
            }*/
        } else {
            unknown++;
            //if (strcmp(prefix, "NULL") != 0)
                fprintf(f6, "%s\t%i\t%i\t%s\t%s\t%c\n", row[0], rstart, rend, row[3], row[4], strand);
        }
    }
    carefulClose(&f);
    lineFileClose(&inBedStream);
    //if (strcmp(prefix, "NULL") != 0){
        carefulClose(&f1);
        carefulClose(&f2);
        carefulClose(&f3);
        carefulClose(&f4);
        carefulClose(&f5);
        carefulClose(&f6);
    //}
    //bedFreeList(&bedList);
    cnt[0] = CCGG;
    cnt[1] = CCGC;
    cnt[2] = GCGC;
    cnt[3] = ACGT;
    cnt[4] = CGCG;
    cnt[5] = unknown;
    return cnt;
}

double calCpGscore (struct mreFrag *mre, unsigned long long int *cnt){
    unsigned long long int norm = 0;
    if (sameWord(mre->site, "CCGG")) {
        norm = cnt[0];
    //} else if (sameWord(mre->site, "CCGC")) {
    } else if (sameWord(mre->site, "CCGC") || sameWord(mre->site, "GCGG")) {
        norm = cnt[1];
    } else if (sameWord(mre->site, "GCGC")){
        norm = cnt[2];
    } else if (sameWord(mre->site, "ACGT")){
        norm = cnt[3];
    } else if (sameWord(mre->site, "CGCG")){
        norm = cnt[4];
    }
    if (norm > 0)
        return (double) mre->reads_count * 1000000.0 / norm;
    else
        return 0;
}

unsigned long long int CpGscorebedGraph(struct hash *hash, unsigned long long int *cnt, char *outfile){
    /* when same CpG site get two score, add them -- disabled temp -- abled
     * for overlapping site, solve by using single base score - seems no this situation?*/
    FILE *f = mustOpen(outfile, "w");
    struct hashEl *hel, *hel2, *hel3; //*hel4;
    struct hashCookie cookie = hashFirst(hash);
    unsigned long long int c = 0; //count for covered CpG site
    double score = 0;
    char key[100]; //key2[100];
    struct hash *cpgHash = newHash(0);
    while ( (hel = hashNext(&cookie)) != NULL ){
        struct mreFrag *mre = (struct mreFrag *) hel->val;
        if (mre->reads_count > 0){
            score = calCpGscore(mre, cnt);
            if (score > 0){
                struct cpgScore *cgC = malloc(sizeof(struct cpgScore));
                //struct cpgScore *cgG = malloc(sizeof(struct cpgScore));
                cgC->score = score;
                //cgG->score = score;
                strcpy(cgC->chr, mre->chr);
                //strcpy(cgG->chr, mre->chr);
                if (mre->head){
                    //fprintf(f, "%s\t%i\t%i\t%.4f\n", mre->chr, mre->start, mre->start + 2, score);
                    cgC->start = mre->start;
                    //cgG->start = mre->start + 1;
                } else {
                    //fprintf(f, "%s\t%i\t%i\t%.4f\n", mre->chr, mre->end - 2, mre->end, score);
                    cgC->start = mre->end - 2;
                    //cgG->start = mre->end - 1;
                }
                if (sprintf(key, "%s:%i", cgC->chr, cgC->start) < 0)
                    errAbort("Mem ERROR");
                //if (sprintf(key2, "%s:%i", cgG->chr, cgG->start) < 0)
                //    errAbort("Mem ERROR");
                hel2 = hashLookup(cpgHash, key);
                if (hel2 == NULL){
                    hashAdd(cpgHash, key, cgC);
                    c++; //CpG site my coverd by more than one MRE site
                } else {
                    struct cpgScore *cg1 = (struct cpgScore *) hel2->val;
                    cg1->score += score;
                }
                //hel4 = hashLookup(cpgHash, key2);
                //if (hel4 == NULL){
                //    hashAdd(cpgHash, key2, cgG);
                //} else {
                //    struct cpgScore *cg4 = (struct cpgScore *) hel4->val;
                //    cg4->score += score;
                //}
            }
        }
    }
    struct hashCookie cookie2 = hashFirst(cpgHash);
    while ((hel3 = hashNext(&cookie2)) != NULL){
        struct cpgScore *cg2 = (struct cpgScore *) hel3->val;
        fprintf(f, "%s\t%i\t%i\t%.4f\n", cg2->chr, cg2->start, cg2->start + 2, cg2->score);
    }
    carefulClose(&f);
    freeHash(&cpgHash);
    return c;
}

struct fragd *fragmentStats(struct hash *hash, unsigned long long int *cnt2, unsigned int mapQ, unsigned long long int *cnt, unsigned long long int cnt1, char *output, int minlen, int maxlen, int win){
    char *outReport;
    if (asprintf(&outReport, "%s.CpG.report", output) < 0)
        errAbort("Preparing output wrong");
    FILE *f = mustOpen(outReport, "w");
    int highend = 0, lowend = 0, solosite = 0, soloreads = 0, pairsite = 0;
    struct slInt *fraglist = NULL, *frag;
    struct hashEl *hel, *hel2;
    struct hashCookie cookie = hashFirst(hash);
    while( (hel = hashNext(&cookie)) != NULL ){
        struct mreFrag *mre = (struct mreFrag *) hel->val;
        if (mre->reads_count > 0){
            if (mre->head){
                hel2 = hashLookup(hash, mre->pair);
                if (hel2 != NULL){
                    struct mreFrag *mre2 = (struct mreFrag *) hel2->val;
                    if (mre2->reads_count > 0){
                        //a valid fragment
                        pairsite++;
                        frag = slIntNew(mre->end - mre->start);
                        //fprintf(stdout, "%i\n", frag->val);
                        slAddHead(&fraglist, frag);
                        if (mre->reads_count >= mre2->reads_count){
                            highend += mre->reads_count;
                            lowend += mre2->reads_count;
                        }else{
                            highend += mre2->reads_count;
                            lowend += mre->reads_count;
                        }
                    } else{
                        //a head solo
                        soloreads += mre->reads_count;
                        solosite++;
                    }
                } else {
                    errAbort("MRE site Key error.\n");
                }
            } else {
                hel2 = hashLookup(hash, mre->pair);
                if (hel2 != NULL){
                    struct mreFrag *mre2 = (struct mreFrag *) hel2->val;
                    if (mre2->reads_count > 0){
                        //a valid fragment, but already processed
                        continue;
                    } else{
                        //a tail solo
                        soloreads += mre->reads_count;
                        solosite++;
                    }
                } else {
                    errAbort("MRE site Key error.\n");
                }
            }
        }
    }
    slReverse(&fraglist);
    //prepare histgraph
    int bins = (maxlen - minlen) / win;
    int histoMin = 0, histoMax = 0, j;
    struct range *scale = malloc(sizeof(struct range) * bins);
    for (j = 0; j < bins; j++){
        scale[j].s = minlen + j*win;
        scale[j].e = minlen + (j+1)*win;
        scale[j].histo = 0;
    }
    //scale[0]->s = 0; scale[0]->e = minlen; scale[0]->histo = 0;
    //scale[j]->s = maxlen; scale[j]->e = 0; scale[j]->histo = 0;

    for ( frag = fraglist; frag != NULL; frag = frag->next){
        if (frag->val <= minlen){
            histoMin++;
        } else if (frag->val > maxlen){
            histoMax++;
        } else{
            for (j = 0; j < bins; j++){
                if (frag->val > scale[j].s && frag->val <= scale[j].e){
                    (scale[j].histo)++;
                    break;
                }
            }
        }
    }
    struct fragd *fragdistro = malloc(sizeof(struct fragd) * (bins + 3));
    for (j=0; j<(bins+2); j++){
        fragdistro[j].label = 0;
        fragdistro[j].pct = 0;
    }
    //fragments distro plot
    char tmpRfile2[50];
    strcpy(tmpRfile2, template_name);
    int fd2 = mkstemp(tmpRfile2);
    if (fd2 == -1)
        errAbort("create temp file error.");
    FILE *fout2 = fdopen(fd2, "w");
    fprintf(fout2, "dat <- data.frame(label=rep(NA, %d), total=rep(0, %d), cov=rep(0, %d), stringsAsFactors=FALSE)\n", bins+2, bins+2, bins+2);
    
    fprintf(f, "total reads (pair): %llu\n", cnt2[0]);
    fprintf(f, "    read ends 1: %llu\n", cnt2[0]);
    fprintf(f, "    read ends 2: %llu\n", cnt2[1]);
    fprintf(f, "    mapped read ends 1: %llu\n", cnt2[2]);
    fprintf(f, "    mapped read ends 2: %llu\n", cnt2[3]);
    fprintf(f, "    used read ends 1: %llu\n", cnt2[4]);
    fprintf(f, "    used read ends 2: %llu\n", cnt2[5]);
    //fprintf(f, "non-redundant reads (pair): %llu\n\n", cnt2[8]);
    fprintf(f, "mapped reads (pair): %llu\n", cnt2[6]);
    fprintf(f, "uniquely mapped reads (pair) (mapQ >= %u): %llu\n", mapQ, cnt2[7]);
    fprintf(f, "MRE filtered reads:\t%llu\n", cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4]);
    fprintf(f, "    CCGG reads:\t%llu\n", cnt[0]);
    fprintf(f, "    CCGC reads:\t%llu\n", cnt[1]);
    fprintf(f, "    GCGC reads:\t%llu\n", cnt[2]);
    fprintf(f, "    ACGT reads:\t%llu\n", cnt[3]);
    fprintf(f, "    CGCG reads:\t%llu\n", cnt[4]);
    fprintf(f, "    Unkown reads:\t%llu\n", cnt[5]);
    fprintf(f, "Sampled CpG sites:\t%llu\n", cnt1);
    fprintf(f, "solo ends:\t%i\n", solosite);
    fprintf(f, "    reads on solo ends:\t%i\n", soloreads);
    fprintf(f, "fragments:\t%i\n", pairsite);
    fprintf(f, "    reads in higher end of fragments:\t%i\n", highend);
    fprintf(f, "    reads in lower end of fragments:\t%i\n", lowend);
    fprintf(f, "fragments size distribution:\n");
    fprintf(f, "%s\t%10s\t%10s\t%c\n", "Scale", "Count", "Percent", '|');
    fprintf(f, "<=%i\t%10d\t%10.2f\t|%s\n", minlen, histoMin, histoMin*100.0/pairsite, print_bar((int)(histoMin*100.0/pairsite)));
    fprintf(fout2, "dat[%d, ] <- c('<=%i', %i, %.2f)\n", 1, minlen, histoMin, histoMin*100.0/pairsite);
    fragdistro[0].label = minlen;
    fragdistro[0].pct = histoMin*100.0/pairsite;
    for (j = 0; j < bins; j++){
        fprintf(f, "%i\t%10d\t%10.2f\t|%s\n", scale[j].e, scale[j].histo, scale[j].histo*100.0/pairsite, print_bar((int)(scale[j].histo*100.0/pairsite)));
        fprintf(fout2, "dat[%d, ] <- c('%i', %i, %.2f)\n", j+2, scale[j].e, scale[j].histo, scale[j].histo*100.0/pairsite);
        fragdistro[j+1].label = scale[j].e;
        fragdistro[j+1].pct = scale[j].histo*100.0/pairsite;
    }
    fprintf(f, ">%i\t%10d\t%10.2f\t|%s\n", maxlen, histoMax, histoMax*100.0/pairsite, print_bar((int)(histoMax*100.0/pairsite)));
    fprintf(fout2, "dat[%d, ] <- c('>%i', %i, %.2f)\n", j+2, maxlen, histoMax, histoMax*100.0/pairsite);
    fragdistro[j+1].label = maxlen;
    fragdistro[j+1].pct = histoMax*100.0/pairsite;
    fragdistro[j+2].label = 999; //end mark
    
    carefulClose(&f);
    fprintf(fout2, "dat$total <- as.numeric(dat$total)\n");
    fprintf(fout2, "dat$cov <- as.numeric(dat$cov)\n");
    fprintf(fout2, "pdf('%s.fragmentsizedistro.pdf')\n", output);
    fprintf(fout2, "op <- par(mar = c(5,7,4,2) + 0.1)\n");
    fprintf(fout2, "barplot(rev(dat$cov), names=rev(dat$label), main=\"MRE fragments size ditribution\", xlab=\"Percentage\", ylab=\"\", col=5, las=1, horiz=TRUE)\n");
    fprintf(fout2, "title(ylab = \"Size range\", cex.lab = 1.5, line = 4.5)\n");
    fprintf(fout2, "par(op)\n");
    fprintf(fout2, "dev.off()\n");
    fclose(fout2);

    char *command2;
    if (asprintf(&command2, "Rscript %s", tmpRfile2) < 0)
        errAbort("Preparing command wrong");
    if (system(command2) == -1)
        fprintf(stderr, "failed to call R for plotting");
    unlink(tmpRfile2);
    
    //mre reads plot
    char tmpRfile[50];
    strcpy(tmpRfile, template_name);
    int fd = mkstemp(tmpRfile);
    if (fd == -1)
        errAbort("create temp file error.");
    FILE *fout = fdopen(fd, "w");
    fprintf(fout, "dat <- data.frame(count=c(%llu, %llu, %llu, %llu, %llu, %llu), label=c('CCGG','CCGC','GCGC', 'ACGT', 'CGCG', 'Unkown'))\n", cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], cnt[5]);
    fprintf(fout, "dat$count <- as.numeric(dat$count)\n");
    fprintf(fout, "pdf('%s.fragments.pdf')\n", output);
    fprintf(fout, "op <- par(mar = c(7,8,4,2) + 0.1)\n");
    fprintf(fout, "bar <- barplot(dat$count, las=1, col=2:7, ylab='',xlab='', main='MRE fragments stats')\n");
    fprintf(fout, "title(ylab = 'Fragments count', cex.lab = 1.5, line = 4.5)\n");
    fprintf(fout, "axis(1, at=bar, labels=dat$label, padj=1, tick=FALSE)\n");
    fprintf(fout, "par(op)\n");
    fprintf(fout, "dev.off()\n");
    fclose(fout);

    char *command;
    if (asprintf(&command, "Rscript %s", tmpRfile) < 0)
        errAbort("Preparing command wrong");
    if (system(command) == -1)
        fprintf(stderr, "failed to call R for plotting");
    unlink(tmpRfile);
    
    //mapping stats plot
    char tmpRfile1[50];
    strcpy(tmpRfile1, template_name);
    int fd1 = mkstemp(tmpRfile1);
    if (fd1 == -1)
        errAbort("create temp file error.");
    FILE *fout1 = fdopen(fd1, "w");
    fprintf(fout1, "dat <- data.frame(count=c(%llu, %llu, %llu, %llu), label=c('total\\nfragments','mapped\\nfragments','uniquely\\nmapped\\nfragments', 'MRE\\nfiltered\\nfragments'))\n", cnt2[0], cnt2[6], cnt2[7], cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4]);
    fprintf(fout1, "dat$count <- as.numeric(dat$count)\n");
    fprintf(fout1, "pdf('%s.mappingStat.pdf')\n", output);
    fprintf(fout1, "op <- par(mar = c(7,8,4,2) + 0.1)\n");
    fprintf(fout1, "bar <- barplot(dat$count, las=1, col=2:5, ylab='',xlab='', main='Mapping stats')\n");
    fprintf(fout1, "title(ylab = 'Fragments count', cex.lab = 1.5, line = 4.5)\n");
    fprintf(fout1, "axis(1, at=bar, labels=dat$label, padj=1, tick=FALSE)\n");
    fprintf(fout1, "par(op)\n");
    fprintf(fout1, "dev.off()\n");
    fclose(fout1);

    char *command1;
    if (asprintf(&command1, "Rscript %s", tmpRfile1) < 0)
        errAbort("Preparing command wrong");
    if (system(command1) == -1)
        fprintf(stderr, "failed to call R for plotting");
    unlink(tmpRfile1);

    //for (j=0; j<(bins+2); j++){
    //    fprintf(stderr, "j:%d,label:%d,pct:%.2f\n",j,fragdistro[j].label,fragdistro[j].pct);
    //}
    return fragdistro;
}

char *print_bar(int x){
    char *s;
    s = malloc(x);
    int i;
    for (i = 0; i < x; i++){
        s[i] = '*';
    }
    s[i] = '\0';
    return s;
}

void sortBedfile(char *bedfile) {
    struct lineFile *lf = NULL;
    FILE *f = NULL;
    struct bedLine *blList = NULL, *bl;
    char *line;
    int lineSize;

    lf = lineFileOpen2(bedfile, TRUE);
    while (lineFileNext(lf, &line, &lineSize)){
        if (line[0] == '#')
            continue;
        bl = bedLineNew(line);
        slAddHead(&blList, bl);
    }
    lineFileClose(&lf);

    slSort(&blList, bedLineCmp);

    f = mustOpen(bedfile, "w");
    for (bl = blList; bl != NULL; bl = bl->next){
        fprintf(f, "%s\t%s\n", bl->chrom, bl->line);
        if (ferror(f)){
    	    perror("Writing error\n");
	    errAbort("%s is truncated, sorry.", bedfile);
	}
    }
    carefulClose(&f);
    bedLineFreeList(&blList);
}

struct hash *calGenomeCovBedGraph(char *chrsize, char *bedgraph){
    struct hash *hash = hashNameIntFile(chrsize);
    struct hash *cov = newHash(0);
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(hash);
    while((hel = hashNext(&cookie)) != NULL) {
        struct gcov *g = malloc(sizeof(struct gcov));
        g->total = ptToInt(hel->val);
        g->cov = 0;
        hashAdd(cov, hel->name, g);
    }
    char *row[20], *line, prechr[100] = "start";
    int start, end, size, cnt = 0;
    struct lineFile *stream = lineFileOpen2(bedgraph, TRUE);
    while ( lineFileNextReal(stream, &line)){
        int numFields = chopByWhite(line, row, ArraySize(row));
        if (numFields < 4)
            errAbort("file %s doesn't appear to be in bedGraph format. At least 4 fields required, got %d", bedgraph, numFields);
        start = (int) strtol(row[1], NULL, 0);
        end = (int) strtol(row[2], NULL, 0);
        size = end - start;
        if ( !sameWord(row[0], prechr)){
            //a new chromosome starts
            struct hashEl *hel2 = hashLookup(cov, prechr);
            if (hel2 != NULL) {
                struct gcov *g = (struct gcov *) hel2->val;
                g->cov = cnt;
            }
            cnt = size;
        }else{
            cnt += size;
        }
        strcpy(prechr, row[0]);
    }
    //last chromosome
    struct hashEl *hel2 = hashLookup(cov, prechr);
    if (hel2 != NULL) {
        struct gcov *g = (struct gcov *) hel2->val;
        g->cov = cnt;
    }
    lineFileClose(&stream);
    return cov;
}

void genMeDIPTex(char *prefix, unsigned long long int *cnt, long long fragbase, int *covCnt, long long *countCnt, struct slInt *slPair, struct hash *chrHash, struct hash *cov, char *optm){
    char *outfile;
    if (asprintf(&outfile, "%s.tex", prefix) < 0)
        errAbort("Preparing output wrong");
    FILE *f = mustOpen(outfile, "w");
    fprintf(f, "\\documentclass[12pt]{article}\n");
    fprintf(f, "\n");
    fprintf(f, "\\usepackage{fullpage}\n");
    fprintf(f, "\\usepackage{amsfonts}\n");
    fprintf(f, "\\usepackage{graphicx}\n");
    fprintf(f, "\\usepackage{hyperref}\n");
    fprintf(f, "\n");
    fprintf(f, "\\begin{document}\n");
    fprintf(f, "\n");
    fprintf(f, "\\title{%s Report}\n", texTitleEscape(prefix));
    fprintf(f, "\\author{methylQA\\footnote{Website: \\url{http://methylQA.sourceforge.net/}} \\ version\\ %s}\n", methylQA_VERSION);
    fprintf(f, "\\date{\\today}\n");
    fprintf(f, "\\maketitle\n");
    fprintf(f, "\\pagenumbering{roman}\n");
    fprintf(f, "\\tableofcontents\n");
    fprintf(f, "\\clearpage\n");
    fprintf(f, "\\setcounter{page}{1}\n");
    fprintf(f, "\\pagenumbering{arabic}\n");
    fprintf(f, "\n");
    fprintf(f, "\\section{Mapping status}\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\begin{tabular}{|c|c|}\n");
    fprintf(f, "\\hline\n");
    fprintf(f, "total fragments & %llu \\\\ \\hline\n", cnt[0]);
    fprintf(f, "mapped fragments & %llu \\\\ \\hline\n", cnt[6]);
    fprintf(f, "uniquely mapped fragments & %llu \\\\ \\hline\n", cnt[7]);
    fprintf(f, "non-redundant uniquely mapped fragments & %llu \\\\ \\hline\n", cnt[9]);
    fprintf(f, "\\end{tabular}\n");
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\includegraphics[width=6.5in]{{%s.mappingStat}.pdf}\n", prefix);
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    if(slPair != NULL){
        fprintf(f, "\\section{Fragments size distribution}\n");
        fprintf(f, "\\begin{center}\n");
        fprintf(f, "\\includegraphics[width=6.5in]{{%s.insertdistro}.pdf}\n", prefix);
        fprintf(f, "\\end{center}\n");
        fprintf(f, "\n");
    }
    fprintf(f, "\\section{Genomic coverage}\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\begin{tabular}{|c|c|c|c|}\n");
    fprintf(f, "\\hline\n");
    fprintf(f, "chromosome & total base & covered base & covered percentage \\\\ \\hline\n");
    //fprintf(f, "chr1 & 1111111 & 111111 & 0.79 \\\\ \\hline\n");
    struct hashEl *hel;
    long long t1 = 0, c1 = 0;
    struct hashCookie cookie = hashFirst(cov);
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct gcov *g = (struct gcov *) hel->val;
        fprintf(f, "%s & %i & %i & %.2f\\%% \\\\ \\hline\n", texTitleEscape(hel->name), g->total, g->cov, (((float)(g->cov))/(g->total)) * 100.0);
        t1 += (long long)g->total;
        c1 += (long long)g->cov;
    }
    fprintf(f, "%s & %lli & %lli & %.2f\\%% \\\\ \\hline\n", "Whole Genome", t1, c1, (((double)c1)/t1) * 100.0);
    fprintf(f, "\\end{tabular}\n");
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\includegraphics[width=6.5in]{{%s.genomeCov}.pdf}\n", prefix);
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    if (optm != NULL){
        fprintf(f, "\\section{CpG status}\n");
        fprintf(f, "\\subsection{CpG coverage}\n");
        fprintf(f, "\\begin{center}\n");
        fprintf(f, "\\begin{tabular}{|c|c|}\n");
        fprintf(f, "\\hline\n");
        fprintf(f, "Range & Count \\\\ \\hline\n");
        //fprintf(f, "0 & 45455455 \\\\ \\hline\n");
        int i;
        for(i=0;i<19;i++){
            if (covCnt[i] != 0)
                fprintf(f, "%s & %i \\\\ \\hline\n", cpglabel[i], covCnt[i]);
        }
        fprintf(f, "\\end{tabular}\n");
        fprintf(f, "\\end{center}\n");
        fprintf(f, "\n");
        fprintf(f, "\\begin{center}\n");
        fprintf(f, "\\includegraphics[width=6.5in]{{%s.cpgCoverage}.pdf}\n", prefix);
        fprintf(f, "\\end{center}\n");
        fprintf(f, "\\subsection{CpG count}\n");
        fprintf(f, "\\begin{center}\n");
        fprintf(f, "\\begin{tabular}{|c|c|}\n");
        fprintf(f, "\\hline\n");
        fprintf(f, "Range & Count \\\\ \\hline\n");
        //fprintf(f, "0 & 45455455 \\\\ \\hline\n");
        for(i=0;i<19;i++){
            if (countCnt[i] != 0)
                fprintf(f, "%s & %lli \\\\ \\hline\n", cpglabel[i], countCnt[i]);
        }
        fprintf(f, "\\end{tabular}\n");
        fprintf(f, "\\end{center}\n");
        fprintf(f, "\n");
        fprintf(f, "\\begin{center}\n");
        fprintf(f, "\\includegraphics[width=6.5in]{{%s.cpgCount}.pdf}\n", prefix);
        fprintf(f, "\\end{center}\n");
        fprintf(f, "\n");
        long long cpgnum = 0;
        for(i=0;i<19;i++){
            cpgnum += (long long)covCnt[i];
        }
        long long genomebase = hashIntSum(chrHash);
        double frac = ((double)countCnt[19] / (double)fragbase) / ((double)cpgnum / (double)genomebase);
        fprintf(f, "\\section{CpG enrichment}\n");
        fprintf(f, "\\begin{eqnarray*}\n");
        fprintf(f, "CpG\\ enrichment &= \\frac{CpG\\ count\\ in\\ fragments\\Big/total\\ base\\ of\\ fragments}{CpG\\ count\\ in\\ genome\\Big/total\\ base\\ of\\ genome} \\\\ \n");
        fprintf(f, "&= \\frac{%lli\\Big/%lli}{%lli\\Big/%lli} \\\\ \n", countCnt[19], fragbase, cpgnum, genomebase);
        fprintf(f, "&= \\textbf{%.2f}\n", frac);
        fprintf(f, "\\end{eqnarray*}\n");
        fprintf(f, "\n");
    }
    fprintf(f, "\\end{document}\n");
    carefulClose(&f);
}

void genMRETex(char *prefix, unsigned long long int *cnt2, unsigned long long int *cnt, unsigned long long int cnt1, struct hash *chrHash, struct hash *cpgHash, long long *cnt3, struct fragd *fragdistro){
    char *outfile;
    if (asprintf(&outfile, "%s.tex", prefix) < 0)
        errAbort("Preparing output wrong");
    FILE *f = mustOpen(outfile, "w");
    fprintf(f, "\\documentclass[12pt]{article}\n");
    fprintf(f, "\n");
    fprintf(f, "\\usepackage{fullpage}\n");
    fprintf(f, "\\usepackage{amsfonts}\n");
    fprintf(f, "\\usepackage{graphicx}\n");
    fprintf(f, "\\usepackage{hyperref}\n");
    fprintf(f, "\n");
    fprintf(f, "\\begin{document}\n");
    fprintf(f, "\n");
    fprintf(f, "\\title{%s Report}\n", texTitleEscape(prefix));
    fprintf(f, "\\author{methylQA\\footnote{Website: \\url{http://methylQA.sourceforge.net/}} \\ version\\ %s}\n", methylQA_VERSION);
    fprintf(f, "\\date{\\today}\n");
    fprintf(f, "\\maketitle\n");
    fprintf(f, "\\pagenumbering{roman}\n");
    fprintf(f, "\\tableofcontents\n");
    fprintf(f, "\\clearpage\n");
    fprintf(f, "\\setcounter{page}{1}\n");
    fprintf(f, "\\pagenumbering{arabic}\n");
    fprintf(f, "\n");
    fprintf(f, "\\section{Mapping status}\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\begin{tabular}{|c|c|}\n");
    fprintf(f, "\\hline\n");
    fprintf(f, "total fragments & %llu \\\\ \\hline\n", cnt2[0]);
    fprintf(f, "mapped fragments & %llu \\\\ \\hline\n", cnt2[6]);
    fprintf(f, "uniquely mapped fragments & %llu \\\\ \\hline\n", cnt2[7]);
    fprintf(f, "MRE filtered fragments & %llu \\\\ \\hline\n", cnt[0]+cnt[1]+cnt[2]+cnt[3]+cnt[4]);
    fprintf(f, "\\end{tabular}\n");
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\includegraphics[width=6.5in]{{%s.mappingStat}.pdf}\n", prefix);
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    fprintf(f, "\\section{MRE fragments stats}\n");
    fprintf(f, "\\subsection{Distribution by enzyme cut site}\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\begin{tabular}{|c|c|}\n");
    fprintf(f, "\\hline\n");
    fprintf(f, " Cut site & Count \\\\ \\hline\n");
    fprintf(f, " CCGG & %llu \\\\ \\hline\n", cnt[0]);
    fprintf(f, " CCGC & %llu \\\\ \\hline\n", cnt[1]);
    fprintf(f, " GCGC & %llu \\\\ \\hline\n", cnt[2]);
    fprintf(f, " ACGT & %llu \\\\ \\hline\n", cnt[3]);
    fprintf(f, " CGCG & %llu \\\\ \\hline\n", cnt[4]);
    fprintf(f, " Unkown & %llu \\\\ \\hline\n", cnt[5]);
    fprintf(f, "\\end{tabular}\n");
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\includegraphics[width=6.5in]{{%s.fragments}.pdf}\n", prefix);
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    fprintf(f, "\\subsection{Distribution by fragment size}\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\begin{tabular}{|c|c|}\n");
    fprintf(f, "\\hline\n");
    fprintf(f, "Range & Percentage \\\\ \\hline\n");
    //fprintf(f, "0 & 45455455 \\\\ \\hline\n");
    int i, size=0;
    for (i = 0;;i++){
        if (fragdistro[i].label == 999){
            break;
        }
        size++;
    }
    char sym[20];
    for(i=0;i<size;i++){
        if (i == 0){
            strcpy(sym, "$ \\le $");
        }else if (i == (size - 1)){
            strcpy(sym, "\\textgreater");
        }else{
            strcpy(sym, " ");
        }
        fprintf(f, "%s%i & %.2f\\%% \\\\ \\hline\n", sym, fragdistro[i].label, fragdistro[i].pct);
    }
    fprintf(f, "\\end{tabular}\n");
    fprintf(f, "\\end{center}\n");
    fprintf(f, "\n");
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\includegraphics[width=6.5in]{{%s.fragmentsizedistro}.pdf}\n", prefix);
    fprintf(f, "\\end{center}\n");
    long long cpgnum = 0;
    struct hashEl *hel;
    struct hashCookie cookie = hashFirst(cpgHash);
    while ( (hel = hashNext(&cookie)) != NULL ) {
        struct binKeeper *bk = (struct binKeeper *) hel->val;
        struct binKeeperCookie becookie = binKeeperFirst(bk);
        struct binElement *be;
        while( (be = binKeeperNext(&becookie)) != NULL ){
            cpgnum += 1;
        }
        binKeeperFree(&bk);
    }
    fprintf(f, "\\section{Sampled CpG sites}\n");
    fprintf(f, "Total \\begin{Huge}%llu\\end{Huge} CpG sites (%.2f\\%%) have been sampled in this dataset.\n", cnt1, (double)cnt1*100.0/cpgnum);
    fprintf(f, "\n");
    long long genomebase = hashIntSum(chrHash);
    double frac = ((double)cnt3[1] / (double)cnt3[0]) / ((double)cpgnum / (double)genomebase);
    fprintf(f, "\\section{CpG enrichment}\n");
    fprintf(f, "\\begin{eqnarray*}\n");
    fprintf(f, "CpG\\ enrichment &= \\frac{CpG\\ count\\ in\\ fragments\\Big/total\\ base\\ of\\ fragments}{CpG\\ count\\ in\\ genome\\Big/total\\ base\\ of\\ genome} \\\\ \n");
    fprintf(f, "&= \\frac{%lli\\Big/%lli}{%lli\\Big/%lli} \\\\ \n", cnt3[1], cnt3[0], cpgnum, genomebase);
    fprintf(f, "&= \\textbf{%.2f}\n", frac);
    fprintf(f, "\\end{eqnarray*}\n");
    fprintf(f, "\n");
    fprintf(f, "\\end{document}\n");
    carefulClose(&f);
}

void tex2pdf(char *prefix){
    char *command;
    if (asprintf(&command, "pdflatex %s >/dev/null 2>&1", prefix) < 0)
        errAbort("Preparing command wrong");
    if (system(command) == -1)
        fprintf(stderr, "failed to call pdflatex for generating PDF report.");
    if (system(command) == -1) //when there is toc, twice pdflatex needed
        fprintf(stderr, "failed to call pdflatex for generating PDF report.");
    //clean stuff
    char *tf1, *tf2, *tf3, *tf4, *tf5;
    if (asprintf(&tf1, "%s.log", prefix) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&tf2, "%s.aux", prefix) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&tf3, "%s.toc", prefix) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&tf4, "%s.tex", prefix) < 0)
        errAbort("Preparing output wrong");
    if (asprintf(&tf5, "%s.out", prefix) < 0)
        errAbort("Preparing output wrong");
    unlink(tf1);
    unlink(tf2);
    unlink(tf3);
    unlink(tf4);
    unlink(tf5);
}

struct hash *chromHashFrombbiFile(char *bbiFile)
/* get chrom sizes hash from bbi file. */
{
struct hash *hash = hashNew(0);
struct bbiFile *bwf = bigWigFileOpen(bbiFile);
struct bbiChromInfo *chrom, *chromList = bbiChromList(bwf);
for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    {
        hashAddInt(hash, chrom->name, chrom->size);
    }
bbiChromInfoFreeList(&chromList);
bbiFileClose(&bwf);
return hash;
}

void bigWigToBedGraph2(char *inFile, char *outFile, float scale, int tominus)
/* bigWigToBedGraph - Convert from bigWig to bedGraph format.. */
{
int minus = 1;
if(tominus)
    minus = -1;
struct bbiFile *bwf = bigWigFileOpen(inFile);
FILE *f = mustOpen(outFile, "w");
struct bbiChromInfo *chrom, *chromList = bbiChromList(bwf);
for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    {
    boolean firstTime = TRUE;
    int saveStart = -1, prevEnd = -1;
    double saveVal = -1.0;

    char *chromName = chrom->name;
    struct lm *lm = lmInit(0);
    int start = 0, end = chrom->size;
    struct bbiInterval *interval, *intervalList = bigWigIntervalQuery(bwf, chromName, 
    	start, end, lm);
    for (interval = intervalList; interval != NULL; interval = interval->next)
	{
	if (firstTime)
	    {
	    saveStart = interval->start;
	    saveVal = (interval->val) * scale * minus;
	    firstTime = FALSE;
	    }
	else
	    {
	    if (!((prevEnd == interval->start) && (saveVal == (interval->val) * scale * minus )))
		{
		fprintf(f, "%s\t%u\t%u\t%g\n", chromName, saveStart, prevEnd, saveVal);
		saveStart = interval->start;
		saveVal = (interval->val) * scale * minus;
		}

	    }
	prevEnd = interval->end;
	}
    if (!firstTime)
	fprintf(f, "%s\t%u\t%u\t%g\n", chromName, saveStart, prevEnd, saveVal);

    lmCleanup(&lm);
    }
bbiChromInfoFreeList(&chromList);
carefulClose(&f);
bbiFileClose(&bwf);
}

void bwScale(char *bwfile, char *outbwfile, char *outbedgraph, float scale, int tominus)
{
    /*
     *naive way of re-scale a bigwig file
     to bedgraph, scale the value, then bedgraph to bigwig again
     * */
    if (!isBigWig(bwfile))
        errAbort("error, %s is not a bigWig file", bwfile);
    fprintf(stderr, "* Applying scale\n");
    struct hash *chrHash = chromHashFrombbiFile(bwfile);
    bigWigToBedGraph2(bwfile, outbedgraph, scale, tominus);
    bigWigFileCreate2(outbedgraph, chrHash, outbwfile);
}
