#include "sam.h"
#include "from_kent.h"

#define methylQA_VERSION "0.2.1 (r050)"
/*! @abstract supplementary (since bwa mem)*/
#define BAM_SUPP        2048


struct mreFrag {
    char pair[100], chr[50];
    unsigned long long int reads_count;
    int head, start, end;
    char site[5]; //dont forget \0
};

struct range {
    int s, e, histo;
};

struct cpgScore {
    int start;
    double score;
    char chr[50];
};

struct cpgC {
    int c; //coverage count
    int mc; //methylated c count
    int umc; //unmethylated c count
    char strand;
};

struct gcov {
    int total;
    int cov;
};

struct fragd {
    int label;
    float pct;
};

int density_usage();
int main_density(int argc, char *argv[]);
int medip_usage();
int main_medip(int argc, char *argv[]);
int cpg_usage();
int main_cpg(int argc, char *argv[]);
int atac_usage();
int main_atac(int argc, char *argv[]);
int genomecov_usage();
int main_genomecov(int argc, char *argv[]);
int bismark_usage();
int main_bismark(int argc, char *argv[]);
int bwscale_usage();
int main_bwscale(int argc, char *argv[]);

char *strrev(char *str);
char *get_filename_without_ext(char *filename);
char *get_filename_ext(char *filename);
char * texTitleEscape(char *title);
bool is_file(const char* path); 
bool is_dir(const char* path);
struct lineFile *lineFileOpen2(char *fileName, bool zTerm);
void plotMappingStat(unsigned long long int *cnt, char *prefix);
void sortBedfile(char *bedfile);
void writeReportDensity(char *outfile, unsigned long long int *cnt, unsigned int mapQ);
void writeReportBismark(char *outfile, unsigned long long int *cnt, unsigned long long int *cnt2, int numFields, char *row[100], int bisMode, long long genomeBase);
void writecpgBismarkLite(struct hash *cpgHash, char *outfilefor, char *outfilerev, int covThres);
void writecpgBismarkLiteHash(struct hash *cpgHash, char *outfilefor, char *outfilerev, int covThres);
void assignCpGcount(struct hash *chrHash, struct hash *cpgHash, struct hash *chgHash, struct hash *chhHash, char *chrom, int start, char *methycall, char strand, int left, int right, unsigned long long int *methyCnt, int fullMode);
unsigned long long int *bismarkBamParse(char *samfile, struct hash *chrHash, struct hash *cpgHash, struct hash *chgHash, struct hash *chhHash, char *forwardread, char *reverseread, int isSam, int addChr, int fullMode, unsigned int iSize);
unsigned long long int *sam2bed(char *samfile, char *outbed, struct hash *chrHash, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat);
struct hash *MREfrag2Hash (char *fragfile, int minlen, int maxlen);
unsigned long long int *filterReadByMREsite(struct hash *hash, char *inBed, char *outBed, int call, char *prefix, int guess);
double calCpGscore (struct mreFrag *mre, unsigned long long int *cnt);
unsigned long long int  CpGscorebedGraph(struct hash *hash, unsigned long long int *cnt, char *outfile);
struct fragd *fragmentStats(struct hash *hash, unsigned long long int *cnt2, unsigned int mapQ, unsigned long long int *cnt, unsigned long long int cnt1, char *outfile, int minlen, int maxlen, int win);
char *print_bar(int x);
boolean binKeeperAnyInclude(struct binKeeper *bk, int start, int end);
int binKeeperCpGstat(struct binKeeper *bk, int start, int end);
void writecpgCount(struct slInt *cpgCount, char *outfile);
long long * plotcpgCount(struct slInt *Count, char *prefix);
void writecpgCov(struct hash *cpgHash, char *outfile);
int * plotcpgCov(struct hash *cpgHash, char *prefix);
long long writeInsertsize(struct slInt *slPair, char *outfile);
long long plotInsertsize(struct slInt *slPair, char *prefix);
struct hash *cpgBed2BinKeeperHash (struct hash *chrHash, char *cpgbedfile);
struct hash *cpgBed2BinKeeperHashBismark (struct hash *chrHash, char *cpgbedfile);
long long * bedCpGstat(struct hash *cpgHash, char *bedfile);
unsigned long long int *sam2bedwithCpGstat(char *samfile, char *outbed, struct hash *chrHash, struct hash *cpgHash, struct slInt **cpgCount, struct slInt **slPair, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat);
struct hash *initGenomeCovHash(struct hash *chrHash);
void writeGenomeCov(struct hash *cov, char *outfile);
unsigned long long int *writecpgBismark(struct hash *cpgHash, char *outfile, char *outcpg, int statsOnly, int covThres);
void plotGenomeCov(struct hash *cov, char *prefix);
struct hash *calGenomeCovBedGraph(char *chrsize, char *bedgraph);
void genMeDIPTex(char *prefix, unsigned int mapQ, unsigned long long int *cnt, long long fragbase, int *covCnt, long long *countCnt, struct slInt *slPair, struct hash *chrHash, struct hash *cov, char *optm);
void genMRETex(char *prefix, unsigned long long int *cnt2, unsigned long long int *cnt, unsigned long long int cnt1, struct hash *chrHash, struct hash *cpgHash, long long *cnt3, struct fragd *fragdistro);
void tex2pdf(char *prefix);
struct hash *chromHashFrombbiFile(char *bbiFile);
void bigWigToBedGraph2(char *inFile, char *outFile, float scale, int tominus);
void bwScale(char *bwfile, char *outbwfile, char *outbedgraph, float scale, int tominus);
unsigned long long int *ATACsam2bed(char *samfile, char *outbed, struct hash *chrHash, struct slInt **slPair, int isSam, unsigned int mapQ, int rmDup, int addChr, int discardWrongEnd, unsigned int iSize, unsigned int extension, int treat, int shift);
