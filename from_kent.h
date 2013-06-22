#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "binRange.h"
#include "basicBed.h"
#include "bigWig.h"
#include "obscure.h"
#include "sqlNum.h"
#include "localmem.h"
#include "dystring.h"
#include "cirTree.h"
#include "sig.h"
#include "zlibFace.h"
#include "bPlusTree.h"
#include "bbiFile.h"
#include "bwgInternal.h"
#include "bigWig.h"

/* define unitSize to be a larger storage class if your counts
 * are overflowing. */
typedef unsigned int unitSize;

#define MAXCOUNT (unitSize)~0
#define MAXMESSAGE "Overflow of overlap counts. Max is %lu.  Recompile with bigger unitSize or use -max option"
#define INCWOVERFLOW(countArray,x) if(countArray[x] == MAXCOUNT) {if(!doMax) errAbort(MAXMESSAGE,(unsigned long)MAXCOUNT);} else countArray[x]++

struct sectionItem
/* An item in a section of a bedGraph. */
    {
    bits32 start, end;			/* Position in chromosome, half open. */
    float val;				/* Single precision value. */
    };

void writeSections(struct bbiChromUsage *usageList, struct lineFile *lf, 
	int itemsPerSlot, struct bbiBoundsArray *bounds, int sectionCount, FILE *f,
	int resTryCount, int resScales[], int resSizes[], 
	boolean doCompress, bits32 *retMaxSectionSize);
struct bbiSummary *bedGraphWriteReducedOnceReturnReducedTwice(struct bbiChromUsage *usageList, 
	int fieldCount, struct lineFile *lf, bits32 initialReduction, bits32 initialReductionCount, 
	int zoomIncrement, int blockSize, int itemsPerSlot, boolean doCompress,
	struct lm *lm, FILE *f, bits64 *retDataStart, bits64 *retIndexStart,
	struct bbiSummaryElement *totalSum);
void bedGraphToBigWig(char *inName, char *chromSizes, char *outName);

void outputCounts(unitSize *counts, char *chrom, unsigned size, FILE *f);
void bedItemOverlapCount(struct hash *chromHash, char *infile, char *outfile);
