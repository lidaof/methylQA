CC=gcc
COPT= -O -g
CFLAGS= -Wall -Wformat -Wimplicit -Wreturn-type -Wuninitialized -Wno-unused-but-set-variable -Wno-unused-result
DFLAGS= -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DCOLOR32

%.o: %.c
	$(CC) $(COPT) $(CFLAGS) $(DFLAGS) -I. -o $@ -c $<


O = aliType.o annoAssembly.o annoFilter.o annoFormatter.o annoFormatTab.o \
    annoGrator.o annoGrateWig.o annoGratorQuery.o annoOption.o annoRow.o \
    annoStreamer.o annoStreamBigBed.o annoStreamBigWig.o annoStreamTab.o annoStreamVcf.o \
    apacheLog.o asParse.o axt.o axtAffine.o bamFile.o base64.o \
    basicBed.o bbiRead.o bbiWrite.o bigBed.o binRange.o bits.o \
    blastOut.o blastParse.o boxClump.o boxLump.o bPlusTree.o \
    bwgCreate.o bwgQuery.o bwgValsOnChrom.o \
    cda.o chain.o chainBlock.o chainConnect.o chainToAxt.o chainToPsl.o \
    cheapcgi.o cirTree.o codebias.o colHash.o common.o correlate.o crTree.o \
    dgRange.o diGraph.o dlist.o dnaLoad.o dnaMarkov.o dnaMotif.o dnaseq.o \
    dnautil.o dtdParse.o dyOut.o dystring.o elmTree.o \
    emblParse.o errCatch.o errabort.o \
    fa.o ffAli.o ffScore.o fieldedTable.o filePath.o fixColor.o flydna.o fof.o \
    fuzzyShow.o \
    gapCalc.o gdf.o gemfont.o genomeRangeTree.o \
    gfNet.o gff.o gff3.o gfxPoly.o gifLabel.o \
    hacTree.o hash.o hex.o histogram.o hmmPfamParse.o hmmstats.o htmlPage.o htmshell.o \
    https.o intExp.o intValTree.o internet.o itsa.o iupac.o \
    jointalign.o jpegSize.o keys.o knetUdc.o kxTok.o linefile.o lineFileOnBigBed.o localmem.o log.o \
    maf.o mafFromAxt.o mafScore.o mailViaPipe.o md5.o memalloc.o memgfx.o meta.o metaWig.o mgCircle.o \
    mgPolygon.o mime.o net.o nib.o nibTwo.o nt4.o numObscure.o \
    obscure.o oldGff.o oligoTm.o options.o osunix.o pairHmm.o paraFetch.o peakCluster.o \
    phyloTree.o pipeline.o portimpl.o pngwrite.o psGfx.o psPoly.o pscmGfx.o \
    psl.o pslGenoShow.o pslShow.o pslTbl.o pslTransMap.o pthreadWrap.o \
    qa.o quickHeap.o quotedP.o \
    ra.o rainbow.o rbTree.o rangeTree.o regexHelper.o repMask.o \
    rle.o rnautil.o rqlEval.o rqlParse.o rudp.o \
    scoreWindow.o seg.o seqOut.o seqStats.o servBrcMcw.o servCrunx.o \
    servcis.o servcl.o servmsII.o servpws.o shaRes.o slog.o snof.o \
    snofmake.o snofsig.o spaceSaver.o spacedColumn.o spacedSeed.o \
    splatAli.o sqlList.o sqlNum.o subText.o sufa.o sufx.o synQueue.o \
    tabRow.o textOut.o tokenizer.o trix.o twoBit.o \
    udc.o vcf.o vcfBits.o vGfx.o vPng.o verbose.o \
    wildcmp.o wormdna.o \
    xAli.o xa.o xap.o xenshow.o xmlEscape.o xp.o zlibFace.o


T = libcuskent.a

$(T): $(O)
	ar rcus $(T) $(O)

clean:
	rm -f $(O) $(T)
