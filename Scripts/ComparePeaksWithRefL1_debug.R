param <- ScanBamParam(what=c('pos', 'qwidth'),
                      which=GRanges(seqnames = Chrom, ranges = IRanges(start = 1, end = ChromLength)),
                                    flag=scanBamFlag(isUnmappedQuery=FALSE))
x <- scanBam(BamFile, param=param)[[1]]


NrChromPieces <- 100

33115910

i <- 2
Chrom       <- names(ChromLengths)[i]
ChromLength <- ChromLengths[i]
Ends <- seq(1, ChromLength, floor(ChromLength/ NrChromPieces))
if (Ends[length(Ends)] < ChromLength) Ends <- c(Ends, ChromLength)
cat("Extracting reads of chromosome", Chrom, "\n")
j <- 14
R1 <- GRanges(seqnames = Chrom, ranges = IRanges(start = Ends[j]+ 1.51*10^6, 
                                                 end = Ends[j] + 1.51*10^6 + 2*10^4))
Reads <- extractReads(bam.file = BamFile, region = R1)
ReadCov <- coverage(Reads)
Islands <- slice(ReadCov, lower = 1)
GRs <- GRanges(seqnames = Chrom, 
                 ranges = Islands@listData[[1]]@ranges,
                 coverTotal = viewSums(Islands)[[1]],
                 coverMax   = viewMaxs(Islands)[[1]],
                 coverMaxPos   = viewWhichMaxs(Islands)[[1]])
GRs

# Results:
# Range 33125910- 33135910: no problem
# Problem around 33135910, chromosome 2
# samtools stats -c 1,10000000,1000000 /srv/gsfs0/projects/levinson/hzudohna/10Xcapture/NA12878_capt10X_aln2hg19.bam chr2:33135910-33145910
# 2878_capt10X_aln2hg19.bam chr2:33135910-33145910