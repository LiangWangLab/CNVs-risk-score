#linux command line for counting reads in genes
#https://www.gencodegenes.org/human/
featureCounts -T 32 -t transcript -a bedtools_genes/hg19.ncbiRefSeq.gtf -o bam.counts bams/*.bam
