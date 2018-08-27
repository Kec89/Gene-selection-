#Tea Kecman

# This is to select ORFs with no overlapping transcriptional unit (on the same strand) 275 nt downstream the PAS
#Protein-coding genes that do not overlap with any other transcription unit (on the opposite strand) in a region from 250 nt upstream to 500 nt downstream of the PAS were selected for the analysis
#of 4105 genes from Eser et al (2016), 1735 ORFs are selected

library(rtracklayer)
#import annotation file
allTUs <- import("/Users/Tea/Desktop/NGSanalysis/Annotation/Pombe/Eser_GFFformat.gff")#input is Eser et al (2016) annotation

minus_allTUs <- subset(allTUs, strand == "-")
plus_allTUs <- subset(allTUs, strand == "+")

m.after_gene <- width(gaps(minus_allTUs))
p.before_gene <- width(gaps(plus_allTUs))

p.after_gene <- as.vector(c(p.before_gene[2:length(p.before_gene)], 20000))

p.old.meta<-as.data.frame(values(plus_allTUs))
p.new.meta<-cbind(p.old.meta, p.after_gene)
names(p.new.meta) <- c(names(p.old.meta), "after_gene")
values(plus_allTUs)<-p.new.meta

m.old.meta<-as.data.frame(values(minus_allTUs))
m.new.meta<-cbind(m.old.meta, m.after_gene)
names(m.new.meta) <- c(names(m.old.meta), "after_gene")
values(minus_allTUs)<-m.new.meta

all_TUs_after <- append(plus_allTUs, minus_allTUs)


# select ORFs with no overlapping transcriptional unit (on the same strand) 275 nt downstream the PAS
after_275 <- subset(all_TUs_after, after_gene > 275)

minus_ORFs <- subset(after_275, strand == "-" & type == "ORF-TU")
plus_ORFs <- subset(after_275, strand == "+" & type == "ORF-TU")

#transcription units were excluded when their PAS overlapped with a region from 250 nt upstream to 500 nt downstream of overlapping genes' PAS.
plus_250TTS <- GRanges(seqnames = seqnames(plus_ORFs), ranges =IRanges(end(plus_ORFs)-249, end(plus_ORFs)+499), strand = "-", mcols =as.data.frame(values(plus_ORFs)))
minus_250TTS <- GRanges(seqnames = seqnames(minus_ORFs), ranges =IRanges(start(minus_ORFs)-499, start(minus_ORFs)+249), strand = "+", mcols =as.data.frame(values(minus_ORFs)))

plus_TTS <- GRanges(seqnames = seqnames(plus_allTUs), ranges =IRanges(end(plus_allTUs), end(plus_allTUs)+1), strand = "+", mcols =as.data.frame(values(plus_allTUs)))
minus_TTS <- GRanges(seqnames = seqnames(minus_allTUs), ranges =IRanges(start(minus_allTUs), start(minus_allTUs)+1), strand = "-", mcols =as.data.frame(values(minus_allTUs)))

plus_no_anno <- plus_250TTS[!(plus_250TTS %over% minus_TTS)]
minus_no_anno <- minus_250TTS[!(minus_250TTS %over% plus_TTS)]

plus_down <- subset(plus_ORFs, mcols(plus_ORFs)$group %in% mcols(plus_no_anno)$mcols.group)
minus_down <- subset(minus_ORFs, mcols(minus_ORFs)$group %in% mcols(minus_no_anno)$mcols.group)



plus_500pA <- GRanges(seqnames = seqnames(plus_down), ranges =IRanges(end(plus_down)-500, end(plus_down)+500), strand = "+", mcols =as.data.frame(values(plus_down)))
minus_500pA <- GRanges(seqnames = seqnames(minus_down), ranges =IRanges(start(minus_down)-500, start(minus_down)+500), strand = "-", mcols =as.data.frame(values(minus_down)))

plus.meta<-as.data.frame(values(plus_500pA))
minus.meta<-as.data.frame(values(minus_500pA))
tot <- rbind(plus.meta, minus.meta)

#output
write.csv(tot, file="/Users/User/Desktop/NGSanalysis/Annotation/list_genes.csv")