source("http://bioconductor.org/biocLite.R")
biocLite("AnnotationDbi")
biocLite("hgu95av2.db")
library(AnnotationDbi)
library(hgu95av2.db)
library(org.Hs.eg.db)

symbolKeys <- c("BRCA1")
cols <- c("ACCNUM","REFSEQ")
brca1 <- select(org.Hs.eg.db, keys=symbolKeys, cols=cols,keytype="SYMBOL")

biocLite("biomaRt")
library(biomaRt)
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene.id <- brca1[1,]$REFSEQ
gene <- getGene(id = gene.id, type = "refseq_mrna", mart = mart)

#seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="cdna", mart = mart)
seq = getSequence(chromosome=gene$chromosome_name,
                  start=gene$start_position, end=gene$end_position , 
                  type="refseq_mrna",
                  seqType="coding", mart = mart)
show(seq)