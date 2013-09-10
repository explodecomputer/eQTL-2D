# Downloaded Gencode16 annotation in ugly gtf format and the annotation eats lots of RAM
# converting to sqlite db via GenomicFeatures
#-------------------------------------------------------------------------------
library(GenomicFeatures)

hg19_chr_info <- getChromInfoFromUCSC("hg19")

gencode16TxDb <- makeTranscriptDbFromGFF("gencode.v16.annotation.gtf.gz",
                                                          format = "gtf",
                                                         species = "Homo sapiens",
                                                      dataSource = "ftp://ftp.sanger.ac.uk/pub/gencode/release_16/gencode.v16.annotation.gtf.gz",
                                                       chrominfo = hg19_chr_info)

saveDb(gencode16TxDb, file = "Gencode16.sqlite")
#-------------------------------------------------------------------------------