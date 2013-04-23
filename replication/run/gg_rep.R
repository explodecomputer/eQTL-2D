library(snpStats)

#' Convert geno matrix to raw for SnpMatrix
#'
#' @param geno Matrix of 0/1/2 numeric
#' @param fam data.frame of plink fam format
#' @param bim data.frame of plink bim format
#' @export
#' @return matrix of raw values
#' @alias
#' @examples \dontrun{
#'
#'}
genoToSnpRaw <- function(geno, fam, bim)
{
    geno <- geno + 1
    a <- matrix(as.raw(geno), nrow(geno), ncol(geno))
    rownames(a) <- fam$iid
    colnames(a) <- bim$snpname
    return(a)
}


#' Write objects to plink format
#'
#' @param geno Matrix of 0/1/2 numeric
#' @param bim data.frame of plink bim format
#' @param fam data.frame of plink fam format
#' @param outfile Path to root filename
#' @export
#' @return NULL
writePlinkFiles <- function(geno, bim, fam, outfile)
{
    snpmatrix <- new("SnpMatrix", .Data=genoToSnpRaw(geno, fam, bim))
    write.plink(
	file.base        = outfile,
	snps             = snpmatrix,
	pedigree         = fam$fid,
	id               = fam$iid,
	father           = fam$father,
	mother           = fam$mother,
	sex              = fam$sex,
	phenotype        = fam$phen,
	chromosome       = bim$chr,
	genetic.distance = bim$gd,
	position         = bim$pd,
	allele.1         = bim$a1,
	allele.2         = bim$a2
    )
}   


load("~/repo/eQTL-2D/data/ggdata.RData")
    gen2 <- gen
    gen2[gen2 == "NC"] <- NA
    gen2[gen2 == "AA"] <- "0"
    gen2[gen2 == "AB"] <- "1"
    gen2[gen2 == "BB"] <- "2"
gen2 <- matrix(as.numeric(gen2), nrow(gen2), ncol(gen2))
gen2 <- t(gen2)

bim <- with(snp, data.frame(chr = Chr, snpname = Name, gd = 0, pd = Position, a1 = "A", a2 = "B"))
fam <- with(id, data.frame(fid = 1:nrow(id), iid = as.character(CHDWB_ID), father = 0, mother = 0, sex = GENDER, phen = -9))
fam$sex <- as.character(fam$sex)
fam$sex[fam$sex == "FEM"] <- 2
fam$sex[fam$sex == "MAL"] <- 1

writePlinkFiles(gen2, bim, fam, "ggdata")

phen <- t(probe)
colnames(phen) <- prinfo$PROBE_ID
rownames(phen) <- NULL
phen <- cbind(fam[,1:2], phen)
colnames(phen)[1:2] <- c("FID", "IID")

write.table(phen, file="ggdata.phen", row=F, col=T, qu=F)


####################################

R --no-save --args ./ggdata ./ggdata.phen interactions_list.RData < replication.R



