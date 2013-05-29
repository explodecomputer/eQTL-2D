sig$type <- "cis-cis"
sig$type[with(sig, chr1 == probechr & chr2 != probechr)] <- "cis-trans"
sig$type[with(sig, chr1 != probechr & chr2 == probechr)] <- "cis-trans"
sig$type[with(sig, chr1 != probechr & chr2 != probechr)] <- "trans-trans"

for(i in 1:nrow(sig))
{
	sig$vc_fehr[[i]] <- sig$vc_fehr[[i]][[1]]
	sig$vc_egcut[[i]] <- sig$vc_egcut[[i]][[1]]
}

bsgs <- as.data.frame(do.call(rbind, sig$vc))
fehr <- as.data.frame(do.call(rbind, sig$vc_fehr))
egcut <- as.data.frame(do.call(rbind, sig$vc_egcut))

bsgs$varA <- bsgs$a. + bsgs$.a
fehr$varA <- fehr$a. + fehr$.a
egcut$varA <- egcut$a. + egcut$.a

bsgs$varD <- bsgs$d. + bsgs$.d
fehr$varD <- fehr$d. + fehr$.d
egcut$varD <- egcut$d. + egcut$.d

bsgs$varI <- bsgs$aa + bsgs$ad + bsgs$da + bsgs$dd
fehr$varI <- fehr$aa + fehr$ad + fehr$da + fehr$dd
egcut$varI <- egcut$aa + egcut$ad + egcut$da + egcut$dd


varA <- data.frame(bsgs = bsgs$varA, fehr = fehr$varA, egcut = egcut$varA)
varD <- data.frame(bsgs = bsgs$varD, fehr = fehr$varD, egcut = egcut$varD)
varI <- data.frame(bsgs = bsgs$varI, fehr = fehr$varI, egcut = egcut$varI)

pairs(varA)
pairs(varD)
pairs(varI)
