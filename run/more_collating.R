nom <- c("chunks1_sub.RData", "chunks101_sub.RData", "chunks201_sub.RData", "chunks301_sub.RData", "chunks_sub.RData")

allsub <- data.frame()
for(i in 1:length(nom))
{
    load(nom[i])
    allsub <- rbind(allsub, subdat)
}

save(allsub, file="allsub.RData")

