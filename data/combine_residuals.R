load("residuals2.RData")
hsq2 <- hsq
probeinfo2 <- probeinfo
resphen2 <- resphen

load("residuals.RData")
hsq <- c(hsq, hsq2)
probeinfo <- rbind(probeinfo, probeinfo2)
resphen <- cbind(resphen, resphen2)


length(hsq)
dim(probeinfo)
dim(resphen)

save(hsq, probeinfo, resphen, file="residuals_all.RData")

