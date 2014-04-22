
# Analayis of regression model function
# Joseph Powell
# 18/04/2014


# load TMEM149 datafile

# Object names TMEM149, contains snp genotypes and probe levels for 846 individual
load("/Users/jpowell/repo/eQTL-2D/Frayling_et_al/For_Jian/TMEM149.RData")

# 4 and 8df tests	

probe <- bsgs_TMEM149$ILMN_1786426 
snp1 <- bsgs_TMEM149$rs8106959
snp2 <- bsgs_TMEM149$rs6926382
inc <- bsgs_TMEM149$IncSeq
adj_probe <- bsgs_TMEM149$adj_ILMN_1786426

# LD
cor(snp1,snp2)^2
cor(snp1,inc)^2
cor(snp2,inc)^2

# boxplots
boxplot(probe~as.numeric(as.matrix(snp1)))
boxplot(probe~as.numeric(as.matrix(snp2)))
boxplot(probe~as.numeric(as.matrix(inc)))

summary(lm(inc~snp1+snp2+snp1:snp2))


# note, the direction of the effect by the inc snp is in the opposit direction

# Original analysis
fullmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2) + as.factor(snp1):as.factor(snp2))
redmod <- lm(probe ~ as.factor(snp1) + as.factor(snp2))
intmod <- anova(redmod, fullmod)	

# Adjusting the model using the genotypes of SNP 1
new_probe <- summary(lm(probe~snp1))$residuals
new_fullmod <- lm(new_probe ~ as.factor(snp1) + as.factor(snp2) + as.factor(snp1):as.factor(snp2))
new_redmod <- lm(new_probe ~ as.factor(snp1) + as.factor(snp2))
newintmod <- anova(new_redmod, new_fullmod)	


# Using the adjusted phenotype
adj_fullmod <- lm(adj_probe ~ as.factor(snp1) + as.factor(snp2) + as.factor(snp1):as.factor(snp2))
adj_redmod <- lm(adj_probe ~ as.factor(snp1) + as.factor(snp2))
adj_intmod <- anova(adj_redmod, adj_fullmod)	


