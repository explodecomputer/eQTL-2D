load("../analysis/data_lists.RData")
replicated_ids <- unique(c(meta_replicated$snp1, meta_replicated$snp2))
write.table(replicated_ids, file = "DATA/replicated_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#---------------------------------------------------------------
# bsgs_imputed_R2_80_cleaned_stage2_chr_all files come from Joseph
#---------------------------------------------------------------
 plink --bfile bsgs_imputed_R2_80_cleaned_stage2_chr_all \
       --r2 \
       --ld-snp-list replicated_ids.txt \
       --ld-window-kb 1000 \
       --ld-window-r2 0.2 \
       --out replicated_ld_1mb.txt