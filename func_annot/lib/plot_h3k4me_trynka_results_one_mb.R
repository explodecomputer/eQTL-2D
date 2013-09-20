# Plot H3K4me3 enrichment results
# blame k.shakhbazov@uq.edu.au

# the data comes in as 3 files(per SNP list): tissues scores, snp scores and 95% for snp scores 

# Lots of boring manipulations to read in data bind all files into a data frame and 
# order all labels as we would like to, plus lots of boring plotting routines
# to make ggplot looks beautiful 

bindMe <- function(tis_score) {
  
  br <- c("Brain_Cingulate_Gyrus",
          "Brain_Anterior_Caudate",
          "Brain_Substantia_Nigra",
          "Brain_Inferior_Temporal_Lobe",
          "Brain_Mid_Frontal_Lobe",
          "Brain_Hippocampus_Middle"
  )

  hemato <- c("CD34_Primary_Cells",
              "Mobilized_CD34_Primary_Cells",
              "CD3_Primary_Cells",
              "CD19_Primary_Cells",
              "CD8_Memory_Primary_Cells",
              "CD8_Naive_Primary_Cells",
              "CD34_Cultured_Cells",
              "CD4_Naive_Primary_Cells",
              "CD4_Memory_Primary_Cells",
              "Treg_Primary_Cells",
              "Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells"
  )
            
  gastro <- c("Colonic_Mucosa",
              "Adult_Liver",
              "Duodenum_Smooth_Muscle",
              "Duodenum_Mucosa",
              "Stomach_Smooth_Muscle",
              "Stomach_Mucosa",
              "Rectal_Smooth_Muscle",
              "Rectal_Mucosa",
              "Colon_Smooth_Muscle"
  )

  other <- c("Pancreatic_Islets",
             "Adipose_Nuclei",
             "Adult_Kidney",
             "Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
             "Muscle_Satellite_Cultured_Cells",
             "Skeletal_Muscle",
             "Chondrocytes_from_Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
             "Mesenchymal_Stem_Cell_Derived_Adipocyte_Cultured_Cells"
  )
   
  tiss_score <- read.table(str_c("DATA/H3K4me3_output_one_mb/", tis_score),  sep = "\t", header = TRUE) 
  dt <- data.frame(tissue = c(br, hemato, gastro, other))
  dt$body_part <- rep(c("brain", "hemato", "gastro", "other"), sapply(list(br, hemato, gastro, other), length))
  merge(tiss_score, dt, by.x = "Tissue", by.y = "tissue", all.x = TRUE)
}

files <- dir(path = "DATA/H3K4me3_output_one_mb/", pattern = "levis", recursive = TRUE)
names(files) <- str_split_fixed(files, "\\.", n = 2)[ ,1]

all_tiss_scores <- lapply(files, bindMe)

all_tiss_scores <- lapply(names(all_tiss_scores), function(x) {
  dt <- all_tiss_scores[[x]]
  dt$samp_ind <- x
  dt
})

all_tiss_scores <- do.call(rbind, all_tiss_scores)

organ_order <- c("brain", "other", "gastro", "hemato")

samp_order <- c("levis_cis", "levis_trans")

tissue_order <-                                   c("Brain_Anterior_Caudate",
                                                     "Brain_Cingulate_Gyrus",
                                                  "Brain_Hippocampus_Middle",
                                              "Brain_Inferior_Temporal_Lobe",
                                                    "Brain_Mid_Frontal_Lobe",
                                                    "Brain_Substantia_Nigra",
                      "Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
                                                            "Adipose_Nuclei",
                                                              "Adult_Kidney",
"Chondrocytes_from_Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
                    "Mesenchymal_Stem_Cell_Derived_Adipocyte_Cultured_Cells",
                                           "Muscle_Satellite_Cultured_Cells",
                                                         "Pancreatic_Islets",
                                                           "Skeletal_Muscle",
                                                               "Adult_Liver",
                                                       "Colon_Smooth_Muscle",
                                                            "Colonic_Mucosa",
                                                           "Duodenum_Mucosa",
                                                    "Duodenum_Smooth_Muscle",
                                                             "Rectal_Mucosa",
                                                      "Rectal_Smooth_Muscle",
                                                            "Stomach_Mucosa",
                                                     "Stomach_Smooth_Muscle",
                  "Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells",
                                                        "CD19_Primary_Cells",
                                                         "CD3_Primary_Cells",
                                                       "CD34_Cultured_Cells",
                                                        "CD34_Primary_Cells",
                                                  "CD4_Memory_Primary_Cells",
                                                   "CD4_Naive_Primary_Cells",
                                                  "CD8_Memory_Primary_Cells",
                                                   "CD8_Naive_Primary_Cells",
                                              "Mobilized_CD34_Primary_Cells",
                                                        "Treg_Primary_Cells")

all_tiss_scores <- within(all_tiss_scores, Tissue <- factor(Tissue, levels = tissue_order))
all_tiss_scores <- within(all_tiss_scores, body_part <- factor(body_part, levels = organ_order))
all_tiss_scores <- within(all_tiss_scores, samp_ind <- factor(samp_ind, levels = samp_order))

levels(all_tiss_scores$body_part) <- c("Brain", "Musculoskeletal,\nendocrine & others", "Gastrointestinal", "Hematopoietic")
levels(all_tiss_scores$samp_ind) <- c("cis SNPs (303)", "trans SNPs (475)")

qplot(x = Tissue, y = -log10(Score + 0.0001), data = all_tiss_scores, geom = "bar", stat = "identity", position = "dodge", fill = body_part) + facet_wrap(~ samp_ind) + coord_flip() + theme_bw() +
scale_x_discrete(labels = c(
Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells = "Mesenchymal stem cells (adipose)",
Adipose_Nuclei = "Adipose nuclei",
Adult_Kidney = "Adult kidney",
Adult_Liver = "Adult liver",
Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells = "Mesenchymal stem cells (bone marrow)",
Brain_Anterior_Caudate = "Anterior caudate",
Brain_Cingulate_Gyrus = "Cingulate gyrus",
Brain_Hippocampus_Middle = "Hippocampus middle",
Brain_Inferior_Temporal_Lobe = "Inferior temporal lobe",
Brain_Mid_Frontal_Lobe = "Mid-frontal lobe",
Brain_Substantia_Nigra = "Substantia nigra",
CD19_Primary_Cells = substitute(paste(CD19^"+", " primary cells")),
CD3_Primary_Cells = substitute(paste(CD3^"+", " primary cells")),
CD34_Cultured_Cells = substitute(paste(CD34^"+", " cultured cells")),
CD34_Primary_Cells = substitute(paste(CD34^"+", " primary cells")),
CD4_Memory_Primary_Cells = substitute(paste(CD4^"+", " memory primary cells")),
CD4_Naive_Primary_Cells = substitute(paste(CD4^"+", " naive primary cells")),
CD8_Memory_Primary_Cells = substitute(paste(CD8^"+", " memory primary cells")),
CD8_Naive_Primary_Cells = substitute(paste(CD8^"+", " naive primary cells")),
Chondrocytes_from_Bone_Marrow_Derived_Mesenchymal_Stem_Cell_Cultured_Cells = "Chondrocytes (mesenchymal stem cells)",
Colon_Smooth_Muscle = "Smooth muscle, colon",
Colonic_Mucosa = "Mucosa, colon",
Duodenum_Mucosa = "Mucosa, duodenum",
Duodenum_Smooth_Muscle = "Duodenum smooth muscle",
Mesenchymal_Stem_Cell_Derived_Adipocyte_Cultured_Cells = "Adipocyte, (mesenchymal stem cells)",
Mobilized_CD34_Primary_Cells = substitute(paste("Mobilized ", CD34^"+", " primary cells")),
Muscle_Satellite_Cultured_Cells = "Muscle satellite cultured cells",
Pancreatic_Islets = "Pancreatic islets",
Rectal_Mucosa = "Mucosa, rectum",
Rectal_Smooth_Muscle = "Rectal smooth muscle",
Skeletal_Muscle = "Skeletal muscle",
Stomach_Mucosa = "Mucosa, stomach",
Stomach_Smooth_Muscle = "Stomach smooth muscle",
Treg_Primary_Cells = substitute(paste(T[reg], " primary cells"))
)) + 
scale_fill_brewer(palette = "Dark2", guide = guide_legend(reverse=TRUE)) + 
theme(
      # Legend
      legend.position="right",
      legend.title=element_blank(),
      legend.background=element_blank(),
      legend.key=element_blank(),
      # Text in general
      text=element_text(family="Helvetica", size=14), 
      # Strip aka facet lables
      strip.background=element_blank(),
      # inside of plot
      panel.grid=element_line(size=0),
      panel.border=element_rect(size=.6),
      panel.grid.minor.x=element_blank(),
      # Give me back my x axis
      axis.line.x=element_line(colour="black", size=2)) + 
ylab("Enrichment for H3K4me3 peaks\n(P value)") +
scale_y_continuous(breaks=c(0, 1, 2, 3, 4), labels=expression(1, 10^-1, 10^-2, 10^-3, 10^-4)) +
geom_hline(y=-log10(0.05/(2*34)), color="darkred", linetype="longdash", size = .3)
ggsave("figures/Trynka_one_mb_rule.pdf")