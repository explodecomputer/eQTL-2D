library(lattice)
library(latticeExtra)
?panel.3dbars

tab <- as.table(matrix(-4:4, 3, 3))
cloud(tab, panel.3d.cloud=panel.3dbars, col="black", col.facet=2:4, xbase=0.9, ybase=0.9)



