library(tidyverse)
library(parallel)

fn <- list.files("../data/scratch", pattern=".rdata$", recursive=T)

res <- mclapply(fn, function(x) {
	load(file.path("../data/scratch", x))
	as_tibble(res)
	}, mc.cores=10) %>% bind_rows()

save(res, file="../data/aggregate.rdata")




