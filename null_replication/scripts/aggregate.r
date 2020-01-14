library(tidyverse)
library(parallel)

args <- commandArgs(T)
suffix <- args[1]
output <- args[2]

fn <- list.files("../data/scratch", pattern=paste0(suffix, "$"), recursive=T)

res <- mclapply(fn, function(x) {
	load(file.path("../data/scratch", x))
	as_tibble(res)
	}, mc.cores=10) %>% bind_rows()

save(res, file=output)
