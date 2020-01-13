library(tidyverse)

fn <- list.files("../data/scratch", pattern=".rdata$", recursive=T)

res <- lapply(fn, function(x) {
	load(x)
	as_tibble(res)
	}) %>% bind_rows()

save(res, file="../data/aggregate.rdata")




