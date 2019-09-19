library(tidyverse)

d <- "../data/scratch/rs67903230"
fn <- list.files(d) %>% grep(".rdata", ., value=TRUE) %>% paste(d, ., sep="/")
fn

res <- lapply(fn, function(x) {
	load(x)
	as_tibble(res)
	}) %>% bind_rows()

save(res, file="../data/rs67903230.rdata")




