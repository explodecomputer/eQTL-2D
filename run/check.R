for(i in 1:1959)
{
	fn <- paste("result_2_", i, ".txt.gz", sep="")
	if(!file.exists(fn)) cat(i, "\n")
}

