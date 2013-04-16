dir list on 16/5/12

-rwxr-xr-x   1 ghemani  DI\Domain Users       3478 23 Apr 16:40 Make_plink_data_final.R
-rw-r--r--   1 ghemani  DI\Domain Users   13648551 15 May 11:13 eqtl2d_objects.RData
-rwxr-xr-x   1 ghemani  DI\Domain Users        231 23 Apr 16:40 extract_geno_plink.sh
-rw-r--r--   1 ghemani  DI\Domain Users  120863770 15 May 11:39 geno.RData
-rw-r--r--   1 ghemani  DI\Domain Users       1379 15 May 11:11 organise_data.R
drwxr-xr-x   7 ghemani  DI\Domain Users        238 15 May 11:59 raw_data
drwxr-xr-x  15 ghemani  DI\Domain Users        510 15 May 11:58 unused

Make_plink_data_final.R by Joseph to create binary ped files and phenotype files and covariate files (these are in 'raw_data')
extract_geno_plink.sh converts plink to 012 format
organise_data.R manipulates raw data to be used by epiGPU and R on the cluster in ../v4/
- this creates the file eqtl2d_objects.RData



