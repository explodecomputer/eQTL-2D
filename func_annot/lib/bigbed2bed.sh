# Convert bigBed to Bed cause rtracklayer
# does not read bb properly at least for me
# get binaries for the Kent utility if no local UCSC mirror present
# wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed

bigBedToBed k562.combined.bb k562.combined.bed