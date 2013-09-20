# You have to have python, numpy and trynka scripts
# lots of requirements you say? well try to figure out how to
# run the scripts yourself first

# Moved all data over to the nemo cluster
# trynka scripts found here www.broadinstitute.org/mpg/epigwas/
# The data comes from folder DATA/final_h3k4me3_lists

# -rw-r--r--  1 kshakhbazov   1.2K 11 Jun 18:48 bonf_all.snpmappings.txt
# -rw-r--r--  1 kshakhbazov   1.0K 11 Jun 18:48 bonf_cis.snpmappings.txt
# -rw-r--r--  1 kshakhbazov   270B 11 Jun 18:48 bonf_trans.snpmappings.txt
# -rw-r--r--  1 kshakhbazov   2.1K 11 Jun 18:48 chr_int_all.snpmappings.txt
# -rw-r--r--  1 kshakhbazov   1.9K 11 Jun 18:48 chr_int_cis.snpmappings.txt
# -rw-r--r--  1 kshakhbazov   272B 11 Jun 18:48 chr_int_trans.snpmappings.txt
# -rw-r--r--  1 kshakhbazov    19K 11 Jun 18:48 levis_all.snpmappings.txt
# -rw-r--r--  1 kshakhbazov   7.8K 11 Jun 18:48 levis_cis.snpmappings.txt
# -rw-r--r--  1 kshakhbazov    11K 11 Jun 18:48 levis_trans.snpmappings.txt

# per list folder
for i in *.snpmappings.txt; do mkdir -p  ${i/.snpmappings.txt/}/ld; done
# move input file into corresponding folder
for i in *.snpmappings.txt; do mv $i  ${i/.snpmappings.txt/}; done
# Actual run; 3 steps aka 3 different scripts
# the code would print command per each list to the terminal 
# once happy pipe into shell i.e. for i in *; do echo "lots of stuff here" | sh; done
# all right lets count 1
for i in *_*; do echo "nohup python /hox/u/uqkshakb/data/gosia/scripts/computeLd/computeLd.py /hox/u/uqkshakb/data/gosia/1kg/beagle/ phase1_integrated_calls.20101123.ALL.panel EUR $i/$i.snpmappings.txt 500000 0.8 $i/ld $i/$i &" ; done
# 2
for i in *_*; do echo "nohup python /hox/u/uqkshakb/data/gosia/scripts/prepFiles_v01/makeFiles.py $i/$i.1kg.map.txt $i/ld/ /hox/u/uqkshakb/data/gosia/H3K4me3/ /hox/u/uqkshakb/data/gosia/H3K4me3/listTissues.txt H3K4me3 0.8 $i/$i &"; done
# 3
for i in *_*; do echo "nohup python /hox/u/uqkshakb/data/gosia/scripts/phenoCellSpec_v01/phenoCellSpecif.py $i/$i.H3K4me3.snpPeak.txt $i/$i.ld.txt $i/$i.H3K4me3.tisIndex.txt /hox/u/uqkshakb/data/gosia/backgroundSNPs/bg.H3K4me3.snpPeak.txt /hox/u/uqkshakb/data/gosia/backgroundSNPs/bg.ld.txt 10000 $i/$i.H3K4me3.10k &"; done
# Results were copied back from the cluster into folder DATA/H3K4me3_output_final