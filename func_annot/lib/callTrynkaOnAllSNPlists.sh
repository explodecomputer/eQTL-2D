for i in *.snpmappings.txt; do mkdir -p  ${i/.snpmappings.txt/}/ld; done
for i in *.snpmappings.txt; do mv $i  ${i/.snpmappings.txt/}; done

for i in *_*; do echo "nohup python /hox/u/uqkshakb/data/gosia/scripts/computeLd/computeLd.py /hox/u/uqkshakb/data/gosia/1kg/beagle/ phase1_integrated_calls.20101123.ALL.panel EUR $i/$i.snpmappings.txt 500000 0.8 $i/ld $i/$i &" ; done
"nohup python /hox/u/uqkshakb/data/gosia/scripts/computeLd/computeLd.py /hox/u/uqkshakb/data/gosia/1kg/beagle/ phase1_integrated_calls.20101123.ALL.panel EUR $i/$i.snpmappings.txt 500000 0.8 $i/ld $i/$i &"


for i in *_*; do echo "nohup python /hox/u/uqkshakb/data/gosia/scripts/prepFiles_v01/makeFiles.py $i/$i.1kg.map.txt $i/ld/ /hox/u/uqkshakb/data/gosia/H3K4me3/ /hox/u/uqkshakb/data/gosia/H3K4me3/listTissues.txt H3K4me3 0.8 $i/$i &"; done


for i in *_*; do echo "nohup python /hox/u/uqkshakb/data/gosia/scripts/phenoCellSpec_v01/phenoCellSpecif.py $i/$i.H3K4me3.snpPeak.txt $i/$i.ld.txt $i/$i.H3K4me3.tisIndex.txt /hox/u/uqkshakb/data/gosia/backgroundSNPs/bg.H3K4me3.snpPeak.txt /hox/u/uqkshakb/data/gosia/backgroundSNPs/bg.ld.txt 10000 $i/$i.H3K4me3.10k &"; done
