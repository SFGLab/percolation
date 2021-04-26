#the name of the loop file in ./loop folder without .bed extension
loop='CTCFplus2' #'CTCF.convergent' 'CTCF.tandem' 'GM12878.RNAPII' 'testloop' 'CTCFplus2.bed' 
echo $loop
mkdir ./waointersect/$loop
mkdir ./output/$loop
mkdir ./waointersect_merged/$loop

#the list of the features with peaks (bed file format)
ARRAY=('A' 'B' 'A1' 'A2' 'B1' 'B2' 'B3' 'B4' 'D001' 'D001.methylated' 'D002' 'D003' 'H001' 'H002' 'H003' 'H006' 'H009' 'H010' 'H011' 'H014' 'H015' 'H016' 'H017' 'O001' 'O002' 'O003' 'O004' 'O005' 'O006' 'S001.active' 'S001.repressed' 'S002.heterochromatin' 'S002.active' 'S002.repressed' 'S003.heterochromatin' 'S003.active' 'S003.repressed' 'S004.heterochromatin' 'S004.active' 'S004.repressed' 'R001' 'R002' 'R003' 'R004' 'R005' 'R006' 'R007' 'R008' 'R010' 'R011' 'R012' 'N001' 'N002' 'T001' 'T002' 'T003' 'T004' 'T005' 'T006' 'T007' 'T008' 'T009' 'T010' 'T011' 'T012' 'T013' 'T014' 'T015' 'T016' 'T017' 'T018' 'T019' 'T020' 'T021' 'T022' 'T023' 'T024' 'T025' 'T026' 'T027' 'T028' 'T029' 'T030' 'T031' 'T032' 'T033' 'T034' 'T035' 'T036' 'T037' 'T038' 'T039' 'T040' 'T041' 'T042' 'T043' 'T044' 'T045' 'T046' 'T047' 'T048' 'T049' 'T050' 'T051' 'T052' 'T053' 'T054' 'T055' 'T056' 'T057' 'T058' 'T059' 'T060' 'T061' 'T062' 'T063' 'T064' 'T065' 'T066' 'T067' 'T068' 'T069' 'T070' 'T071' 'T072' 'T073' 'T074' 'T075' 'T076' 'T077' 'T078' 'T079' 'T080' 'T081' 'T082' 'T083' 'T084' 'T085' 'T086' 'T087' 'T088' 'T089' 'T090' 'T091' 'T092' 'T093' 'T094' 'T095' 'T096' 'T097' 'T098' 'T099' 'T100' 'T101' 'T102' 'P001')

ELEMENTS=${#ARRAY[@]}
for (( i=0;i<$ELEMENTS;i++)); do
    date
    echo ${ARRAY[${i}]}
    #calculate intersection between loops and feature peaks
    bedtools intersect -wao -a ./loop/$loop.bed -b ./4columns/${ARRAY[${i}]}.bed > ./waointersect/$loop/${ARRAY[${i}]}.$loop.waointersect.txt
    #calculates mean sum min max from input in bed format
    python pipe.py ./waointersect/$loop/${ARRAY[${i}]}.$loop.waointersect.txt ${ARRAY[${i}]} ./output/$loop/${ARRAY[${i}]}.$loop.output.txt
    #merge overlaps between peaks to calculate fraction corectlly
    bedtools merge -i ./4columns/${ARRAY[${i}]}.bed > ./merged/${ARRAY[${i}]}.bed
    bedtools intersect -wao -a ./loop/$loop.bed -b ./merged/${ARRAY[${i}]}.bed > ./waointersect_merged/$loop/${ARRAY[${i}]}.$loop.waointersect.txt
    #calculates fraction from input in merged peaks bed format
    python pipe.py ./waointersect_merged/$loop/${ARRAY[${i}]}.$loop.waointersect.txt ${ARRAY[${i}]}.merged ./output/$loop/${ARRAY[${i}]}.merged.$loop.output.txt
done

#the list of the features in bigwig format (signal for whole genome not only peaks)
ARRAY=('H001.bw' 'H002.bw' 'H003.bw' 'H006.bw' 'H009.bw' 'H010.bw' 'H011.bw' 'H014.bw' 'H015.bw' 'H016.bw' 'H017.bw' 'O002s1.bw' 'O002s2.bw' 'O006s1.bw' 'O006s2.bw' 'GCPe')
ARRAY2=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
#ARRAY=('GCPe')
#ARRAY2=('chr19')
ELEMENTS=${#ARRAY[@]}
ELEMENTS2=${#ARRAY2[@]}
for (( i=0;i<$ELEMENTS;i++)); do
    echo -e 'chr\tstart\tend\t'${ARRAY[${i}]}'_mean\t'${ARRAY[${i}]}'_sum\t'${ARRAY[${i}]}'_fraction\t'${ARRAY[${i}]}'_min\t'${ARRAY[${i}]}'_max\t'${ARRAY[${i}]}'_mean.pbp' > ./output/$loop/${ARRAY[${i}]}.$loop.output.txt 
    for (( j=0;j<$ELEMENTS2;j++)); do
        date
        echo ${ARRAY[${i}]} ${ARRAY2[${j}]} 
        #writes mean sum fraction from input in bigwig format (to avoid memory error for each chromosome separately)
        ./bigWigToBedGraph -chrom=${ARRAY2[${j}]} ./bigwig/${ARRAY[${i}]}.bigwig ./bigwig/bedgraph/${ARRAY[${i}]}.${ARRAY2[${j}]}.bedgraph
        grep -P "${ARRAY2[${j}]}\t" ./loop/$loop.bed > ./loop/tmp_$loop.bed 
        bedtools intersect -wao -a ./loop/tmp_$loop.bed -b ./bigwig/bedgraph/${ARRAY[${i}]}.${ARRAY2[${j}]}.bedgraph > ./waointersect/$loop/${ARRAY[${i}]}.$loop.waointersect.${ARRAY2[${j}]}.txt
        rm ./loop/tmp_$loop.bed
        rm ./bigwig/bedgraph/${ARRAY[${i}]}.${ARRAY2[${j}]}.bedgraph
        python pipe.py ./waointersect/$loop/${ARRAY[${i}]}.$loop.waointersect.${ARRAY2[${j}]}.txt ${ARRAY[${i}]} ./output/$loop/${ARRAY[${i}]}.$loop.${ARRAY2[${j}]}.output.txt
        rm ./waointersect/$loop/${ARRAY[${i}]}.$loop.waointersect.${ARRAY2[${j}]}.txt
        #merge the output from each chromosome to one file for a chromatin feature
        tail +2 ./output/$loop/${ARRAY[${i}]}.$loop.${ARRAY2[${j}]}.output.txt | cat >> ./output/$loop/${ARRAY[${i}]}.$loop.output.txt
        rm ./output/$loop/${ARRAY[${i}]}.$loop.${ARRAY2[${j}]}.output.txt
    done
done
     
#Selection of only those feature values that is reasonable considering the experiments from where thay origin
paste ./output/$loop/*.$loop.output.txt | tr -s '\t' ',' | csvcut -c chr,start,end,A.merged_fraction,B.merged_fraction,A1.merged_fraction,A2.merged_fraction,B1.merged_fraction,B2.merged_fraction,B3.merged_fraction,B4.merged_fraction,GCPe_mean.pbp,O003_sum,O003.merged_fraction,O005_sum,O005.merged_fraction,N001_sum,N001.merged_fraction,D001_mean,D001_sum,D002_mean,D002_sum,D003_mean,D003_sum,P001_mean,D001.methylated.merged_fraction,S001.active.merged_fraction,S001.repressed.merged_fraction,S002.heterochromatin.merged_fraction,S002.active.merged_fraction,S002.repressed.merged_fraction,S003.heterochromatin.merged_fraction,S003.active.merged_fraction,S003.repressed.merged_fraction,S004.heterochromatin.merged_fraction,S004.active.merged_fraction,S004.repressed.merged_fraction,N002_mean,N002_sum,N002.merged_fraction,O001_mean,O001_sum,O001.merged_fraction,O002_mean,O002_sum,O002.merged_fraction,O004_mean,O004_sum,O004.merged_fraction,O006_mean,O006_sum,O006.merged_fraction,R001_mean,R001_sum,R001.merged_fraction,R002_mean,R002_sum,R002.merged_fraction,R003_mean,R003_sum,R003.merged_fraction,R004_mean,R004_sum,R004.merged_fraction,R005_mean,R005_sum,R005.merged_fraction,R006_mean,R006_sum,R006.merged_fraction,R007_mean,R007_sum,R007.merged_fraction,R008_mean,R008_sum,R008.merged_fraction,R010_mean,R010_sum,R010.merged_fraction,R011_mean,R011_sum,R011.merged_fraction,R012_mean,R012_sum,R012.merged_fraction,T001_mean,T001_sum,T001.merged_fraction,T002_mean,T002_sum,T002.merged_fraction,T003_mean,T003_sum,T003.merged_fraction,T004_mean,T004_sum,T004.merged_fraction,T005_mean,T005_sum,T005.merged_fraction,T006_mean,T006_sum,T006.merged_fraction,T007_mean,T007_sum,T007.merged_fraction,T008_mean,T008_sum,T008.merged_fraction,T009_mean,T009_sum,T009.merged_fraction,T010_mean,T010_sum,T010.merged_fraction,T011_mean,T011_sum,T011.merged_fraction,T012_mean,T012_sum,T012.merged_fraction,T013_mean,T013_sum,T013.merged_fraction,T014_mean,T014_sum,T014.merged_fraction,T015_mean,T015_sum,T015.merged_fraction,T016_mean,T016_sum,T016.merged_fraction,T017_mean,T017_sum,T017.merged_fraction,T018_mean,T018_sum,T018.merged_fraction,T019_mean,T019_sum,T019.merged_fraction,T020_mean,T020_sum,T020.merged_fraction,T021_mean,T021_sum,T021.merged_fraction,T022_mean,T022_sum,T022.merged_fraction,T023_mean,T023_sum,T023.merged_fraction,T024_mean,T024_sum,T024.merged_fraction,T025_mean,T025_sum,T025.merged_fraction,T026_mean,T026_sum,T026.merged_fraction,T027_mean,T027_sum,T027.merged_fraction,T028_mean,T028_sum,T028.merged_fraction,T029_mean,T029_sum,T029.merged_fraction,T030_mean,T030_sum,T030.merged_fraction,T031_mean,T031_sum,T031.merged_fraction,T032_mean,T032_sum,T032.merged_fraction,T033_mean,T033_sum,T033.merged_fraction,T034_mean,T034_sum,T034.merged_fraction,T035_mean,T035_sum,T035.merged_fraction,T036_mean,T036_sum,T036.merged_fraction,T037_mean,T037_sum,T037.merged_fraction,T038_mean,T038_sum,T038.merged_fraction,T039_mean,T039_sum,T039.merged_fraction,T040_mean,T040_sum,T040.merged_fraction,T041_mean,T041_sum,T041.merged_fraction,T042_mean,T042_sum,T042.merged_fraction,T043_mean,T043_sum,T043.merged_fraction,T044_mean,T044_sum,T044.merged_fraction,T045_mean,T045_sum,T045.merged_fraction,T046_mean,T046_sum,T046.merged_fraction,T047_mean,T047_sum,T047.merged_fraction,T048_mean,T048_sum,T048.merged_fraction,T049_mean,T049_sum,T049.merged_fraction,T050_mean,T050_sum,T050.merged_fraction,T051_mean,T051_sum,T051.merged_fraction,T052_mean,T052_sum,T052.merged_fraction,T053_mean,T053_sum,T053.merged_fraction,T054_mean,T054_sum,T054.merged_fraction,T055_mean,T055_sum,T055.merged_fraction,T056_mean,T056_sum,T056.merged_fraction,T057_mean,T057_sum,T057.merged_fraction,T058_mean,T058_sum,T058.merged_fraction,T059_mean,T059_sum,T059.merged_fraction,T060_mean,T060_sum,T060.merged_fraction,T061_mean,T061_sum,T061.merged_fraction,T062_mean,T062_sum,T062.merged_fraction,T063_mean,T063_sum,T063.merged_fraction,T064_mean,T064_sum,T064.merged_fraction,T065_mean,T065_sum,T065.merged_fraction,T066_mean,T066_sum,T066.merged_fraction,T067_mean,T067_sum,T067.merged_fraction,T068_mean,T068_sum,T068.merged_fraction,T069_mean,T069_sum,T069.merged_fraction,T070_mean,T070_sum,T070.merged_fraction,T071_mean,T071_sum,T071.merged_fraction,T072_mean,T072_sum,T072.merged_fraction,T073_mean,T073_sum,T073.merged_fraction,T074_mean,T074_sum,T074.merged_fraction,T075_mean,T075_sum,T075.merged_fraction,T076_mean,T076_sum,T076.merged_fraction,T077_mean,T077_sum,T077.merged_fraction,T078_mean,T078_sum,T078.merged_fraction,T079_mean,T079_sum,T079.merged_fraction,T080_mean,T080_sum,T080.merged_fraction,T081_mean,T081_sum,T081.merged_fraction,T082_mean,T082_sum,T082.merged_fraction,T083_mean,T083_sum,T083.merged_fraction,T084_mean,T084_sum,T084.merged_fraction,T085_mean,T085_sum,T085.merged_fraction,T086_mean,T086_sum,T086.merged_fraction,T087_mean,T087_sum,T087.merged_fraction,T088_mean,T088_sum,T088.merged_fraction,T089_mean,T089_sum,T089.merged_fraction,T090_mean,T090_sum,T090.merged_fraction,T091_mean,T091_sum,T091.merged_fraction,T092_mean,T092_sum,T092.merged_fraction,T093_mean,T093_sum,T093.merged_fraction,T094_mean,T094_sum,T094.merged_fraction,T095_mean,T095_sum,T095.merged_fraction,T096_mean,T096_sum,T096.merged_fraction,T097_mean,T097_sum,T097.merged_fraction,T098_mean,T098_sum,T098.merged_fraction,T099_mean,T099_sum,T099.merged_fraction,T100_mean,T100_sum,T100.merged_fraction,T101_mean,T101_sum,T101.merged_fraction,T102_mean,T102_sum,T102.merged_fraction,H001_mean,H001_sum,H001.merged_fraction,H001.bw_mean.pbp,H002_mean,H002_sum,H002.merged_fraction,H002.bw_mean.pbp,H003_mean,H003_sum,H003.merged_fraction,H003.bw_mean.pbp,H006_mean,H006_sum,H006.merged_fraction,H006.bw_mean.pbp,H009_mean,H009_sum,H009.merged_fraction,H009.bw_mean.pbp,H010_mean,H010_sum,H010.merged_fraction,H010.bw_mean.pbp,H011_mean,H011_sum,H011.merged_fraction,H011.bw_mean.pbp,H014_mean,H014_sum,H014.merged_fraction,H014.bw_mean.pbp,H015_mean,H015_sum,H015.merged_fraction,H015.bw_mean.pbp,H016_mean,H016_sum,H016.merged_fraction,H016.bw_mean.pbp,H017_mean,H017_sum,H017.merged_fraction,H017.bw_mean.pbp,O002s1.bw_mean.pbp,O002s2.bw_mean.pbp,O006s1.bw_mean.pbp,O006s2.bw_mean.pbp | tr ',' '\t' > ./output/$loop.featuretable.txt

#Feature values considered for each experiment---------------------------------
#O003, O005, N001 (peak sum, fraction)
#D001, D002, D003 (peak mean, sum)
#P001 (peak mean)
#D001.methylated S001.active, S001.repressed, S002-S004.heterochromatin, S002-S004.active, S002-S004.repressed (peak fraction)
#T001-T102, O001, O002, O004, O006, R001-R008, R010-R012, N002 (peak mean, sum, fraction)
#H001-H003, H006, H009-H011, H014-H017 (peak mean, sum, fraction, bigwigmean)
#GCPe, O002s1, O002s2, O006s1, O006s2 (bigwigmean)




