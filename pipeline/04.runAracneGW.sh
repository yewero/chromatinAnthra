#!/bin/bash

jar_dir="/home/tianf/lustre/Tools/ARACNe-AP/dist"

# threshold
java -Xmx5G -jar $jar_dir/Aracne.jar -e ../data/aracne_APtumor.exp  -o ../data/ARACNE_all_res --tfs ../data/aracne_APgeneList.exp --pvalue 1E-8 --seed 1 --calculateThreshold --threads 20

# bootstrap
for i in {1..100}
do
echo $i
java -Xmx150G -jar $jar_dir/Aracne.jar  -e ../data/aracne_APtumor.exp -o ../data/ARACNE_all_res --tfs ../data/aracne_APgeneList.exp --pvalue 1E-8 --threads 20 --seed $i
done
# consolidate
java -Xmx25G -jar $jar_dir/Aracne.jar -o ../data/ARACNE_all_res/ --consolidate
