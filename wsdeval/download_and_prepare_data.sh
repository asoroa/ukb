#!/bin/bash

# Change me
# path to the UKB binaries
path_to_ukb=../bin

# Do not change anything below
abspath_to_ukb=$(realpath ${path_to_ukb})
root_dir=$(realpath .)
script_dir=${root_dir}/../scripts

###########################################
# Create WN30 WordNet graph and dictionary
###########################################
install -d wn30
cd wn30
# download and extract
wget -q http://wordnetcode.princeton.edu/3.0/WordNet-3.0.tar.gz
tar xzf WordNet-3.0.tar.gz
wget -q http://wordnetcode.princeton.edu/glosstag-files/WordNet-3.0-glosstag.tar.bz2
tar xjf WordNet-3.0-glosstag.tar.bz2 
# create ukb relation files from source (including gloss relations)
perl ${script_dir}/wnet2graph.pl WordNet-3.0/dict/* > wn30_rel.txt
perl ${script_dir}/wnetgloss2graph.pl WordNet-3.0/dict/index.sense  WordNet-3.0/glosstag/merged > wn30_gloss_rel.txt
# compile graph
cat wn30_rel.txt wn30_gloss_rel.txt | ${abspath_to_ukb}/compile_kb -o wn30g.bin --note "wn30_rel.txt wn30_gloss_rel.txt" -
# create dictionary
perl ${script_dir}/wnet2dict.pl WordNet-3.0/dict/index.sense > wn30_dict.txt
# create a mapping between wordnet database ids and lexicographic codes
perl ${root_dir}/create_wnidmap.pl WordNet-3.0/dict/index.sense > id2lc.map
cd ..

####################################
# Download WSDEval and prepare input
####################################
install -d wsdeval_src
cd wsdeval_src
# download and extract
wget -q http://lcl.uniroma1.it/wsdeval/data/WSD_Unified_Evaluation_Datasets.zip
unzip -o WSD_Unified_Evaluation_Datasets.zip >& /dev/null
# create 'raw' contexts
perl ${root_dir}/wsdeval2ukb.pl WSD_Unified_Evaluation_Datasets/ALL/ALL.data.xml > wsdeval_raw.txt
# create input file by grouping sentences until contexts have at least 20 words
perl ${root_dir}/ctx20words.pl wsdeval_raw.txt > wsdeval.txt
# compile evaluation script
cd WSD_Unified_Evaluation_Datasets
javac Scorer.java
cd ${root_dir}

