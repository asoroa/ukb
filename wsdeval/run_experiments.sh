#!/bin/bash

allexp=$1
rootdir=$(realpath .)
graph_fname=
dict_fname=/

# Change me
# path to the UKB binaries
path_to_ukb=../bin

# Do not change anything below
#
# the files below are created by the 'download_and_prepare_data.sh' script
input_fname=wsdeval_src/wsdeval.txt
graph_fname=wn30/wn30g.bin
dict_fname=wn30/wn30_dict.txt

# Directory for storing keyfiles
keydir=${rootdir}/Keys/raw
install -d ${keydir}

function run {
	ukb_args="$1"
	output_fname=${keydir}/$2
	echo ""
    echo -n "${path_to_ukb}/ukb_wsd ${ukb_args} -D ${dict_fname} -K ${graph_fname} ${input_fname} > ${output_fname}"
    time ${path_to_ukb}/ukb_wsd ${ukb_args} -D ${dict_fname} -K ${graph_fname} ${input_fname} > ${output_fname}
}

run "--ppr_w2w" "pprw2w.key"
# comment line below for additional experiments
if [ $# -gt 0 ]; then
	run "--ppr" "ppr.key"
	run "--dgraph_dfs --dgraph_rank ppr" "dfs.key"
	run "--ppr_w2w --nodict_weight" "pprw2w_nf.key"
	run "--ppr --nodict_weight" "ppr_nf.key"
	run "--dgraph_dfs --dgraph_rank ppr --nodict_weight" "dfs_nf.key"
	run "--static" "static.key"
fi
