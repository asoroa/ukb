#!/bin/bash

if [ $# -gt 0 ] ; then
    ver=$1
else
    ver=$(../../compile_kb --version)
fi

echo $ver
rootdir=../results/v${ver}
dir=${rootdir}/dfs
install -d $dir
gbin=$rootdir/graph.bin
dict=../input/dict.txt
ctx=../input/ctx.txt
graphSrc=../input/test_graph.txt
../../compile_kb -o $gbin ${graphSrc}
../../ukb_wsd --all --nodict_weight --dgraph_dfs --dgraph_rank static -D ${dict} -K $gbin ${ctx} > $dir/dfs_static.txt
../../ukb_wsd --all --nodict_weight --dgraph_dfs --dgraph_rank ppr -D ${dict} -K $gbin ${ctx} > $dir/dfs_ppr.txt
../../ukb_wsd --all --nodict_weight --dgraph_dfs --dgraph_rank ppr_w2w -D ${dict} -K $gbin ${ctx} > $dir/dfs_pprw2w.txt
../../ukb_wsd --all --nodict_weight --dgraph_dfs --dgraph_rank degree -D ${dict} -K $gbin ${ctx} > $dir/dfs_degree.txt

