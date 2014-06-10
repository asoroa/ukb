#!/bin/bash

if [ $# -gt 0 ] ; then
    ver=$1
else
    ver=$(../../compile_kb --version)
fi

echo $ver
rootdir=../results/v${ver}
dir=${rootdir}/main_wsd
install -d $dir
gbin=$rootdir/graph.bin
dict=../input/dict.txt
ctx=../input/ctx_nodict.txt
ctx_dos=../input/ctx_dos.txt
graphSrc=../input/test_graph.txt
../../compile_kb -o $gbin ${graphSrc}
../../ukb_wsd --all --nopos --ppr_w2w -D ${dict} -K $gbin ${ctx} > $dir/wsd_w2w_nopos_dictweight_nodict.txt
../../ukb_wsd --all --ppr -D ${dict} -K $gbin ${ctx} > $dir/wsd_ppr_nopos_dictweight_nodict.txt
