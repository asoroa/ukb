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
ctx=../input/ctx.txt
ctx_dos=../input/ctx_dos.txt
graphSrc=../input/test_graph.txt
../../compile_kb -o $gbin ${graphSrc}
../../ukb_wsd --all --ppr_w2w -D ${dict} -K $gbin ${ctx} > $dir/wsd_w2w.txt
../../ukb_wsd --all --ppr -D ${dict} -K $gbin ${ctx} > $dir/wsd_ppr.txt
../../ukb_wsd --all --ppr -D ${dict} -K $gbin ${ctx_dos} > $dir/wsd_dos_ppr_ctxweight.txt
../../ukb_wsd --all --ppr_w2w -D ${dict} -K $gbin ${ctx_dos} > $dir/wsd_dos_pprw2w_ctxweight.txt
../../ukb_wsd --dict_weight --all --ppr -D ${dict} -K $gbin ${ctx} > $dir/wsd_ppr_dictweight.txt
../../ukb_wsd --all --static -D ${dict} -K $gbin ${ctx} > $dir/wsd_static.txt
#../ukb_wsd --all --ppr_sim -D ${dict} -K $gbin ${ctx} > $dir/wsd_sim.txt
#../ukb_wsd --all --ppr_qd -D ${dict} -K $gbin ${ctx} > $dir/wsd_qd.txt
