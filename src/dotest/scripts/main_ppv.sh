#!/bin/bash

if [ $# -gt 0 ] ; then
    ver=$1
else
    ver=$(../../compile_kb --version)
fi

echo $ver
rootdir=../results/v${ver}
dir=${rootdir}/main_ppv
install -d $dir
gbin=$rootdir/graph.bin
dict=../input/dict.txt
ctx=../input/ctx.txt
graphSrc=../input/test_graph.txt
../../compile_kb -o $gbin ${graphSrc}
../../ukb_ppv --nodict_weight --variants --static -D ${dict} -K $gbin > ${dir}/static.ppv
../../ukb_ppv --nodict_weight --variants --prefix nopos_C_ --nopos -C -O $dir -D ${dict} -K $gbin ${ctx}
../../ukb_ppv --nodict_weight --variants --prefix nopos_G_ --nopos -G -O $dir -D ${dict} -K $gbin ${ctx}
../../ukb_ppv --nodict_weight --variants --prefix pos_ -O $dir -D ${dict} -K $gbin ${ctx}
../../ukb_ppv --nodict_weight --variants --prefix pos_C_ -C -O $dir -D ${dict} -K $gbin ${ctx}
../../ukb_ppv --nodict_weight --variants --prefix pos_G_ -G -O $dir -D ${dict} -K $gbin ${ctx}

rm -f ${dir}/*_sorted.ppv >& /dev/null
for ppv_fi in ${dir}/*.ppv; do
	j=$(basename $ppv_fi \\\.ppv)
	sort -k 2 -g -r ${ppv_fi} > ${dir}/${j}_sorted.ppv
done
