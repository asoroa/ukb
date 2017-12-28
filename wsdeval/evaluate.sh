#!/bin/bash

kdir=$1
if [ $# -ne 1 ]; then
    echo "Evaluating in Keys"
	kdir="Keys"
fi

function split_kfile {
	fname=$1
	oDir=$2
	method=$3
	install -d $oDir
	cut -d '.' -f 1 $fname | sort | uniq | while read -r line; do
		dataset=$(echo $line | cut -d " " -f 1)
		ofname="$oDir/$dataset.$method.key"
		grep $dataset $fname | sed -e "s/^[^\.]\+\.//" > $ofname
	done
}


for i in ${kdir}/raw/*.key; do
    method=`basename $i .key`
    perl create_keyfile.pl wn30/id2lc.map $i > ${kdir}/ALL.${method}.key
    split_kfile ${kdir}/ALL.${method}.key ${kdir} ${method}
	echo "* $method"
	for dataset in $(echo "ALL semeval2007  semeval2013  semeval2015  senseval2 senseval3"); do
		sys=${kdir}/${dataset}.${method}.key
		gs=wsdeval_src/WSD_Unified_Evaluation_Datasets/${dataset}/${dataset}.gold.key.txt
        r=$(java -cp wsdeval_src/WSD_Unified_Evaluation_Datasets Scorer $gs ${sys})
		rr=$(echo $r | sed -e 's/\n/ /g')
		echo "	$dataset $rr"
    done
    echo
done
