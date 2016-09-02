
if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` result1_dir result2_dir"
  exit 65 # E_BADARGS
fi

dirA=$(echo $1 | sed -e "s/\/$//")
dirB=$(echo $2 | sed -e "s/\/$//")

function ddif {
	fileA=$1
	bname=$(basename ${fileA})
	fileB=$(echo $fileA | sed -e "s|^${dirA}|${dirB}|")
	tmpA=$(mktemp)
	egrep -v "^!!" $1 > ${tmpA}
	tmpB=$(mktemp)
	egrep -v "^!!" $fileB > ${tmpB}
	dout=$(diff -pu ${tmpA} ${tmpB} | wc -l)
	if [ $dout -gt 0 ] ; then
		echo ${bname} " " ${dout}
	fi
	rm ${tmpA}
	rm ${tmpB}
}

for i in $(find ${dirA} -name "*.txt" -o -name "*.ppv" ); do
	ddif $i
done
