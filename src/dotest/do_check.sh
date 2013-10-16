
if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` result1_dir result2_dir"
  exit 65 # E_BADARGS
fi

dirA=$1
dirB=$2

function ddif {
	fileA=$1
	fileB=$(echo $fileA | sed -e "s|^${dirA}|${dirB}|")
	tmpA=$(mktemp)
	egrep -v "^!!" $1 > ${tmpA}
	tmpB=$(mktemp)
	egrep -v "^!!" $fileB > ${tmpB}
	#echo -n "${fileA} ${fileB} "
	echo -n $(basename ${fileA}) " "
	diff -pu ${tmpA} ${tmpB} | wc -l
	rm ${tmpA}
	rm ${tmpB}
}

for i in $(find ${dirA} -name "*.txt" -o -name "*.ppv" ); do
	ddif $i
done
