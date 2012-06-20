#!/bin/bash

for scr in scripts/*.sh; do
	echo $scr
	j=$(basename ${scr})
	(cd scripts; bash ./${j})
done