#!/bin/bash

## This short script moves things to "trash"
while read -r line;
do
	echo "$line";
	trash_dir=$(dirname ${line})
	mkdir -p trash/${trash_dir}
	mv $line trash/$line
done < trash.txt
