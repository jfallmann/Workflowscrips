#!/usr/bin/env bash

inf=$1

if [[ -s $inf ]]
then
	if [[ "$inf" == *.gz* ]]
	then
		zcat $inf |samtools faidx -
	else
		samtools faidx $inf
	fi
else
	touch $inf.fai
fi
