#!/usr/bin/env bash

inf=$1

if [[ -s $inf ]]
then
	samtools faidx $inf
else
	touch $inf.fai
fi
