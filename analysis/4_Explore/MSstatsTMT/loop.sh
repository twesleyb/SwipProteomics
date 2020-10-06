#!/usr/bin/env bash

# desperate times call for desperate measures
# loop through all proteins, fit with lmer and 
# assess differences between conditioned means
# append to results.csv

for i in {1..8586}
do
	echo "iter: $i"
	./fit1.R $i
done
