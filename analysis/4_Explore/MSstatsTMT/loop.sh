#!/usr/bin/env bash

# desperate times call for desperate measures
for i in {1..8586}
do
	echo "iter: $i"
	./fit1.R $i
done
