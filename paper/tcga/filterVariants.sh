#!/bin/bash

tcgaMAFFile=$1
outFile=$2

less $tcgaMAFFile |\
grep -Ff filterWords |\
cut -f 1,16 -d $'\t' |\
sort -u > $outFile

