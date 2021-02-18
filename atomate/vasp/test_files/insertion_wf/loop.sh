#!/bin/bash
cur=$(pwd)
for f in *static*
do
    cd "$f/outputs"
    rm AECCAR?.gz
    gzip AECCAR0
    gzip AECCAR2
    cd $cur
done
