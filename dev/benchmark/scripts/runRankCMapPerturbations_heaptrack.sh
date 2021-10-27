#!/bin/bash

# Before running this script, open /usr/bin/R and comment lines 291-297 and 299

threads=1
diffExpr="input/diffExprStat.rds"

for loadZscores in TRUE FALSE; do
    #for type in knockdown overexpression compound; do
    for type in all; do
        log=logs/heaptrack_rankCMap_${type}_load-${loadZscores}_${threads}-threads.log
        R -d heaptrack -f R/rankCMapPerturbations.R \
                              --args ${diffExpr} ${type} \
                                     ${loadZscores} ${threads} > ${log} 2>&1 &
    done
done
