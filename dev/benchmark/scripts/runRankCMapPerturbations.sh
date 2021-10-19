#!/bin/bash

threads=1
diffExpr="input/diffExprStat.rds"

for loadZscores in TRUE FALSE; do
    for type in all; do
    #for type in knockdown overexpression compound; do
        log=logs/noDebugger_rankCMap_${type}_load-${loadZscores}_${threads}-threads.log
        R -f R/rankCMapPerturbations.R --args ${diffExpr} ${type} \
                                              ${loadZscores} ${threads} > ${log} 2>&1 &
    done
done
