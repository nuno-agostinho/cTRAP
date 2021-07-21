#!/bin/bash

diffExpr="input/diffExprStat.rds"

for loadZscores in TRUE FALSE; do
    for type in knockdown overexpression compound; do
        log=logs/profvis_rankCMapPerturbations_${type}_loadZscores-${loadZscores}.log
        R -f R/rankCMapPerturbations_profvis.R \
                        ${diffExpr} ${type} \
                        ${loadZscores} > ${log} 2>&1 &
    done
done
