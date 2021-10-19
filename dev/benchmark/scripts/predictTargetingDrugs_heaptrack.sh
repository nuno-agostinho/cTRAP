#!/bin/bash

# Before running this script, open /usr/bin/R and comment lines 291-297 and 299

threads=1
diffExpr="input/diffExprStat.rds"

log=logs/heaptrack_predictTargetingDrugs.log
R -d heaptrack -f R/predictTargetingDrugs.R --args ${diffExpr} > ${log} 2>&1
