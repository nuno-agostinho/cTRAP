#!/bin/sh

# Convert heaptrack files to massif
# Details memory allocated at specific times

for i in heaptrack*gz; do
    echo $i
    heaptrack_print --massif-detailed-freq 0 -M ${i}.massif ${i} &
done
