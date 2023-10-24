#!/bin/bash

rm tasks.txt

echo "Remember to configure the desired file paths in task_dispatcher.py"

python3 task_dispatcher.py

TMPFILE=$(mktemp /tmp/grids_output-XXXXX)

echo "Now starting parallel jobs"

cat tasks.txt | parallel -j4 ./process_grids >> $TMPFILE

rm $TMPFILE
