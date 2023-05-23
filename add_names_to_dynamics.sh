#!/bin/bash
#Usage: bash add_names_to_dynamics.sh <interaction dynamics file from R script>
paste <(cut -f 2 $1 | while read line; do grep -w -m1 $line names.txt || echo -e "$line\t$line"; done | cut -f 2) <(cut -f 3 $1 | while read line; do grep -w -m1 $line names.txt || echo -e "$line\t$line"; done | cut -f 2) <(cut -f 6 $1) | tail -n +2 | sort | uniq > named_$1

