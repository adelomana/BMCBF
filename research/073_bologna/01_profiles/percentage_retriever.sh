#!/usr/bin/env bash
cd messages || exit 1
for f in *.err.txt; do
    awk '/\[quant\] processed/{p=$3; a=$5; gsub(/,/,"",p); gsub(/,/,"",a)}
         END{if(p) printf "%s\t%.2f%%\n", FILENAME, 100*a/p}' "$f"
done