#!/usr/bin/env bash

# Summary stats
(echo -e "sample\tcategory\tcount" && for file in `find -type f -name "*stats.txt"`; do sample=`echo $file | cut -f2 -d"/"`; awk -v sample=$sample '($0 ~ /^#/){split($0,arr,":"); print sample"\t"substr(arr[1],2,100)"\t"arr[2]}' $file; done;) > stats_summary.txt

# Summary counts:
(head -n1 `ls */*counts.txt | head -n1` | awk '{print "sample\t"$0}' && ls */*counts.txt | xargs -I file awk '{split(FILENAME,arr,"/"); fn=arr[1]; if(NR>1){print fn"\t"$0}}' file) > counts_summary.txt
