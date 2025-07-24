#!/bin/bash
if [[ -z "$2" ]]
then
echo "joinbedgraph file.rep1.bedgraph file.rep2.bedgraph [output.bedgraph]"
exit
fi
if [[ -z "$3" ]]
then
    if [[ $1 =~ (.+)\.[^.]+\.([^.]+)$ ]]
    then
        new="${BASH_REMATCH[1]}.txt"
    fi
else
    new=$3
fi
TMPBG=table.bg
bedtools intersect -a $1 -b $2 -wa -wb | awk 'BEGIN { OFS="\t" } {print $1, $2, $3, $4} {print $5, $6, $7, $8}' | sort -k1,1 -k2,2n -k3,3n -o $TMPBG
bedops --partition $TMPBG | bedmap --echo --echo-map-id-uniq --delim '\t' - $TMPBG | awk '{ n = split($4, a, ";"); sum = 0; for(i = 1; i <= n; i++) { sum += a[i]; } meann = sum/n; print $1"\t"$2"\t"$3"\t"meann; }' - > $new
rm -f $TMPBG
