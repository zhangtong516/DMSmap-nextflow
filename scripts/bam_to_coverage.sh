#!/usr/bin/env bash

GENOMESIZE=$3
IS5P=$4 ## either "T" or "F"
NUMHIT=$5
MINCOV=$6
IS5P=${IS5P:="T"}
BEDGRAPH="F"
if [ "$IS5P" = "T" ]; then
        FIVEP="-5"
else
        FIVEP=""
fi

if [ "$BEDGRAPH" = "F" ]; then
        DEPTHFLAG="-dz"
else
        DEPTHFLAG="-bg"
fi

samtools view -h $1 | mawk -v numhit=$NUMHIT '{if($0 ~ /^@/){print;}else{split($1,a,"XN:");if(a[2]>0 && a[2]<=numhit){print}};}' | samtools view  -S -b - | bedtools genomecov -strand + $DEPTHFLAG $FIVEP -ibam - -g  $GENOMESIZE  | mawk -v id=$2 -v mincov=$MINCOV 'BEGIN{FS="\t";OFS="\t";}{$1=$1":Pos";$2=$2+1; if($3>=mincov){print $0,id}}'

samtools view -h $1 | mawk -v numhit=$NUMHIT '{if($0 ~ /^@/){print;}else{split($1,a,"XN:");if(a[2]>0 && a[2]<=numhit){print}};}' | samtools view -S -b - | bedtools genomecov -strand - $DEPTHFLAG $FIVEP -ibam - -g  $GENOMESIZE  | mawk -v id=$2 -v mincov=$MINCOV 'BEGIN{FS="\t";OFS="\t";}{$1=$1":Neg";$2=$2+1; if($3>=mincov){print $0,id}}'
