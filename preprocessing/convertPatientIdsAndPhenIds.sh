#!/bin/bash
awk 'BEGIN{FS=","; OFS=","}NR == 1; NR > 1 {print $0 | "sort -n -t, -k2,2"}' $1 \
  | awk 'BEGIN{FS=","; OFS=","; i=-1; last =-1;} {if(NR==1){print $0; next;}if($2=="NA"){next;} if(last!=$2){last=$2; ++i} $2=i; print $0;}'  \
  | awk  'BEGIN{FS=","; OFS=","} NR == 1; NR > 1 {print $0 | "sort -n -t, -k1,1"}'  \
  | awk 'BEGIN {FS=",";OFS=","; i=-1; last=""}{if(NR==1){print $0;next;} if(last!=$1) {last=$1;++i} $1=i; print $0;} ' \
  > $2