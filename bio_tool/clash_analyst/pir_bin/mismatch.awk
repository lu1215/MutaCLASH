#!/bin/awk -f
{
    if($8)
        mismatch=gsub(",","",$8)+1
    else
        mismatch=0
    printf "%s,%s,%s-%s,%d\n",$1,$3,$4+1,$4+length($6),mismatch

}
