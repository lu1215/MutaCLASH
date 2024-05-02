#!/bin/awk -f

BEGIN {
	OFS=","
}
{
    if(NR==1){
        next
    }

    key_arr[$key]=key_arr[$key]"/"$value
    value_sum[$key]+=1
}
END{
    for(i in value_sum){
        print i,value_sum[i],key_arr[i] 
        }
}
