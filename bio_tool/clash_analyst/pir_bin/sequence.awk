#!/bin/awk -f

BEGIN {
	OFS=","
    btn=0
}
{
    if ( btn == 1){
        print $0 
        btn++
    }else if(btn == 3){
        btn=0
    }else{
        btn++
    }

}
#/^@/
