#!/bin/awk -f

BEGIN {
	FS=","
	OFS=","
}
{
	if (NR==1){
		# changer reference0,input0 to reference_id,input_id	
		#gsub("0","_id",$0)
		print $0,"remain0","remain_seq","remain_len","remain_pos" 

		#print "remain0","sequence","remain_len" > remain_result 
		#close(remain_result)
	}else{

		split($target_pos,pos,"-")
		if(pos[1]-1 > $ref_len-pos[2]){
			remain_len=pos[1]-1
			remain_seq=substr($ref_seq,1,remain_len)
			remain_pos="1-"(remain_len)
		}else{
			remain_len=$ref_len-pos[2]
			remain_seq=substr($ref_seq,pos[2]+1)
			remain_pos=(pos[2]+1)"-"$ref_len

		}
        if(remain_len >= len){
            print $0,"remain"(NR-1),remain_seq,remain_len,remain_pos 
            # for step 2 remain.fasta
            #printf ">%s\n%s\n","remain"(NR-1),remain_seq >> remain_fasta
            #print "remain"(NR-1),remain_seq,remain_len >> remain_result 
            #close(remain_result)
        }
	}
}
