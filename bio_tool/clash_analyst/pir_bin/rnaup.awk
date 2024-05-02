#!/bin/awk -f

function append(seq,ch,count,frontBack){
    #print"append"
    if(frontBack==0){
        for(i=0 ; i<count ; i++){
            seq=ch""seq
        }
    }else{
        for(i=0 ; i<count ; i++){
            seq=seq""ch
        }
    }
    return seq
}

function seqCompare(seq1,seq2,pth1,pth2,sC_result,__ARGVEND__,i,seq1_len,ch1,ch2,seq2_len){
    #print "seqCompare"
    # because awk index start from 1
    i=1
    seq1_len=length(seq1)
    while(i<=seq1_len){
        ch1=substr(pth1,i,1)
        ch2=substr(pth2,i,1)
        if(ch1 != ch2){
            if( ch1 == "." ){

                seq2=substr(seq2,1,i-1)"-"substr(seq2,i)
                pth2=substr(pth2,1,i-1)"."substr(pth2,i)
            }else{
                seq1=substr(seq1,1,i-1)"-"substr(seq1,i)
                pth1=substr(pth1,1,i-1)"."substr(pth1,i)
                seq1_len++
            }
        }
        i++
    }

    # because seq2 also extend it's length
    seq2_len=length(seq2)

    if(seq1_len < seq2_len){

        seq2_len-=seq1_len
        seq1=append(seq1,"-",seq2_len,1)
    }
    sC_result[1]=seq1
    sC_result[2]=seq2
}
function reverse(seq,__ARGVEND__,rev_seq,i){
    #print "reverse"
    rev_seq=""
    i=length(seq)
    for(i ; i>0 ; i--){
        rev_seq=rev_seq""substr(seq,i,1)
    }
    return rev_seq
}

function rnaup(RNA_slice_seq, input_seq, RNA_start, RNA_seq, RNA_len, input_len, RNA_slice_len, old,__ARGVEND__,sC_result,sequence,cmd,count,var,res1,pth_seq,seq1_pos,seq2_pos,up_seq,len_btn,seq_arr,pth_arr,swap,pos_arr1,pos_arr2,etd_len1,etd_len2,lim_len1,lim_len2,etd_more1,etd_more2,up_result){
    #print "rnaup"
    sequence=RNA_slice_seq"\\&"input_seq
    #print sequence
    cmd="echo "sequence" | RNAup -o -b -d2 --noLP -c 'S' -"
    count=1

    # get RNAup result
    while(cmd|getline var){
        #print var
        if(count==1){
            split(var,res1," ")
            pth_seq=res1[1]
            seq1_pos=res1[2]
            seq2_pos=res1[4]

        }else if(count==2){
            up_seq=var
            break
        } 
        count++
    }
    close(cmd)

    if(up_seq=="&"){
        #up_result[1]="0-0"
        #up_result[2]=RNA_slice_seq
        #up_result[3]=input_seq
        print old,"100","0-0",RNA_slice_seq,input_seq
    }else{
        # judge length for RNAup may change pos
        len_btn=(RNA_slice_len < input_len) ? 1 : 0

        split(up_seq,seq_arr,"&")
        split(pth_seq,pth_arr,"&")

        # reverse seq1
        seq_arr[1]=reverse(seq_arr[1])
        pth_arr[1]=reverse(pth_arr[1])

        # change seq1 "(" to ")" for seqCompare() 
        gsub("\\(",")",pth_arr[1])

        # seqCompare(seq1,seq2,pth1,pth2,seq1_len,seq2_len,sC_result)
        seqCompare(seq_arr[1],seq_arr[2],pth_arr[1],pth_arr[2],sC_result)
        
        seq_arr[1]=reverse(sC_result[1])
        seq_arr[2]=sC_result[2]

        # if seq1 less then seq2     
        # swap pos and seq
        if(len_btn == 1){
            swap=seq_arr[1]
            seq_arr[1]=seq_arr[2]
            seq_arr[2]=swap

            swap=seq1_pos
            seq1_pos=seq2_pos
            seq2_pos=swap
        }

        split(seq1_pos,pos_arr1,",")
        split(seq2_pos,pos_arr2,",")

        #transfer relative pos to absolute pos
        #print RNA_start
        pos_arr1[1]=RNA_start+pos_arr1[1]-1
        pos_arr1[2]=RNA_start+pos_arr1[2]-1

        # extend input_seq to full
        etd_len1=pos_arr2[1]-1
        etd_len2=input_len-pos_arr2[2]
        
        # calculate the limits that can be extended
        lim_len1=pos_arr1[1]-1
        lim_len2=RNA_len-pos_arr1[2]

        # calculate the length can't not extend more (for front and back gap "-") 
        etd_more1=0
        etd_more2=0
        if(lim_len1 < etd_len2){
            etd_more1=etd_len2-lim_len1
            etd_len2=lim_len1
        }
        if(lim_len2 < etd_len1){
            etd_more2=etd_len1-lim_len2
            etd_len1=lim_len2
        }

        # extend seq1 front pos
        #print "---"pos_arr1[1],etd_len2
        pos_arr1[1]=pos_arr1[1]-etd_len2
        # pos_arr1[2]=pos_arr1[2]+etd_len1

        #extend seq2 fron pos
        pos_arr2[1]=pos_arr2[1]-etd_len1

        #print "---"pos_arr1[1]
        #print substr(RNA_seq,pos_arr1[1],etd_len2)
        seq_arr[1]=substr(RNA_seq,pos_arr1[1],etd_len2)""seq_arr[1]""substr(RNA_seq,pos_arr1[2]+1,etd_len1)
        # etd_len1(now)+etd_more2=etd_len1(origin)
        seq_arr[2]=substr(input_seq,1,etd_len1+etd_more2)""seq_arr[2]""substr(input_seq,pos_arr2[2]+1)


       # transfer U to T 
        gsub("U","T",seq_arr[1])
        gsub("U","T",seq_arr[2])


        # extend seq1 back pos
        pos_arr1[2]=pos_arr1[2]+etd_len1
        
        # append(seq,ch,count,frontBack)
        # append seq1 that can't expend more by "-"
        #print etd_more1,pos_arr1[1],"-------"
        if(etd_more1 > 0){
            seq_arr[1]=append(seq_arr[1],"-",etd_more1,0)
            pos_arr1[1]=pos_arr1[1]-etd_more1
        }
        if(etd_more2 > 0){
            seq_arr[1]=append(seq_arr[1],"-",etd_more2,1)
            pos_arr1[2]=pos_arr1[2]+etd_more2
        }

        up_result[1]=pos_arr1[1]"-"pos_arr1[2]
        up_result[2]=seq_arr[1]
        up_result[3]=seq_arr[2]

        rnaup_score(RNA_seq, RNA_len, pos_arr1[1], pos_arr1[2], input_seq, old,up_result)

    }
}

function rnaup_score(RNA_seq, RNA_len, RNA_start, RNA_end, input_seq, old,up_result,__ARGVEND__,sequence,cmd,count,res,seq1,seq2){
    #print "rnaup_score"

    RNA_start=(RNA_start < 1)? 1 : RNA_start
    RNA_end=(RNA_end > RNA_len)? RNA_len : RNA_end
    slice_RNA=substr(RNA_seq,RNA_start,RNA_end-RNA_start+1)


    sequence=slice_RNA"\\&"input_seq
    cmd="echo "sequence" | RNAup -o -b -d2 --noLP -c 'S' -"
    count=1

    # get RNAup result
    cmd|getline var
    close(cmd)
    split(var,res," ")
    up_result4=substr(res[5],2)


    # remove bulge at regulator front and end
    # remember seq2 need to be reverse
    seq1=up_result[2]
    seq2=up_result[3]

    front=index(seq2,"-")
    while(front == 1){
        seq1=substr(seq1,1,length(seq1)-1)
        seq2=substr(seq2,2)
        front=index(seq2,"-")
    }

    back=substr(seq2,length(seq2))
    while (back == "-"){
        seq1=substr(seq1,2)
        seq2=substr(seq2,1,length(seq2)-1)
        back=substr(seq2,length(seq2))
    }


    print old,up_result4,up_result[1],seq1,seq2

}


function etd_seq(RNA_seq, input_seq, remain_seq, remain_RNA_pos, RNA_len, input_len, remain_len, old,__ARGVEND__,pos,pos1,pos2,final_remain_len,final_remain_seq,etd1,etd2,lim_len1,lim_len2){
    #print "etd_seq"
    split(remain_RNA_pos,pos,"-")
    pos1=pos[1]
    pos2=pos[2]
    final_remain_seq=remain_seq
    final_remain_len=remain_len
	if (remain_len < input_len){
		etd1=input_len-remain_len
		etd2=etd1

        lim_len1=pos1-1
        lim_len2=RNA_len-pos2

        etd1=(lim_len1<etd1) ? lim_len1 : etd1
        etd2=(lim_len2<etd2) ? lim_len2 : etd2

        pos1-=etd1
        #pos2+=etd2

        final_remain_len=remain_len+etd1+etd2

        # awk string index start from "1"
        final_remain_seq=substr(RNA_seq,pos1,final_remain_len)

	}
    # call rnaup.sh
    #cmd="bash rnaup.sh rnaup TCTTGTTCATCACGACTTCACG TAAAGTCGTGATGACAAATGA 213"
    #rnaup(RNA_slice_seq,input_seq,RNA_start,RNA_seq,RNA_len,input_len,RNA_slice_len,up_result)

    rnaup(final_remain_seq,input_seq,pos1,RNA_seq,RNA_len,input_len,final_remain_len,old)
}
BEGIN {
	FS=","
	OFS=","
    
}

{
    etd_seq($RNA_seq,$input_seq,$remain_seq,$remain_RNA_pos,$RNA_len,$input_len,$remain_len,$0)
}
