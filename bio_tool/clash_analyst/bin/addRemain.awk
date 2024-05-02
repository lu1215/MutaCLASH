BEGIN{
    FS=","
    OFS=","
}
{
    if(NR==1){
        print $0,"remain0","remain_seq","remain_len","hybrid_len"
    }else{
        split($rem_col,pos,"-");
        len=pos[2]-pos[1]+1
        seq=substr($hyb_col,pos[1],len)
        print $0,"remain"NR,seq,len,length($hyb_col)
    }

}
