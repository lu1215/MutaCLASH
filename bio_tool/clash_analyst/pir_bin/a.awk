BEGIN{
    seq1="--AAAGACUAUU--U-CACGGUUGAGGBB"
    seq2="---UUUCUGUA---ACGUUAACUCCUU--"


    print seq1
    print seq2
    front=index(seq2,"-")
    while(front == 1){
        seq1=substr(seq1,2)
        seq2=substr(seq2,2)
        front=index(seq2,"-")
    }

    print seq1
    print seq2

    back=substr(seq2,length(seq2))
    while (back == "-"){
        seq1=substr(seq1,1,length(seq1)-1)
        seq2=substr(seq2,1,length(seq2)-1)
        back=substr(seq2,length(seq2))
    }
    print seq1
    print seq2
}
