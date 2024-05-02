#!/bin/awk -f


BEGIN {
	OFS=","
    print "Gene_name","counts_of_transcripts","transcript_name","counts_of_regulators","regulator_name"

}
FILENAME==ARGV[1]{
    gene_RNA[$1]=$3
}
FILENAME==ARGV[2]{
    split($3,RNA,"/")
    value=""
    for(i in RNA){
        if(RNA[i]!=""){
            if(gene_RNA[RNA[i]]){
                value=value""gene_RNA[RNA[i]]
            }
        }
    }
    if(value){
        end_value=""
        count=0
        split(value,value_arr,"/")
        for(j in value_arr){
            if(value_arr[j])
                unique_arr[value_arr[j]]++
        }
        for(k in unique_arr){
            end_value=end_value"/"k
            count++

        }
        
        print $1,$2,$3,count,end_value
        delete value_arr
        delete unique_arr
    }else{

      
        print $1,$2,$3,0,0

     }

}
END{
}
