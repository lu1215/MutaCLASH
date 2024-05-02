BEGIN{
    FS="\t"
    OFS=","
    print "hyb_seq","regulator0","reg_hyb_target_pos","on_reg_pos","transcript0","remain_pos","rem_tran_target_pos"
}
{
    if (index($4,"tran")){
        print $2,$10,$11"-"$12,$13"-"$14,$4,$5"-"$6,$7"-"$8
    }else{
        print $2,$4,$5"-"$6,$7"-"$8,$10,$11"-"$12,$13"-"$14
    }

}

