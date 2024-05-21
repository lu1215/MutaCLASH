if [ $1 = "chira_chimeras" ] || [ $1 = "chira_single" ]
then
    # merge regulator and target reference
    python merge_ref.py --reg $4 --tar $5

    # find "MD" tag
    samtools calmd $2 reference.fa > sorted.sam
    (printf '0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\n'; cat sorted.sam) > MD.sam

    python bwa_find.py --MD MD.sam --input $3 --ref_col transcript_name
    rm reference.fa reference.fa.fai sorted.sam MD.sam
else
    python bowtie2_find.py --inputname $3 --trans ${5%.*}.csv
fi
