# merge regulator and target reference
python merge_ref.py --reg $3 --tar $4

# find "MD" tag
samtools calmd $1 reference.fa > sorted.sam
(printf '0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\n'; cat sorted.sam) > MD.sam

python bwa_find.py --MD MD.sam --input $2 --ref_col transcript_name
rm reference.fa reference.fa.fai sorted.sam MD.sam
