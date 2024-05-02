name=$1
f=$2
#for csv
#grep -q ,$name, $f

#for tsv
grep -q "^$name"$'\t' $f

find=$?
#echo $find
if [ $find -eq 1 ];then
    # have find
    echo $name
fi 
