
#!/bin/awk -f
BEGIN {
	OFS=","
}
{
    if(NR==1){
        next
    }
    print $hyb,"("$reg":"$tran")"

}
