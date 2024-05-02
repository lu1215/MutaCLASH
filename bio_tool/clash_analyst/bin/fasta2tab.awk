#!/usr/bin/awk -f
BEGIN{
    FS=","
    OFS="\t";
    }
/^>/{
        a=substr($1,2);
        getline;
        print a,$1;
    }

