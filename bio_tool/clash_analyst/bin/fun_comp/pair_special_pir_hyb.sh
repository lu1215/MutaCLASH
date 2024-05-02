# (hyb or pir) not clan
cat pir.tab hyb.tab|cut -d$'\t' -f 2,3|sort -u|comm -23 - <(cut -d$'\t' -f 2,3 clan.tab|sort -u)|wc -l 
# all three
cat pir.tab hyb.tab clan.tab|cut -d$'\t' -f 2,3|sort -u|wc -l 
# (hyb and pir ) not clan
comm -12 <(cut -d$'\t' -f 2,3 pir.tab|sort -u) <(cut -d$'\t' -f 2,3 hyb.tab|sort -u)| comm -23 - <(cut -d$'\t' -f 2,3 clan.tab|sort -u)|wc -l 


#three common
#comm -12 <(cut -d$'\t' -f 2,3 pir.tab|sort -u) <(cut -d$'\t' -f 2,3 hyb.tab|sort -u)| comm -12 - <(cut -d$'\t' -f 2,3 clan.tab|sort -u)|wc -l 
