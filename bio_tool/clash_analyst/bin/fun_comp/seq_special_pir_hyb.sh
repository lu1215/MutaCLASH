# (hyb or pir) not clan
cat pir.tab hyb.tab|sort -u|comm -23 - <(sort -u clan.tab)|wc -l 
# all three
cat pir.tab hyb.tab clan.tab|sort -u|wc -l 
# (hyb and pir ) not clan
comm -12 <(cat pir.tab|sort -u) <(cat hyb.tab|sort -u)| comm -23 - <(cat clan.tab|sort -u)|wc -l 


#three common
#comm -12 <(cut -d$'\t' -f 2,3 pir.tab|sort -u) <(cut -d$'\t' -f 2,3 hyb.tab|sort -u)| comm -12 - <(cut -d$'\t' -f 2,3 clan.tab|sort -u)|wc -l 
