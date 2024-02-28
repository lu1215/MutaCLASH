if [ $5 = "region" ] || [ $5 = "site" ] || [ $5 = "up" ] || [ $5 = "abu" ]
then
    python abundance.py --inputname $1 --reg $2 --trans $3 --ex $4 --region_type $5
fi