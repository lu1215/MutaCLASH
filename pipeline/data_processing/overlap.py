import pandas as pd
import time
import argparse
from collections import Counter
import re

parser = argparse.ArgumentParser()
parser.add_argument("--inputname", help="type input filename", type=str)
args = parser.parse_args()
inputname = args.inputname

data = pd.read_csv('tmp/{}'.format(inputname))
print(len(data))
over_list = []
for s1, s2 in zip(data['remain_pos'], data['reg_hyb_target_pos']):
    s1 = s1.split('-')
    s2 = s2.split('-')
    start1 = int(s1[0])
    stop1 = int(s1[1])
    start2 = int(s2[0])
    stop2 = int(s2[1])
    if start2 <= stop1 and stop2 >= start1:
        center1 = (stop1+start1)/2
        center2 = (stop2+start2)/2
        if center1 < center2:
            over_list.append('{}-{}'.format(str(start2), str(stop1)))
        elif center2 < center1:
            over_list.append('{}-{}'.format(str(start1), str(stop2)))
        else:
            print(start1, stop1, center1)
            print(start2, stop2, center2)
            break
            #over_list,append('')
    else:
        over_list.append('0')
data['overlap'] = over_list
outputname = inputname.split('/')[-1].replace('.csv', '')
data.to_csv('tmp/{}_with_overlap.csv'.format(outputname), index=False)
