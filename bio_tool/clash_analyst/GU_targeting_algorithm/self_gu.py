import pandas as pd
import sys
import numpy as np
from tran_fa import judge,gu_fun

inp = sys.argv[1]
out = inp.replace('.csv','_gu.csv')
data=pd.read_csv(inp)
new = data.apply(lambda x: pd.Series(judge(x['regulator_seq'],x['transcript_seq'],x['rem_tran_target_pos'])),axis=1)
new.columns = ['GU_target_position','GU_targeting_score','xGU_inseed','GU_inseed','xGU_outseed','GU_outseed','total_mismatch','xGU_misPos','GU_misPos']
data = data.join(new)
data.to_csv(out,index=0)

