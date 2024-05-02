import pandas as pd

def save(df1,filt,jobID,way):
    fold5=df1["RNAup_score"] <= filt
    pdata_f5=df1[(fold5)]
    out=pdata_f5.loc[:,["regulator_name","transcript_name"]]

    out.to_csv("{}/{}_up_{}.tab".format(jobID,way,abs(filt)),sep="\t",header=0,index=0)

def Affect(jobID,way,filterArray):
    df1 = pd.read_csv("/home/bba753951/Django/master_project/media/uploadfile/{}/{}/hyb_file_step5.csv".format(jobID,way),usecols=["regulator_name","transcript_name","RNAup_score"])

    # df=df.drop_duplicates(subset=None, keep='first', inplace=False)
    for i in filterArray:
        save(df1,i,jobID,way)


for i in range(87,93):
    Affect("GSM12194{}".format(i),"pir",[-5,-10,-15])
    Affect("GSM12194{}".format(i),"clan",[-5,-10,-15])
    Affect("GSM12194{}".format(i),"hyb",[-5,-10,-15])

