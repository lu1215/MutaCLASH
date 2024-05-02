import pandas as pd
import matplotlib.pyplot as plt

def drawAllCount():
    df = pd.read_csv("all_count.txt", header=None,index_col=0) 
    df.columns=['pir', 'hyb','clan']

    print(df)

    plt.figure()
    df.iloc[0:6,:].plot(kind='bar',colormap='cividis')
    plt.savefig("human_count.png")
    plt.figure()
    df.iloc[6:,:].plot(kind='bar',colormap='cividis')
    plt.savefig("elegan_count.png")

def drawHits(columns,kind):
    df = pd.read_csv("hit_{}_count.csv".format(kind), header=None,index_col=0) 

    df.columns=columns

    print(df)

    plt.figure()
    df.iloc[0:6,:].plot(kind='bar',colormap='cividis')
    plt.xlabel(kind+"(human)")
    plt.ylabel("hits percent(%)")
    plt.savefig("human_hit_{}.png".format(kind))

    plt.figure()
    df.iloc[6:,:].plot(kind='bar',colormap='cividis')
    plt.xlabel(kind+"(c.elegans)")
    plt.ylabel("hits percent(%)")
    plt.savefig("elegan_hit_{}.png".format(kind))

readCount_col=["0","1","2","3","4"]
fold_col=["0","-5","-10","-15","-20"]

drawHits(readCount_col,"readCount")
drawHits(fold_col,"fold")


