#coding:utf-8
#python 3
#author: Joshua Chou

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

####
FILENAME="G7_data.xls"
GENES=["GAPDH","COL"]
NORMOLIZE_TO="18sRNA"
####

df=pd.read_excel(FILENAME)
dt=df[["Sample Name","Target Name","Cт Mean"]].fillna(0)
dt=dt.drop_duplicates()
dt.columns=["sample","target","ct"]
dt=dt[dt.loc[:,"sample"]!="NTC"]
dt=dt.reset_index()
dt=dt.drop("index",axis=1)

# calculate dct
dcts=[]
for i in range(len(dt["sample"])):
    pc_sample=str(dt.loc[i,"sample"])
    pc_target=str(dt.loc[i,"target"])
    pc_ct=float(dt.loc[i,"ct"])
    nc_sample=pc_sample
    nc_target=NORMOLIZE_TO
    nc_ct=dt[(dt["sample"]==nc_sample)&(dt["target"]==nc_target)].ct.values[0]
    dcts.append(pc_ct-nc_ct)
dt["dct"]=dcts

# calculate ddct
ddcts=[]
for i in range(len(dt["sample"])):
    pc_sample=str(dt.loc[i,"sample"])
    pc_target=str(dt.loc[i,"target"])
    pc_dct=float(dt.loc[i,"dct"])
    nc_sample="NC-1"
    nc_target=pc_target
    nc_dct=dt[(dt["sample"]==nc_sample)&(dt["target"]==nc_target)].dct.values[0]
    ddcts.append(pc_dct-nc_dct)
dt["ddct"]=ddcts

# calculate 2-ddct
dt["2-ddct"]=np.power(0.5,dt.ddct)

# calculate average
averages=[]
for i in range(len(dt["sample"])):
    pc_sample=str(dt.loc[i,"sample"])
    pc_target=str(dt.loc[i,"target"])
    pc_2ddct=float(dt.loc[i,"2-ddct"])
    average=[float(dt[(dt["target"]==pc_target)&(dt["sample"]==x)].loc[:,"2-ddct"]) for x in ["NC-1","NC-2","NC-3"]]
    averages.append(np.mean(average))
dt["average"]=averages

# calculate 2-ddct/average
dt["2-ddct/average"]=dt["2-ddct"]/dt["average"]

# filter 18sRNA
dt=dt[dt.target!=NORMOLIZE_TO]

# sort
dt=dt.sort_values("target")
dt=dt.reset_index()
dt=dt.drop("index",axis=1)

# output
dt.to_csv("output.csv",index=False)

# plot
dt["Treatments"]=[x.split("-")[0] for x in list(dt.loc[:,"sample"])]
for i in GENES:
    sns.boxplot(x="target",y="2-ddct/average",hue="Treatments",data=dt[dt.target==i],palette="muted")
    plt.xlabel("")
    plt.ylabel("(2^-ddct)/average")
    plt.show()
    plt.close()
