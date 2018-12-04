import pandas as pd
import os
### resolve tissue later on

def De(path):
    res = pd.DataFrame(columns=['Drug','Tissue','#p<0.05','#q<0.05','list_p','list_q'])
    for tissue in ["liver", "kidney"]:
        fnames = os.listdir("{}{}/DE".format(path, tissue))
        for f in fnames:
            drug = f.strip(".txt")
            df = pd.read_csv("{}{}/DE/{}".format(path, tissue ,f), sep="\t")
            df1 = df[df["P.Value"] < 0.05]
            n_p = df1.shape[0]
            gene_list = df['Unnamed: 0'].tolist()
            df2 = df[df["adj.P.Val"] < 0.05]
            n_q = df2.shape[0]
            gene_list_q = df2['Unnamed: 0'].tolist()
            temp = pd.DataFrame([[drug,tissue,n_p,n_q,gene_list,gene_list_q]], columns=['Drug','Tissue','#p<0.05','#q<0.05','list_p','list_q'])
            res = res.append(temp)
    res.to_csv("{}DE_tablepq.txt".format(path), sep="\t", header=True)
def Gse(path):
    res = pd.DataFrame(columns=['Drug','Tissue','up<0.05','down<0.05','nondir<0.05','list_up','list_down','list_non'])
    for tissue in ["liver", "kidney"]:
        fnames = os.listdir("{}{}/piano".format(path, tissue))
        if 'desktop.ini' in fnames:
            fnames.remove('desktop.ini')
        for f in fnames:
            drug = f.strip("_piano.txt")
            df = pd.read_csv("{}{}/piano/{}".format(path, tissue ,f), sep="\t")
            df1 = df[df["p adj (dist.dir.up)"] < 0.05]
            n_1 = df1.shape[0]
            gene_list = df['Name'].tolist()
            df2 = df[df["p adj (dist.dir.dn)"] < 0.05]
            n_2 = df2.shape[0]
            gene_list_1 = df2['Name'].tolist()
            df3 = df[df["p adj (non-dir.)"] < 0.05]
            n_3 = df3.shape[0]
            gene_list_2 = df2['Name'].tolist()
            temp = pd.DataFrame([[drug,tissue,n_1,n_2,n_3,gene_list,gene_list_1,gene_list_2]], columns=['Drug','Tissue','up<0.05','down<0.05','nondir<0.05','list_up','list_down','list_non'])
            res = res.append(temp)
    res.to_csv("{}GSE_tableq.txt".format(path), sep="\t", header=True)
    print(res['up<0.05'],res['down<0.05'],res['nondir<0.05'])

if __name__ == "__main__":
    #De("./new_dir/data/")
    Gse("./new_dir/data/")