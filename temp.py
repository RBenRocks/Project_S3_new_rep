import pandas as pd
import os
### resolve tissue later on

def De(path):
    res = pd.DataFrame(columns=['Drug','Tissue','q<0.05'])
    for tissue in ["liver", "kidney"]:
        fnames = os.listdir("{}{}/DE".format(path, tissue))
        for f in fnames:
            drug = f.strip(".txt")
            df = pd.read_csv("{}{}/DE/{}".format(path, tissue ,f), sep="\t")
            df2 = df[df["adj.P.Val"] < 0.05]
            n_q = df2.shape[0]
            temp = pd.DataFrame([[drug,tissue,n_q]], columns=['Drug','Tissue','q<0.05'])
            res = res.append(temp)
        fnames = os.listdir("{}/common_de/".format(path))
    for f in fnames:
        drug = f.strip("_intersect.txt")
        df = pd.read_csv("{}/common_de/{}".format(path,f), sep="\t")
        n_q = df.shape[0]
        temp = pd.DataFrame([[drug,"intersect",n_q]], columns=['Drug','Tissue','q<0.05'])
        res = res.append(temp)
    res.to_csv("{}DE_tablerep.txt".format(path), sep="\t", header=True, index=False)
def Gse(path):
    res = pd.DataFrame(columns=['Drug','Tissue','up<0.05','down<0.05','nondir<0.05'])
    for tissue in ["liver", "kidney"]:
        fnames = os.listdir("{}{}/piano/".format(path, tissue))
        if 'desktop.ini' in fnames:
            fnames.remove('desktop.ini')
        for f in fnames:
            print(f)
            drug = f.strip("_piano.txt")
            df = pd.read_csv("{}{}/piano/{}".format(path, tissue ,f), sep="\t")
            df1 = df[df["p adj (dist.dir.up)"] < 0.05]
            n_1 = df1.shape[0]
            df2 = df[df["p adj (dist.dir.dn)"] < 0.05]
            n_2 = df2.shape[0]
            df3 = df[df["p adj (non-dir.)"] < 0.05]
            n_3 = df3.shape[0]
            temp = pd.DataFrame([[drug,tissue,n_1,n_2,n_3]], columns=['Drug','Tissue','up<0.05','down<0.05','nondir<0.05'])
            res = res.append(temp)
            print(drug)
    fnames = os.listdir("{}/5padj/UP_gse/".format(path))
    for f in fnames:
        drug = f.strip(".txt")
        df = pd.read_csv("{}/5padj/UP_gse/{}".format(path, f), sep="\t")
        n_1 = df.shape[0]
        df = pd.read_csv("{}/5padj/DOWN_gse/{}".format(path, f), sep="\t")
        n_2 = df.shape[0]
        df = pd.read_csv("{}/5padj/NDIR_gse/{}".format(path, f), sep="\t")
        n_3 = df.shape[0]
        temp = pd.DataFrame([[drug,"intersect",n_1,n_2,n_3,]], columns=['Drug','Tissue','up<0.05','down<0.05','nondir<0.05'])
        res = res.append(temp)
    res.to_csv("{}GSE_tabl_rep.txt".format(path), sep="\t", header=True, index=False)

if __name__ == "__main__":
    De("../data/")
    #Gse("../data/")