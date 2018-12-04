import pandas as pd
import csv

def Data_for_two(infile, tissue1, tissue2):
    """
    Returns metadata for drugs bresent in both kidney and liver
    """
    df = pd.read_csv(infile, sep="\t")
    drugs_t1 = Drug_for_tissue(df, tissue1)
    drugs_t2 = Drug_for_tissue(df, tissue2)
    common_drugs = Intercept(drugs_t1, drugs_t2)
    df1 = df[df["Parameter Value[Compound]"].isin(common_drugs)]
    df2 = df1[df1["Parameter Value[DoseLevel]"].isin(["Control", "High"])]
    df3 = df2[df2["Parameter Value[TimeOfSacrifice]"].isin(["24"])]
    Write_to_files(df3, tissue1, tissue2)


def Drug_for_tissue(df, tissue):
    """
    list of drugs for which complete the data is present in database
    """
    df1 = df[df["Characteristics[CellType]"].isin([tissue])]
    return set(df1["Parameter Value[Compound]"])
    

def Intercept(list1, list2):
    lst3 = [value for value in list1 if value in list2] 
    return lst3 

def Write_to_files(df, tissue1, tissue2):
    df1 = df[df["Characteristics[CellType]"].isin([tissue1])]
    df2 = df[df["Characteristics[CellType]"].isin([tissue2])]
    df1.to_csv("../data/{}.txt".format(tissue1), sep="\t", header=True)
    df2.to_csv("../data/{}.txt".format(tissue2), sep="\t", header=True)



if __name__ == "__main__":
    Data_for_two("../new_metadata.txt", "liver", "kidney")