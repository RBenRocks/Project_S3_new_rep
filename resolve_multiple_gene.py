import pandas as pd

def Main(infile):
    """
    solves problem of multiple gene symbols
    """
    df = pd.read_csv(infile, sep="\t")
    for i, row in df.iterrows():
        #entry = df[i, row]
        entry = str(row['Gene Symbol'])
        if "///" in entry:
            symbol = parse(entry)
        else:
            if entry == "nan":
                symbol = ""
            else:
                symbol = entry
        df.at[i,'Gene Symbol'] = symbol.upper()
    df.to_csv("../data/Array_resolved.txt", sep="\t", header=True)

def parse(entry):
    """
    entry a string to be modified(chosen)
    """
    entry = entry.split("///")
    l = []
    for i in entry:
        j = i.replace(" ", "")
        l.append(j)
    k = []
    for i in l:
        if not i.startswith("LOC"):
            k.append(i)
    if len(k) > 0:
        return k[0]
    else:
        return ""
        
if __name__ == "__main__":
    #print(parse("LOC001 /// gto9"))
    print(Main("../data/Array_design.txt"))
    