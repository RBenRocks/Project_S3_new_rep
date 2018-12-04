import pandas as pd
import shutil

def New_move_f(DB, tissue):
    df = pd.read_csv("../data/{}.txt".format(tissue), sep="\t")
    flist = df["Array Data File"].tolist()
    for file in flist:
        shutil.copy("{}{}".format(DB,file), "../data/{}".format(tissue))

if __name__ == "__main__":
    New_move_f("../data/raw/", "liver")
    