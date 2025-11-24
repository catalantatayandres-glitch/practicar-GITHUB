from io import StringIO
import pandas as pd
import numpy as np

TSV_data = """chrom   pos     ref     alt     gene
chr1    1050    A       G       BRCA1
chr1    2040    C       T       BRCA1
chr1    2040    C       T       BRCA1
chr2    3300    G       A       TP53
chr2    3310    T       C       TP53
chr2    5000    A       T       MYC
"""

StringIO(TSV_data)



df= pd.read_csv(StringIO(TSV_data),sep="\s+")
def mut_type(row):
    mutation_pos= {("A","G"),("G","A"),("C","T"),("T","C")}
    if (row["ref"],row["alt"]) in mutation_pos:
        return "Transition"
    else:
        return "Transversion"
    
df["mutation_type"]= df.apply(mut_type,axis=1)
diff=df.groupby("gene")["pos"].apply(lambda x: x.sort_values().diff())
diff_mean=diff.groupby(level=0).mean()


df["diff_mean"]= df.apply(lambda x: diff_mean[x["gene"]],axis=1)

diff_mean = diff_mean.rename("DIFF_MEAN").reset_index()
diff_mean.columns = ["GENE", "DIFF_MEAN"]

df["is_hotspot"] = df["pos"].duplicated(keep=False)

print(df["mutation_type"].value_counts())
