from io import StringIO
import pandas as pd
import numpy as np

TSV_data = """chrom	pos	ref	alt	gene	quality
chr1	1050	A	G	BRCA1	42
chr1	1050	A	G	BRCA1	42
chr1	2040	C	T	BRCA1	12
chr1	3000	G	A	BRCA1	55
chr2	3300	G	A	TP53	18
chr2	3310	T	C	TP53	90
chr3	5000	A	T	MYC	5
"""
MUT= {("A","G"),("G","A"),("C","T"),("T","C")}

StringIO(TSV_data)

df= pd.read_csv(StringIO(TSV_data),sep="\t")
df=df[df["quality"]>=40]
df["mut_type"]=df.apply(lambda row: "TRANSITION" if (row["ref"],row["alt"]) in MUT else "TRANSVERSION",axis=1)

df["quality_mean"]= df["quality"].mean()
transition_count= df[df["mut_type"]=="TRANSITION"].groupby("gene").size()
transversion_count= df[df["mut_type"]=="TRANSVERSION"].groupby("gene").size()
df["is_hostpot"]=df["pos"].duplicated(keep=False)
print(df)

mut_pos_dicc= {}


def create_dicc(row):
    pos=row["pos"]
    gen=row["gene"]
    if gen not in mut_pos_dicc:
        mut_pos_dicc[gen]=[]
    mut_pos_dicc[gen].append(pos)


df.apply(create_dicc,axis=1)
print(mut_pos_dicc)


