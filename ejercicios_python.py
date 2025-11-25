from io import StringIO
import pandas as pd
import numpy as np

TSV_data = """chrom	pos	ref	alt	gene	qual	depth
chr1	1050	A	G	BRCA1	42	30
chr1	2040	C	T	BRCA1	12	80
chr1	3000	G	A	BRCA1	55	40
chr2	3300	G	A	TP53	18	100
chr2	3310	T	C	TP53	90	50
chr3	5000	A	T	MYC	5	200
"""
StringIO(TSV_data)
df= pd.read_csv(StringIO(TSV_data),sep="\t")
df=df[(df["qual"]>=20)&(df["depth"]>=20)]
MUT = {("A","G"),("G","A"),("C","T"),("T","C")}
df["transition"]=df.apply(lambda row: (row["ref"],row["alt"]) in MUT,axis=1)

dicc_df=df.groupby("gene")["pos"].apply(list).to_dict()

df["mut_count"]=df["gene"].map(df["gene"].value_counts())
print(df)

