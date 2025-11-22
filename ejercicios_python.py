from io import StringIO
import pandas as pd
import numpy as np

TSV_data = """gene_id	sequence	ref_base	mut_base
G1	ATGCCCG	A	T
G2	ATGCGTATCG	C	G
G3	TTTACGTT	A	A
G4	CGCGCGATAT	T	C
G5	TGATGATGA	G	A
"""

StringIO(TSV_data)
transitions = {("A","G"), ("G","A"), ("C","T"), ("T","C")}
df = pd.read_csv(StringIO(TSV_data), sep="\t")
df["sequence"]= df["sequence"].apply(lambda seq: np.array(list(seq)))
df["length"]= df["sequence"].apply(len)
df["pos_ref"]= [np.where(seq==ref)[0] for seq,ref in zip(df["sequence"],df["ref_base"])]
df["mutation_possible"]= [any(seq[pos]==mut for pos in ref_pos) for seq,mut,ref_pos in zip(df["sequence"],df["mut_base"],df["pos_ref"]) ]
df["mutation_type"]= df.apply(lambda row: "synonymous"
                                if row["ref_base"] == row["mut_base"] 
                                else "transition"
                                if (row["ref_base"],row["mut_base"])in transitions
                                else "transversion",
                                axis=1)
df["GC_content"] = df["sequence"].apply(lambda seq: ((seq=="G").sum() + (seq=="C").sum()) / seq.size * 100)
mean_len= np.mean(df["length"])
std_len= np.std(df["length"])
mean_gc=np.mean(df["GC_content"])
std_gc= np.std(df["GC_content"])
df["norm_len"]=(df["length"] - mean_len) /std_len
df["gc_norm"]= (df["GC_content"]-mean_gc)/std_gc
df["gene_interesant"]= df.apply(lambda row: "YES"
                                if "".join(row["sequence"]).startswith("ATG")
                                 or row["mutation_possible"]
                                 or row["GC_content"] > mean_gc
                                else "NO",
                                axis=1)

def escoger_genes(row):
    seq= "".join(row["sequence"])
    if seq.startswith("ATG"):
        return "YES"
    elif row["mutation_possible"]:
        return "YES"
    elif row["GC_content"] > mean_gc:
        return "YES"
    else:
        return "NO"
    
df["gene_interesantee"]= df.apply( escoger_genes,axis=1)

print(df.columns)
print(df["gene_interesant"].unique())

print(df.groupby("gene_interesant").size())
print(df.groupby("gene_interesant")["GC_content"].mean())

