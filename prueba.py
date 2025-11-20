import pandas as pd
import numpy as np
fasta= """>GeneA
ATGCGTACGTAGCTAGCTAG
>GeneB
TTTACGATGCGTATATATCG
>GeneC
CGCGCGATATATATATGCGC
"""
def read_fasta(fasta_file):
    lines= fasta_file.strip().split("\n")
    secuencias= {}
    name=""
    seq=""

    for line in lines:
        if line.startswith(">"):
            if name != "":

                secuencias[name]=seq
            name= line[1:]
            seq=""

        else:
            seq+= line.upper()
    secuencias[name]=seq
    return secuencias

data= read_fasta(fasta)
print(data)

df= pd.DataFrame({
    "gene_id": list(data.keys()),
    "sequence": list(data.values())
}
)
print(df)

df["sequence"]= df["sequence"].apply(lambda seq: np.array(list(seq)))
df["GC_content"]= df["sequence"].apply(lambda seq: ((seq=="G").sum() + (seq=="C").sum()/seq.size *100))
df["length"]= df["sequence"].apply(len)
df["ATG_times"]=df["sequence"].apply(lambda seq: "".join(seq).count("ATG"))
print(df)