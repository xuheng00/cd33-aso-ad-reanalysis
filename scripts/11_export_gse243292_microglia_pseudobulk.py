#!/usr/bin/env python3
import csv, gzip, os, shutil, sys
import numpy as np
import h5py

root = os.path.abspath(os.path.join(os.getcwd(), "AD_CD33_reanalysis"))
src = os.path.join(root, "data", "human_raw", "GSE243292_Microglia_GEO.h5ad.gz")
tmp = os.path.join("/tmp", "GSE243292_Microglia_GEO.h5ad")
with gzip.open(src, "rb") as fi, open(tmp, "wb") as fo:
    shutil.copyfileobj(fi, fo)

with h5py.File(tmp, "r") as f:
    categories = f["obs/__categories"]
    def categorical(name):
        codes = f[f"obs/{name}"][:]
        cats = [x.decode() if isinstance(x, bytes) else str(x) for x in categories[name][:]]
        return np.array([cats[int(i)] if int(i) >= 0 else "" for i in codes], dtype=object)
    sample = categorical("sampleID")
    atscore = categorical("atscore")
    apoe = categorical("apoe")
    trem2 = categorical("trem2")
    genes = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/genename"][:]]
    unique_samples = sorted(set(sample), key=lambda x: int(x))
    sample_index = {s:i for i,s in enumerate(unique_samples)}
    pseudobulk = np.zeros((len(genes), len(unique_samples)), dtype=np.float64)
    data = f["X/data"][:]
    indices = f["X/indices"][:]
    indptr = f["X/indptr"][:]
    for row, sid in enumerate(sample):
        lo, hi = int(indptr[row]), int(indptr[row+1])
        np.add.at(pseudobulk[:, sample_index[sid]], indices[lo:hi], data[lo:hi])
    meta = []
    for sid in unique_samples:
        ix = np.where(sample == sid)[0]
        def one(values):
            vals = sorted(set(values[ix]))
            return ";".join(vals)
        meta.append([sid, one(atscore), one(apoe), one(trem2), str(len(ix))])

outdir = os.path.join(root, "data", "human_processed")
os.makedirs(outdir, exist_ok=True)
with gzip.open(os.path.join(outdir, "GSE243292_microglia_pseudobulk_counts.csv.gz"), "wt", newline="") as fh:
    w = csv.writer(fh); w.writerow(["gene"] + unique_samples)
    for gene, values in zip(genes, pseudobulk):
        w.writerow([gene] + [f"{x:.8g}" for x in values])
with open(os.path.join(outdir, "GSE243292_microglia_sample_metadata.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sampleID","atscore","apoe","trem2","n_nuclei"]); w.writerows(meta)
print("GSE243292 microglia pseudobulk exported")
