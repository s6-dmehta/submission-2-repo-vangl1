# bio312-vangl1-ntd
Analysis for VANGL1 variant–conservation paper (BIO 312): alignment, conservation scoring, enrichment tests, figures

# VANGL1 Neural Tube Defect (NTD) Variant Analysis  
BIO 312 – Molecular Evolution & Bioinformatics  
Dalisha Mehta

## Overview
This project investigates evolutionary conservation, domain structure, and variant distribution in the VANGL1 gene, a planar cell polarity protein implicated in neural tube defects (NTDs). 

All analyses were performed on AWS EC2 using Python, pandas, matplotlib, and scripts.

---

## Project Objectives
1. **Characterize VANGL1 domain architecture** using curated transmembrane and C-terminal annotations.  
2. **Quantify evolutionary conservation** across the full protein and within functional regions.  
3. **Compare constraint at NTD-associated vs. benign variants** using gnomAD data.  
4. **Assess enrichment** of NTD variants in high-constraint residues and the PDZ-binding motif.  
5. **Visualize VANGL1 structure** using AlphaFold-predicted 3D coordinates.  
6. **Reconstruct phylogeny** of VANGL1 orthologs to place sequence patterns in evolutionary context.

## 1.1  Start the lab, make sure your instance is running on EC2 and log in via ssh.


## 1.2 Clone Lab 13

On the command line, clone the repository.

```bash
git clone git@github.com:s6-dmehta/submission-2-repo-vangl1.git
```
This command downloads a full copy of the repository from your GitHub account onto the current machine using SSH authentication.

```bash
cd submission-2-repo-vangl1
```
This moves you into the project directory so that all later commands run inside the cloned repository.

## 1.3 Create Folders 

```bash
mkdir -p figures tables scripts 
```
This command creates all 3 directories at once and the -p flag ensures the command doesn’t error if a folder already exists.

## 1.4 Existing data files and their summaries

### `VANGL1_aligned.fasta`
**What it is:** Multiple sequence alignment of VANGL1 homologs in FASTA format.  
**How I got it:** I produced this file in Lab 4 by taking BLASTP homolog hits (from Lab 3), filtering them, and aligning them using MAFFT. The original file was `VANGL1.homologs.al.fas` and I renamed it for clarity.

---

### `VANGL1_alignment.phy`
**What it is:** The same alignment as above, but in PHYLIP format.  
**How I got it:** Created during Lab 5 by converting the MAFFT alignment to PHYLIP, removing duplicate sequences, and preparing it for phylogenetic inference. The original file was `uniqueseq.phy`, renamed here.

---

### `VANGL1_conservation.tsv`
**What it is:** Cleaned per-residue conservation scores for human VANGL1.  
**How I got it:** I imported output into Python and extracted only the numeric columns needed for Figures 1 and 4.

---

### `VANGL1_domains.tsv`
**What it is:** Domain annotation table for VANGL1.  
**How I got it:** I manually looked up VANGL1 domain boundaries on UniProt, then typed the coordinates into a TSV file using `nano`.

---

### `VANGL1_variants.tsv`
**What it is:** Combined table of pathogenic NTD variants and benign population variants.  
**How I got it:**  
- I manually extracted NTD variants from the literature used for the final paper.  
- I used the VANGL1 gnomAD missense variants, filtered for coding changes, and saved them in TSV format.

---

### `VANGL1_AF2.pdb`
**What it is:** AlphaFold2 3D structural model of human VANGL1.  
**How I got it:** **How I got it:** In Lab 6, I mapped VANGL1 RefSeq accessions to UniProt IDs using the provided `refseq_uniprot.tsv` file, then ran Josh’s `afetch_by_uniprot.sh` script on my VANGL1 mapping file. This script downloaded the AlphaFold/SWISS-MODEL PDBs into `lab06-$MYGIT/VANGL1/`. I copied the human VANGL1 model into this repo and renamed it `VANGL1_AF2.pdb`.

 ---

### `VANGL1.treefile`
**What it is:** Maximum-likelihood phylogenetic tree for VANGL1/VANGL2 homologs.   
**How I got it:** I ran IQ-TREE on the VANGL1 homolog alignment from Lab 5 to infer a maximum-likelihood tree, then copied the resulting `VANGL1.treefile` into the `data/` folder for this reproducible repo.


## 2. Creating Figures (Figure 1)

## 2.1 Creating the File
In the terminal 

```bash
cd ~/submission-2-repo-vangl1/scripts
cat > fig1.py
```
cat: opens a blank input so whatever you paste gets written out.
fig1.py: this is the file where your pasted Python code will be saved.

Paste these code blocks in order and make sure to not click enter (I've provided a large code block at the end of this section to make it much easier and efficient, use that if you don't want to go through each chunk)

### 2.2 Create the script

```bash
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
```
pandas: loads table data.
plt: makes the plot.
Path: handles file paths.

```bash
DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)
```
DATA_DIR: where TSV files live.
OUTPUT_DIR: where the PNG will be saved.

```bash
cons_file = DATA_DIR / "VANGL1_conservation.tsv"
domains_file = DATA_DIR / "VANGL1_domains.tsv"
```
These point to the specific TSV files needed to create Figure 1.

```bash
cons_df = pd.read_csv(cons_file, sep="\t")
pos_col = cons_df.columns[0]
score_col = cons_df.columns[1]
```
Loads conservation scores and identifies the position and score columns.

```bash
cons_df = cons_df.sort_values(pos_col)
cons_df["smooth"] = (
    cons_df[score_col]
    .rolling(window=7, center=True, min_periods=1)
    .mean()
)
```
Sorts by residue position and creates a 7-residue smoothed trend line.

```bash 
domains_df = pd.read_csv(domains_file, sep="\t")
start_col = "start"
end_col = "end"
name_col = "name"
```
Loads TM/PDZ positions from the domain file.

```bash
plt.style.use("default")
fig, ax = plt.subplots(figsize=(10,4))
```
Creates a clean plotting canvas.

```bash 
tm_rows = domains_df[domains_df[name_col].str.contains("TM", case=False, na=False)]
for _, r in tm_rows.iterrows():
    ax.axvspan(r[start_col], r[end_col], color="lightgray", alpha=0.4)
```
Highlights transmembrane regions in light grey.

```bash 
pdz_rows = domains_df[domains_df[name_col].str.contains("PDZ", case=False, na=False)]
for _, r in pdz_rows.iterrows():
    ax.axvspan(r[start_col], r[end_col], color="orange", alpha=0.5)
```
Marks the PDZ-binding region in orange.

```bash 
ax.plot(cons_df[pos_col], cons_df[score_col], color="steelblue", lw=1, alpha=0.6)
ax.plot(cons_df[pos_col], cons_df["smooth"], color="navy", lw=2)
```
Thin blue line = raw scores.
Thick navy line = smoothed trend.

```bash 
ax.set_xlabel("Residue position")
ax.set_ylabel("Conservation score")
ax.set_title("VANGL1 conservation landscape")

ax.set_ylim(0.0, 1.05)
ax.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()
```
Adds labels/title, sets y-axis from 0–1, and formats layout.

```bash 
outfile = OUTPUT_DIR / "Figure1_conservation.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)
```
Saves the figure into the `figures` folder

### Large code block for effiency : 

```bash 
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

cons_file = DATA_DIR / "VANGL1_conservation.tsv"
domains_file = DATA_DIR / "VANGL1_domains.tsv"

cons_df = pd.read_csv(cons_file, sep="\t")
pos_col = cons_df.columns[0]
score_col = cons_df.columns[1]

cons_df = cons_df.sort_values(pos_col)
cons_df["smooth"] = (
    cons_df[score_col]
    .rolling(window=7, center=True, min_periods=1)
    .mean()
)

domains_df = pd.read_csv(domains_file, sep="\t")
start_col = "start"
end_col = "end"
name_col = "name"

plt.style.use("default")
fig, ax = plt.subplots(figsize=(10,4))

tm_rows = domains_df[domains_df[name_col].str.contains("TM", case=False, na=False)]
for _, r in tm_rows.iterrows():
    ax.axvspan(r[start_col], r[end_col], color="lightgray", alpha=0.4)

pdz_rows = domains_df[domains_df[name_col].str.contains("PDZ", case=False, na=False)]
for _, r in pdz_rows.iterrows():
    ax.axvspan(r[start_col], r[end_col], color="orange", alpha=0.5)

ax.plot(cons_df[pos_col], cons_df[score_col], color="steelblue", lw=1, alpha=0.6)
ax.plot(cons_df[pos_col], cons_df["smooth"], color="navy", lw=2)

ax.set_xlabel("Residue position")
ax.set_ylabel("Conservation score")
ax.set_title("VANGL1 conservation landscape")

ax.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()

outfile = OUTPUT_DIR / "Figure1_conservation.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)

```
Then Press CTRL + D 

Now, in the terminal type:

```bash 
python fig1.py
```
You've created the first figure !!!

### 3. Figure 2
(Once again a larger code block is at the end) 

### 3.1 Create the File

```bash
cd ~/submission-2-repo-vangl1/scripts
cat > fig2.py
```

### 3.2 Create the Script

For the sake of length I will only explain the most import code chunks and the full code block that needs to be pasted will be below

 ```bash
import numpy as np
```
adds jitter for scatter points.

```bash 
groups = ["NTD", "benign"]
data = [merged.loc[merged["group"] == g, "score"] for g in groups]
```
Splits conservation scores into two lists: NTD-associated vs benign variants.

```bash 
box = ax.boxplot(
    data,
    positions=[1, 2],
    widths=0.5,
    labels=["NTD variants", "Benign variants"],
    patch_artist=True,
    showfliers=False
)
```
Boxplots for the two groups, with no outlier points shown.

```bash
for b in box["boxes"]:
    b.set(facecolor="white", edgecolor="black", linewidth=1.2)
...
```
Makes the boxplot clean, white and sharp edged

```bash for i, g in enumerate(groups, start=1):
    y = merged.loc[merged["group"] == g, "score"]
    x = np.random.normal(i, 0.04, size=len(y))
    ax.scatter(x, y, color="gray", alpha=0.6, s=18, zorder=3)
```
Adds lightly jittered gray dots to show individual variants at each position.

```bash 
outfile = OUTPUT_DIR / "Figure2_variant_constraint_boxplot.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)
```
Writes the figure to the `figure` folder

### Large Code Block 

```bash 
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

cons_file = DATA_DIR / "VANGL1_conservation.tsv"
vars_file = DATA_DIR / "VANGL1_variants.tsv"

cons_df = pd.read_csv(cons_file, sep="\t")
pos_col = cons_df.columns[0]
score_col = cons_df.columns[1]

cons_df = cons_df.rename(columns={pos_col: "pos", score_col: "score"})
cons_df["pos"] = cons_df["pos"].astype(int)

vars_df = pd.read_csv(vars_file, sep="\t")
vars_df["pos"] = vars_df["pos"].astype(int)

merged = pd.merge(vars_df, cons_df, on="pos", how="left")
merged = merged.dropna(subset=["score"])

groups = ["NTD", "benign"]
data = [merged.loc[merged["group"] == g, "score"] for g in groups]

plt.style.use("default")
fig, ax = plt.subplots(figsize=(5, 4))

box = ax.boxplot(
    data,
    positions=[1, 2],
    widths=0.5,
    labels=["NTD variants", "Benign variants"],
    patch_artist=True,
    showfliers=False
)

for b in box["boxes"]:
    b.set(facecolor="white", edgecolor="black", linewidth=1.2)
for whisker in box["whiskers"]:
    whisker.set(color="black", linewidth=1.0)
for cap in box["caps"]:
    cap.set(color="black", linewidth=1.0)
for median in box["medians"]:
    median.set(color="black", linewidth=1.4)

for i, g in enumerate(groups, start=1):
    y = merged.loc[merged["group"] == g, "score"]
    x = np.random.normal(i, 0.04, size=len(y))
    ax.scatter(x, y, color="gray", alpha=0.6, s=18, zorder=3)

ax.set_ylabel("Constraint score")
ax.set_title("Constraint at VANGL1 variant sites")
ax.set_xlim(0.5, 2.5)
ax.grid(axis="y", linestyle="--", alpha=0.4)
plt.tight_layout()

outfile = OUTPUT_DIR / "Figure2_variant_constraint_boxplot.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)
```
Then press CTRL + D 

Now run it in the terminal: 

```bash 
python fig2.py
```
You've created figure 2!

### Figure 3 

### 4.1 Create the file
```bash 
cd ~/submission-2-repo-vangl1/scripts
cat > fig3.py
```
Once again I will only explain the important chunks, large code block is at the end. 
### 4.2 Scripts 

```bash 
import math
```
`math` is used for odds ratio calculations.

```bash 
cons_df = pd.read_csv(cons_file, sep="\t")
pos_col = cons_df.columns[0]
score_col = cons_df.columns[1]
cons_df = cons_df.rename(columns={pos_col: "pos", score_col: "score"})
cons_df["pos"] = cons_df["pos"].astype(int)
```
Reads conservation scores, standardizes column names to pos and score, and makes positions integers.

```bash 
vars_df = pd.read_csv(vars_file, sep="\t")
vars_df["pos"] = vars_df["pos"].astype(int)
```
Loads VANGL1 variants and makes sure variant positions are integers as well. 

```bash 
domains_df = pd.read_csv(domains_file, sep="\t")

sites = vars_df[["pos", "group"]].drop_duplicates()
```
Reads TM/PDZ annotations and keeps unique (position, group) pairs for variants.

```bash 
q90 = cons_df["score"].quantile(0.90)
q95 = cons_df["score"].quantile(0.95)
q80 = cons_df["score"].quantile(0.80)

top10_pos = set(cons_df.loc[cons_df["score"] >= q90, "pos"])
top5_pos = set(cons_df.loc[cons_df["score"] >= q95, "pos"])
top20_pos = set(cons_df.loc[cons_df["score"] >= q80, "pos"])
```
Defines which positions fall into the top 5%, 10%, and 20% most constrained residues.

```bash 
pdz_rows = domains_df[domains_df["name"].str.contains("PDZ", case=False, na=False)]
pdz_pos = set()
for _, r in pdz_rows.iterrows():
    pdz_pos.update(range(int(r["start"]), int(r["end"]) + 1))
```
Finds rows labeled PDZ and builds a set of all residue positions inside the PDZ region.

```bash categories = [
    ("Top 10% constraint", top10_pos),
    ("Top 5% constraint", top5_pos),
    ("Top 20% constraint", top20_pos),
    ("PDZ motif", pdz_pos),
]
```
Each category has a label and a set of positions where to test NTD vs benign enrichment.

```bash 
def compute_or_ci(cat_pos):
    a = ((sites["group"] == "NTD") & (sites["pos"].isin(cat_pos))).sum()
    b = ((sites["group"] == "NTD") & (~sites["pos"].isin(cat_pos))).sum()
    c = ((sites["group"] == "benign") & (sites["pos"].isin(cat_pos))).sum()
    d = ((sites["group"] == "benign") & (~sites["pos"].isin(cat_pos))).sum()

    a1, b1, c1, d1 = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    or_val = (a1 * d1) / (b1 * c1)
    log_or = math.log(or_val)
    se = math.sqrt(1/a1 + 1/b1 + 1/c1 + 1/d1)
    ci_low = math.exp(log_or - 1.96 * se)
    ci_high = math.exp(log_or + 1.96 * se)
    return or_val, ci_low, ci_high, (a, b, c, d)
```
Builds a 2×2 NTD vs benign table for a category, computes odds ratio, and adds a 95% confidence interval using log(OR) ± 1.96·SE.

```bash 
results = []
for name, pos_set in categories:
    or_val, ci_low, ci_high, counts = compute_or_ci(pos_set)
    results.append((name, or_val, ci_low, ci_high, counts))
```
Runs compute_or_ci for each category and stores the label, OR, CI, and raw counts.

```bash 
print("\nCategory\tOR\tCI_low\tCI_high\ta\tb\tc\td")
for name, or_v, lo, hi, (a, b, c, d) in results:
    print(f"{name}\t{or_v:.2f}\t{lo:.2f}\t{hi:.2f}\t{a}\t{b}\t{c}\t{d}")
```
Prints a mini table of odds ratios, confidence intervals, and raw counts.

### Large Code Block 

```bash 
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import math

DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

cons_file = DATA_DIR / "VANGL1_conservation.tsv"
vars_file = DATA_DIR / "VANGL1_variants.tsv"
domains_file = DATA_DIR / "VANGL1_domains.tsv"

cons_df = pd.read_csv(cons_file, sep="\t")
pos_col = cons_df.columns[0]
score_col = cons_df.columns[1]
cons_df = cons_df.rename(columns={pos_col: "pos", score_col: "score"})
cons_df["pos"] = cons_df["pos"].astype(int)

vars_df = pd.read_csv(vars_file, sep="\t")
vars_df["pos"] = vars_df["pos"].astype(int)

domains_df = pd.read_csv(domains_file, sep="\t")

sites = vars_df[["pos", "group"]].drop_duplicates()

q90 = cons_df["score"].quantile(0.90)
q95 = cons_df["score"].quantile(0.95)
q80 = cons_df["score"].quantile(0.80)

top10_pos = set(cons_df.loc[cons_df["score"] >= q90, "pos"])
top5_pos = set(cons_df.loc[cons_df["score"] >= q95, "pos"])
top20_pos = set(cons_df.loc[cons_df["score"] >= q80, "pos"])

pdz_rows = domains_df[domains_df["name"].str.contains("PDZ", case=False, na=False)]
pdz_pos = set()
for _, r in pdz_rows.iterrows():
    pdz_pos.update(range(int(r["start"]), int(r["end"]) + 1))

categories = [
    ("Top 10% constraint", top10_pos),
    ("Top 5% constraint", top5_pos),
    ("Top 20% constraint", top20_pos),
    ("PDZ motif", pdz_pos),
]

def compute_or_ci(cat_pos):
    a = ((sites["group"] == "NTD") & (sites["pos"].isin(cat_pos))).sum()
    b = ((sites["group"] == "NTD") & (~sites["pos"].isin(cat_pos))).sum()
    c = ((sites["group"] == "benign") & (sites["pos"].isin(cat_pos))).sum()
    d = ((sites["group"] == "benign") & (~sites["pos"].isin(cat_pos))).sum()

    a1, b1, c1, d1 = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    or_val = (a1 * d1) / (b1 * c1)
    log_or = math.log(or_val)
    se = math.sqrt(1/a1 + 1/b1 + 1/c1 + 1/d1)
    ci_low = math.exp(log_or - 1.96 * se)
    ci_high = math.exp(log_or + 1.96 * se)
    return or_val, ci_low, ci_high, (a, b, c, d)

results = []
for name, pos_set in categories:
    or_val, ci_low, ci_high, counts = compute_or_ci(pos_set)
    results.append((name, or_val, ci_low, ci_high, counts))

plt.style.use("default")
fig, ax = plt.subplots(figsize=(6,4))

y = np.arange(len(results))[::-1]
or_vals = [r[1] for r in results]
ci_lows = [r[2] for r in results]
ci_highs = [r[3] for r in results]

xmin = min(ci_lows) * 0.6
xmax = max(ci_highs) * 1.6

colors = ["tab:blue", "tab:blue", "tab:blue", "tab:orange"]

for yi, or_v, lo, hi, col in zip(y, or_vals, ci_lows, ci_highs, colors):
    ax.plot([lo, hi], [yi, yi], color=col, linewidth=1.8)
    ax.scatter(or_v, yi, color=col, s=40, zorder=3)

ax.axvline(1.0, color="gray", linestyle="--", linewidth=1)

ax.set_yticks(y)
ax.set_yticklabels([r[0] for r in results])

ax.set_xscale("log")
ax.set_xlim(xmin, xmax)
ax.set_xlabel("Odds ratio (NTD vs benign)")
ax.set_title("NTD variant enrichment in key categories")

ax.grid(axis="x", linestyle="--", alpha=0.3)
plt.tight_layout(rect=[0.18, 0.08, 0.98, 0.95])

outfile = OUTPUT_DIR / "Figure3_enrichment_forest.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)

print("\nCategory\tOR\tCI_low\tCI_high\ta\tb\tc\td")
for name, or_v, lo, hi, (a, b, c, d) in results:
    print(f"{name}\t{or_v:.2f}\t{lo:.2f}\t{hi:.2f}\t{a}\t{b}\t{c}\t{d}")
```
Then Press CTRL + D 

Now run it: 

```bash 
python fig3.py
```
You've created Figure 3! 

### 5 Figure 4 

### 5.1 Create the file 
```bash 
cd ~/submission-2-repo-vangl1/scripts
cat > fig4.py
```
### 5.2 Scripts

```bash 
from collections import Counter
```
Counter counts amino acids at each position.

```bash 
seqs = []
current = []
with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        seqs.append("".join(current))
```
Reads the FASTA file line-by-line, collects sequence lines for each header, and stores all aligned sequences in seqs.

```bash 
window = 40
start = max(0, L - window)
end = L

sub_seqs = [s[start:end] for s in seqs]
sub_L = end - start
```
Slices out the last ~40 alignment columns to focus on the C-terminal region where the PDZ motif lives.


```bash 
consensus = []
frac_cons = []
```
Consensus will store the consensus amino acid at each column; frac_cons stores how frequent that residue is.

```bash 
for i in range(sub_L):
    col = [s[i] for s in sub_seqs]
    nongap = [aa for aa in col if aa not in ("-", ".")]
    if not nongap:
        consensus.append("-")
        frac_cons.append(0.0)
        continue
    counts = Counter(nongap)
    aa, count = counts.most_common(1)[0]
    consensus.append(aa)
    frac_cons.append(count / len(nongap))
```
For each column in the C-terminal window, it finds the most common non-gap amino acid and computes the fraction of sequences that carry that residue.

```bash 
pdz_len = 4
pdz_start_idx = max(0, sub_L - pdz_len)
pdz_mask = np.zeros(sub_L, dtype=bool)
pdz_mask[pdz_start_idx:] = True
```
Marks the last 4 positions of the 40-aa window as the PDZ-binding motif. (I extracted the last ~40 alignment columns because the PDZ-binding motif occurs at the extreme C-terminus of VANGL1. A 40-residue window allows us to zoom in on this functional region while still providing surrounding sequence context for comparison.) 

```bash
ax.set_xticks(x)
ax.set_xticklabels(consensus, fontsize=8)
ax.set_xlim(0.5, sub_L + 0.5)
ax.set_ylim(0, 1.05)
```
Labels each position with its consensus amino acid and keeps the y-axis between 0 and just above 1.0.

### 5.3 Large Code Block 
``` bash 
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from collections import Counter

DATA_DIR = Path("data")
OUTPUT_DIR = Path("figures")
OUTPUT_DIR.mkdir(exist_ok=True)

fasta_file = DATA_DIR / "VANGL1_aligned.fasta"

seqs = []
current = []
with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        seqs.append("".join(current))

if not seqs:
    raise ValueError("No sequences found in FASTA")

L = len(seqs[0])
for s in seqs:
    if len(s) != L:
        raise ValueError("Sequences are not all the same length, alignment may be broken")

window = 40
start = max(0, L - window)
end = L

sub_seqs = [s[start:end] for s in seqs]
sub_L = end - start

consensus = []
frac_cons = []

for i in range(sub_L):
    col = [s[i] for s in sub_seqs]
    nongap = [aa for aa in col if aa not in ("-", ".")]
    if not nongap:
        consensus.append("-")
        frac_cons.append(0.0)
        continue
    counts = Counter(nongap)
    aa, count = counts.most_common(1)[0]
    consensus.append(aa)
    frac_cons.append(count / len(nongap))

consensus = np.array(consensus)
frac_cons = np.array(frac_cons)

pdz_len = 4
pdz_start_idx = max(0, sub_L - pdz_len)
pdz_mask = np.zeros(sub_L, dtype=bool)
pdz_mask[pdz_start_idx:] = True

plt.style.use("default")
fig, ax = plt.subplots(figsize=(8,3))

x = np.arange(1, sub_L + 1)

ax.bar(x[~pdz_mask], frac_cons[~pdz_mask],
       color="tab:blue", alpha=0.8, label="C-terminal region")
ax.bar(x[pdz_mask], frac_cons[pdz_mask],
       color="tab:orange", alpha=0.9, label="PDZ-binding motif")

ax.set_xticks(x)
ax.set_xticklabels(consensus, fontsize=8)
ax.set_xlim(0.5, sub_L + 0.5)
ax.set_ylim(0, 1.05)

ax.set_ylabel("Fraction consensus")
ax.set_xlabel("Aligned C-terminal positions (last ~40 aa)")
ax.set_title("C-terminal VANGL1 alignment highlighting PDZ motif")
ax.legend(frameon=True, fontsize=8)

ax.grid(axis="y", linestyle="--", alpha=0.4)
plt.tight_layout()

outfile = OUTPUT_DIR / "Figure4_C_terminal_PDZ_region.png"
plt.savefig(outfile, dpi=300)
print("Saved:", outfile)
```

Then click CTRL + D 

Now run it 
```bash 
python fig4.py 
```
You've Created Figure 4 !! 

### 6. Figure 5 

### 6.1 Create the file 

```bash
cd ~/submission-2-repo-vangl1/scripts
cat > fig5.py
```

### 6.2 Create the Script 

```bash
records = []
with open(pdb_file) as f:
    for line in f:
        line = line.strip("\n")
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            resnum = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        records.append((resnum, x, y, z))
```
Reads the PDB file line by line, keeps only ATOM records for C-alpha (CA) atoms, and stores residue index plus 3D coordinates for each residue.

```bash
resnums = np.array([r[0] for r in records])
xs = np.array([r[1] for r in records])
ys = np.array([r[2] for r in records])
zs = np.array([r[3] for r in records])
```
Separates residue numbers and coordinates into four numpy arrays so they can be easily plotted and normalized.

```bash
colors = (resnums - resnums.min()) / (resnums.max() - resnums.min())
```
Creates a color gradient 

```bash
plt.style.use("default")
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection="3d")
```
Creates a 3D axes object where the VANGL1 backbone will be drawn.

```bash
ax.plot(xs, ys, zs, linewidth=1.5, alpha=0.7)
sc = ax.scatter(xs, ys, zs, c=colors, cmap="viridis", s=10)
```
Draws a line through the Cα atoms to trace the backbone and overlays small points colored from N-terminus to C-terminus using the viridis colormap.`

### Large Code Block 

```bash
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
import numpy as np

DATA_DIR = Path("data")
OUT = Path("figures")
OUT.mkdir(exist_ok=True)

pdb_file = DATA_DIR / "VANGL1_AF2.pdb"

records = []
with open(pdb_file) as f:
    for line in f:
        line = line.strip("\n")
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            resnum = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        records.append((resnum, x, y, z))

if not records:
    raise RuntimeError("No C-alpha atoms found in PDB")

resnums = np.array([r[0] for r in records])
xs = np.array([r[1] for r in records])
ys = np.array([r[2] for r in records])
zs = np.array([r[3] for r in records])

colors = (resnums - resnums.min()) / (resnums.max() - resnums.min())

plt.style.use("default")
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111, projection="3d")

ax.plot(xs, ys, zs, linewidth=1.5, alpha=0.7)
sc = ax.scatter(xs, ys, zs, c=colors, cmap="viridis", s=10)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("VANGL1 AlphaFold structure (Cα trace, colored N→C)")

cb = fig.colorbar(sc, ax=ax, shrink=0.6)
cb.set_label("Residue index (N-terminus → C-terminus)")

plt.tight_layout()
out_file = OUT / "Figure5_VANGL1_structure.png"
plt.savefig(out_file, dpi=300)
print("Saved Figure 5 to:", out_file)
```
Then press CRTL + D

Now Run it 
```bash 
python fig5.py
```

You've Created Figure 5!

### 7. Figure 6

### 7.1 Create the File 

```bash
cd ~/submission-2-repo-vangl1/scripts
cat > fig6_simple_tree.R
```
We're calling this by a different name to remind you that we're going to use R to run this and it will be a pdf rather than png. 

### 7.2 Create the Script 

```bash 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript fig6_simple_tree.R input_treefile output_pdf")
}
```
This grabs whatever you pass after Rscript fig6_simple_tree.R ... and makes sure there are at least two arguments.

```bash
in_tree <- args[1]
out_pdf <- args[2]
```
`in_tree` is the path to the Newick tree file

```bash
suppressPackageStartupMessages(library(ape))
```
Loads the ape package, which has functions for phylogenetic trees

```bash
tree <- read.tree(in_tree)
tree <- ladderize(tree)
```
Reads a phylogenetic tree from the Newick file into tree, then “ladderizes” it so the branches are ordered in a more visually tidy, left-to-right way.

### 7.3 Large Code Block 
```bash
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript fig6_simple_tree.R input_treefile output_pdf")
}

in_tree <- args[1]
out_pdf <- args[2]

suppressPackageStartupMessages(library(ape))

tree <- read.tree(in_tree)
tree <- ladderize(tree)

pdf(out_pdf, width = 6, height = 7)
par(mar = c(1, 1, 1, 1))

plot(tree,
     cex = 0.5,
     no.margin = TRUE)

add.scale.bar()
dev.off()

cat("Saved simple tree figure to", out_pdf, "\n")
``` 
Then press CTRL + D 

Now run it using this command line : 

```bash
Rscript fig6_simple_tree.R data/VANGL1.treefile \
    figures/Figure6_VANGL1_simple_phylogeny.pdf
```
You've Created figure 6!

### 8. Table 1

### 8.1 Creating the file 

For this we will use `nano`

```bash 
cd ~/submission-2-repo-vangl1/tables
nano VANGL1_domains.tsv
```

PASTE THIS EXACTLY 

```bash 
Region	Residue_start	Residue_end	Domain_length
TM1	118	138	21
TM2	152	172	21
TM3	183	203	21
TM4	223	243	21
C_terminal_tail	424	524	101
PDZ_binding_motif	521	524	4
```
Then press : 
CTRL + O
then: 
Enter
then: 
CTRL + X

To get a clean image of the table download the tsv by left click -> download 

You've Created table 1 !!


9 Push any new files to the remote Git repository.
Add all of your results to the repository, commit them, and push them to the remote repository at GitHub.com

First, ensure you are in the main directory:

```bash
cd submission-2-repo-vangl1
```
 
Now, commit your changes. Make sure each line works before you continue with the next line. 

```bash
git status
git add -A
git commit -a -m "Final Push of All Figures."
``` 
Did all of the above work? Make sure! If so, then run the following command:

```bash
git push
```

### The end. 
