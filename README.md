# 1kg_hg38_LDscore

A fully reproducible pipeline to **build LD Score files on GRCh38/hg38** from 1000 Genomes data and validate them via **heritability (h²) estimation** with LD Score Regression (LDSC).

---

## Overview

This repository provides:

- Scripts to generate LD Score files for **AFR, AMR, EAS, EUR, SAS** (1000G) on **hg38**.
- Integration of **HapMap3 SNPs** and **Beagle GRCh38 genetic maps** for cM-based LD windows.
- Optional addition of **MAF** columns (via PLINK `--freq`) and per-population **merged** LD score outputs.
- An end-to-end **validation** using a European BMI GWAS (GRCh38).

---

## Dependencies

Example environment setup in BC4:

```bash
module load languages/python/ldsc-2.0.1
module load plink/1.9-beta6.27-openblas-ibxp
unset PYTHONPATH
export PYTHONNOUSERSITE=1
```

Required in BC4:
Use a locally patched `ldsc.py` in `~/ldsc_patched/ldsc.py` to fix Python 3 `.M` write issues — the scripts automatically detect it.

---

## Input Resources

### 1. HapMap3 SNP list

```bash
mkdir -p resources && cd resources
wget -O w_hm3.snplist.gz https://zenodo.org/records/7773502/files/w_hm3.snplist.gz
gunzip -f w_hm3.snplist.gz
```

### 2. Beagle GRCh38 genetic maps

```bash
cd /user/home/xd14188/repo/1kg_hg38_LDscore/resources
mkdir -p genetic_maps_hg38 && cd genetic_maps_hg38
wget -q https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip -o plink.GRCh38.map.zip
# -> plink.chr{1..22}.GRCh38.map
```

---

## Pipeline — LD Score Construction

Main script: **`src.sh`**

At the top, set which stages to run:

```bash
DO_HM3=1
DO_BUILD_MAPS=1
DO_PAINT_CM=1
DO_LDSC=1
USE_CM=1
DO_ADD_MAF=0
DO_MERGE=1
FORCE_MERGE=0
```

Paths inside the script:

- `MERGED_DIR` — merged PLINK per-population 1000G data  
- `RES_DIR` — where your HapMap3 & maps live [RESOUCE] 
- `OUTROOT` — output directory (e.g. `ldsc_out_NoMAF`)

Outputs per ancestry:

```
ldsc_out_NoMAF/<POP>_w_hm3/
 ├── chr1.l2.ldscore.gz
 ├── …
 ├── chr22.l2.ldscore.gz
 ├── <POP>.l2.ldscore.gz
 ├── <POP>.l2.M
 └── <POP>.l2.M_5_50
```

To include chr 1-22:
```bash
CHRS=({1..22})
```

---

## Validation: Heritability Estimation on BMI (GRCh38)
Test script: **`/test/src.sh`**
### Input GWAS
| Field | Description |
|-------|--------------|
| **Study ID** | GCST90029007 |
| **Trait** | Body Mass Index |
| **Sample size** | 532,396 (European) |
| **Genome build** | GRCh38 |
| **File** | 29892013-GCST90029007-EFO_0004340.h.tsv.gz |
| **Modified** | 2025-01-14 |

---

### LDSC Run Example

```bash
ldsc.py \
  --h2 gwas_hg38.sumstats.gz \
  --ref-ld /user/home/xd14188/repo/1kg_hg38_LDscore/ldsc_out_NoMAF/EUR_w_hm3/EUR \
  --w-ld   /user/home/xd14188/repo/1kg_hg38_LDscore/ldsc_out_NoMAF/EUR_w_hm3/EUR \
  --out ldsc_hg38_check_merged
```

---

### Example Result

```
Total Observed scale h2: 0.255 (0.0077)
Lambda GC: 2.5641
Mean Chi^2: 3.5141
Intercept: 1.1837 (0.0172)
Ratio: 0.0731 (0.0068)
```

### Interpretation

Consistent with:
- **Wainschtein et al., Science Advances (2019)** — BMI h²ₛₙₚ = 0.27
- **Neale Lab UKB LDSC** — h²ₛₙₚ ≈ 0.25  
Thus, 0.255 (±0.0077) is fully reasonable.

---

## References

- Bulik-Sullivan et al., *Nat Genet* (2015)
- Wainschtein et al., *Sci Adv* (2019)
- Beagle GRCh38 maps — [https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)
- HapMap3 SNPs — [Zenodo 7773502](https://zenodo.org/records/7773502)
- LDSC — [https://github.com/bulik/ldsc](https://github.com/bulik/ldsc)
