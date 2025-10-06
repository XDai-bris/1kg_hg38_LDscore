# clean-ish start (slurm stays)
module --force purge
module load languages/python/ldsc-2.0.1

unset PYTHONPATH
export PYTHONNOUSERSITE=1

# sanity checks
which python
which ldsc.py
python - <<'PY'
import sys, numpy, scipy
print(sys.executable)
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
PY

# now run one LDSC job
ldsc.py \
  --bfile /user/home/xd14188/repo/1kg_hg38_LDscore/resources/hm3_filtered/AFR_hm3_chr1 \
  --l2 --ld-wind-kb 1000 \
  --out ldsc_out/AFR_w_hm3/chr1




hm3 data download
mkdir -p resources && cd resources
wget -O w_hm3.snplist.gz https://zenodo.org/records/7773502/files/w_hm3.snplist.gz
gunzip -f w_hm3.snplist.gz


# Genetic maps for hg38 from Beagle
https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/?utm_source=chatgpt.com
cd /user/home/xd14188/repo/1kg_hg38_LDscore/resources
mkdir -p genetic_maps_hg38 && cd genetic_maps_hg38

# Beagle provides PLINK-format maps for GRCh38
wget -q https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip -o plink.GRCh38.map.zip   # yields files like plink.chr1.GRCh38.map ... plink.chr22.GRCh38.map

1. run LDSC, using the new ldscore to estimate the heri r2 to see if it works
2. git and upload to the 1kg repo 