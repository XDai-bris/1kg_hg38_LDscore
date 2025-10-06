#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

############################################
# EDIT THESE
GWAS_GZ="/user/home/xd14188/repo/gwas2vcf/test/gwas_b38/29892013-GCST90029007-EFO_0004340.h.tsv.gz"

# Merged ldscore base (NO suffix); expects:
#   ${LD_BASE}.l2.ldscore.gz  (+ optional ${LD_BASE}.l2.M, ${LD_BASE}.l2.M_5_50)
LD_BASE="/user/home/xd14188/repo/1kg_hg38_LDscore/ldsc_out_NoMAF/EUR_w_hm3/EUR"

# If your input has per-SNP 'n', set N_FROM_COL=1; else 0 and set CONST_N.
N_FROM_COL=1
CONST_N=300000

OUTDIR="ldsc_check_out"
JOBNAME="ldsc_hg38_check_merged"
############################################

# --- Environment (ldsc 2.0.1 conda env) ---
module --force purge
module load languages/python/ldsc-2.0.1
unset PYTHONPATH
export PYTHONNOUSERSITE=1

# --- Paths to (possibly) patched scripts ---
if [[ -f "${HOME}/ldsc_patched/ldsc.py" ]]; then
  LDSC_RUN=(python "${HOME}/ldsc_patched/ldsc.py")
else
  echo "[WARN] ~/ldsc_patched/ldsc.py not found; using system ldsc.py"
  LDSC_RUN=(ldsc.py)
fi

if [[ -f "${HOME}/ldsc_patched/munge_sumstats.py" ]]; then
  MUNGE_RUN=(python "${HOME}/ldsc_patched/munge_sumstats.py")
else
  echo "[WARN] ~/ldsc_patched/munge_sumstats.py not found; using system munge_sumstats.py"
  MUNGE_RUN=(munge_sumstats.py)
fi

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "[ENV] python=$(which python)"
echo "[ENV] LDSC = ${LDSC_RUN[*]}"
echo "[ENV] MUNGE= ${MUNGE_RUN[*]}"

# ==============================================================
# 1) Build minimal TSV for munge_sumstats
# ==============================================================
echo "[STEP] Preparing slim file from ${GWAS_GZ}"

zcat "${GWAS_GZ}" | awk -F'\t' '
BEGIN{OFS="\t"}
NR==1{
  for(i=1;i<=NF;i++){h[$i]=i}
  need="hm_rsid,effect_allele,other_allele,beta,standard_error,p_value"
  split(need, arr, ",")
  for(k in arr){ if(!(arr[k] in h)) miss=miss arr[k]" " }
  if(miss!=""){ print "[ERROR] Missing columns: "miss >"/dev/stderr"; exit 2 }
  hasN = ("n" in h)?1:0
  print "SNP","A1","A2","BETA","SE","P",(hasN?"N":"N")
  next
}
{
  snp=$h["hm_rsid"]; a1=$h["effect_allele"]; a2=$h["other_allele"]
  b=$h["beta"]; se=$h["standard_error"]; p=$h["p_value"]
  if("n" in h) n=$h["n"]; else n="."
  if(snp!="" && p!="" && p!="NA" && snp!="NA")
     print snp,a1,a2,b,se,p,n
}' > gwas_for_ldsc.tsv

echo "[INFO] Wrote gwas_for_ldsc.tsv ($(wc -l < gwas_for_ldsc.tsv) lines)"

# ==============================================================
# 2) Munge to LDSC format
# ==============================================================
echo "[STEP] Munge summary statistics"

if [[ "${N_FROM_COL}" -eq 1 ]]; then
  "${MUNGE_RUN[@]}" \
    --sumstats gwas_for_ldsc.tsv \
    --snp SNP --a1 A1 --a2 A2 --p P \
    --signed-sumstats BETA,0 \
    --N-col N \
    --out gwas_hg38
else
  "${MUNGE_RUN[@]}" \
    --sumstats gwas_for_ldsc.tsv \
    --snp SNP --a1 A1 --a2 A2 --p P \
    --signed-sumstats BETA,0 \
    --N "${CONST_N}" \
    --out gwas_hg38
fi

echo "[INFO] Created gwas_hg38.sumstats.gz"
ls -lh gwas_hg38.sumstats.gz

# ==============================================================
# 3) Run LDSC h2 using your **merged** hg38 LD scores
# ==============================================================
echo "[STEP] Run LDSC heritability (merged ref-ld)"

# Sanity: merged files present?
[[ -f "${LD_BASE}.l2.ldscore.gz" ]] || { echo "[ERR] Missing ${LD_BASE}.l2.ldscore.gz"; exit 1; }
# Print merged header (for confirmation it’s 4 cols, no MAF)
echo -n "[INFO] Merged LD header: "
zcat "${LD_BASE}.l2.ldscore.gz" | head -1 || true

"${LDSC_RUN[@]}" \
  --h2 gwas_hg38.sumstats.gz \
  --ref-ld "${LD_BASE}" \
  --w-ld   "${LD_BASE}" \
  --out "${JOBNAME}"

echo
echo '--- Key results ---'
grep -E "Total Observed scale h2|Intercept|Mean Chi|lambda GC" "${JOBNAME}.log" || true