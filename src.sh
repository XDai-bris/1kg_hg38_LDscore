#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

################################
# STAGE TOGGLES (0=skip, 1=run)
DO_HM3=1         # HM3 filter + split per chr from merged POP PLINKs
DO_BUILD_MAPS=1  # Build 3-col headered maps from Beagle/PLINK GRCh38 maps
DO_PAINT_CM=1    # Add cM positions to per-chr HM3 files
DO_LDSC=1        # Run LDSC per chr (uses *_cm if USE_CM=1 else *_hm3)
USE_CM=1         # If DO_LDSC=1: 1 => --ld-wind-cm 1 (needs *_cm); 0 => --ld-wind-kb 1000
DO_ADD_MAF=0     # <<< NEW: add MAF into existing ldscore files by joining PLINK --freq
DO_MERGE=1       # Merge per-chr outputs into POP.l2.ldscore.gz + .M + .M_5_50
FORCE_MERGE=0    # Overwrite merged outputs if they already exist
################################

# --- Environment (use the conda env that ships with ldsc-2.0.1) ---
module --force purge
module load languages/python/ldsc-2.0.1
module load plink/1.9-beta6.27-openblas-ibxp 2>/dev/null || true
unset PYTHONPATH
export PYTHONNOUSERSITE=1

# --- Paths (edit if needed) ---
MERGED_DIR="/user/home/xd14188/repo/1kg_hg38/merged_genome"
RES_DIR="/user/home/xd14188/repo/1kg_hg38_LDscore/resources"
OUTROOT="ldsc_out_NoMAF"   
HM3_SNPLIST="${RES_DIR}/w_hm3.snplist"            # input 3-col file with header
HM3_IDS="${RES_DIR}/w_hm3.ids"                    # derived 1-col rsID list

# genetic maps (source → derived)
MAP_SRC_DIR="${RES_DIR}/genetic_maps_hg38"        # contains plink.chr*.GRCh38.map
MAP_2COL_DIR="${MAP_SRC_DIR}/clean_2col"          # intermediate (BP cM)
MAP_3COL_DIR="${MAP_SRC_DIR}/three_col"           # final 3-col with header (position, rate, cM)

# outputs
HM3_DIR="${RES_DIR}/hm3_filtered"
HM3_CM_DIR="${RES_DIR}/hm3_filtered_cm"

# Patched ldsc.py to avoid the Py3 .M write bug
if [[ -f "${HOME}/ldsc_patched/ldsc.py" ]]; then
  LDSC_RUN=(python "${HOME}/ldsc_patched/ldsc.py")
else
  echo "[WARN] ~/ldsc_patched/ldsc.py not found; using env ldsc.py (Py3 build may crash on some clusters)."
  LDSC_RUN=(ldsc.py)
fi

# Cohorts / chromosomes
POPS=(AFR AMR EAS EUR SAS)
#CHRS=(X)                     # X only
CHRS=({1..22})                # autosomes (add X by changing to: CHRS=({1..22} X))

mkdir -p "${OUTROOT}" logs "${HM3_DIR}" "${HM3_CM_DIR}"

echo "[env] python = $(which python)"
echo "[env] ldsc   = ${LDSC_RUN[*]}"
echo "[env] OUTROOT= ${OUTROOT}"

trap 'echo "[ERR] line ${LINENO} exit=${?}" >&2' ERR
log(){ echo "[$(date +%H:%M:%S)] $*"; }

# ---------- 0) prep HM3 id list ----------
prep_hm3_ids() {
  if [[ ! -f "${HM3_IDS}" ]]; then
    [[ -f "${HM3_SNPLIST}" ]] || { echo "[ERROR] Missing ${HM3_SNPLIST}"; exit 1; }
    log "Create HM3 1-col IDs -> ${HM3_IDS}"
    tail -n +2 "${HM3_SNPLIST}" | cut -f1 > "${HM3_IDS}"
    log "HM3 count: $(wc -l < "${HM3_IDS}")"
  else
    log "HM3 IDs exist: ${HM3_IDS} (n=$(wc -l < "${HM3_IDS}"))"
  fi
}

# ---------- 1) HM3 filter & split per chr ----------
do_hm3_filter_split() {
  prep_hm3_ids
  for POP in "${POPS[@]}"; do
    IN="${MERGED_DIR}/${POP}_genome"
    OUT="${HM3_DIR}/${POP}_hm3"
    if [[ ! -f "${IN}.bed" ]]; then
      echo "[WARN] Missing ${IN}.bed — skip ${POP}"; continue
    fi
    if [[ ! -f "${OUT}.bed" ]]; then
      log "[HM3] filter ${POP}"
      plink --bfile "${IN}" --extract "${HM3_IDS}" --make-bed --out "${OUT}"
    else
      log "[HM3] ${POP}: filtered exists"
    fi
    for CHR in "${CHRS[@]}"; do
      CHR_OUT="${HM3_DIR}/${POP}_hm3_chr${CHR}"
      if [[ ! -f "${CHR_OUT}.bed" ]]; then
        log "[SPLIT] ${POP} chr${CHR}"
        plink --bfile "${OUT}" --chr "${CHR}" --make-bed --out "${CHR_OUT}"
      fi
    done
  done
}

# ---------- 2) Build 3-col headered maps from Beagle/PLINK GRCh38 ----------
build_maps_3col() {
  mkdir -p "${MAP_2COL_DIR}" "${MAP_3COL_DIR}"
  for CHR in "${CHRS[@]}"; do
    src="${MAP_SRC_DIR}/plink.chr${CHR}.GRCh38.map"
    [[ -f "${src}" ]] || { echo "[WARN] missing ${src} — skip chr${CHR}"; continue; }
    out2="${MAP_2COL_DIR}/chr${CHR}.bp_cm.map"
    awk 'NR==1 && ($0 ~ /position|Genetic_Map/) {next}
         NF>=4 {print $4, $3; next}   # CHR ID cM BP
         NF==3 {print $1, $3; next}   # BP rate cM
         NF==2 {print $1, $2; next}' "${src}" > "${out2}"
    out3="${MAP_3COL_DIR}/chr${CHR}.bp_rate_cm.map"
    { echo -e "position\tCOMBINED_rate\tGenetic_Map"; awk '{printf "%s\t0\t%s\n",$1,$2}' "${out2}"; } > "${out3}"
    log "[MAP] chr${CHR} -> $(wc -l < "${out3}") rows"
  done
}

# ---------- 3) Paint cM into per-chr HM3 files ----------
paint_cm() {
  mkdir -p "${HM3_CM_DIR}"
  for POP in "${POPS[@]}"; do
    for CHR in "${CHRS[@]}"; do
      IN="${HM3_DIR}/${POP}_hm3_chr${CHR}"
      OUT="${HM3_CM_DIR}/${POP}_hm3_chr${CHR}"
      [[ -f "${IN}.bed" ]] || { echo "[WARN] ${IN}.bed missing — skip"; continue; }
      if [[ -f "${OUT}.bed" ]]; then log "[cM-map] skip ${POP} chr${CHR} (exists)"; continue; fi
      log "[cM-map] ${POP} chr${CHR}"
      plink --bfile "${IN}" --cm-map "${MAP_3COL_DIR}/chr@.bp_rate_cm.map" --make-bed --out "${OUT}"
    done
  done
}

# ---------- 4) Run LDSC per chr ----------
run_ldsc_per_chr() {
  for POP in "${POPS[@]}"; do
    OUTDIR="${OUTROOT}/${POP}_w_hm3"
    mkdir -p "${OUTDIR}"

    for CHR in "${CHRS[@]}"; do
      if [[ "${USE_CM}" -eq 1 ]]; then
        BFILE="${HM3_CM_DIR}/${POP}_hm3_chr${CHR}"
        WINDOW_FLAG=(--ld-wind-cm 1)
      else
        BFILE="${HM3_DIR}/${POP}_hm3_chr${CHR}"
        WINDOW_FLAG=(--ld-wind-kb 1000)
      fi

      OUTPRE="${OUTDIR}/chr${CHR}"

      if [[ ! -f "${BFILE}.bed" ]]; then
        echo "[WARN] Missing ${BFILE}.bed — skip ${POP} chr${CHR}"
        continue
      fi
      if [[ -f "${OUTPRE}.l2.ldscore.gz" ]]; then
        echo "[SKIP] ${OUTPRE}.l2.ldscore.gz"
        continue
      fi

      log "[LDSC] ${POP} chr${CHR} ${WINDOW_FLAG[*]}"
      "${LDSC_RUN[@]}" \
        --bfile "${BFILE}" \
        --l2 "${WINDOW_FLAG[@]}" \
        --out "${OUTPRE}"
    done
  done
}
# ---------- 5) Add MAF into each chr*.l2.ldscore.gz by joining PLINK --freq ----------
add_maf_to_ldscore() {
  module load plink/1.9-beta6.27-openblas-ibxp 2>/dev/null || true

  for POP in "${POPS[@]}"; do
    for CHR in "${CHRS[@]}"; do
      # BFILE must match what you used for LDSC l2
      if [[ "${USE_CM}" -eq 1 ]]; then
        BFILE="${HM3_CM_DIR}/${POP}_hm3_chr${CHR}"
      else
        BFILE="${HM3_DIR}/${POP}_hm3_chr${CHR}"
      fi
      LDSC_FILE="${OUTROOT}/${POP}_w_hm3/chr${CHR}.l2.ldscore.gz"

      [[ -f "${BFILE}.bed" ]] || { echo "[WARN] no BFILE ${BFILE}.bed"; continue; }
      [[ -f "${LDSC_FILE}"  ]] || { echo "[WARN] no LDscore ${LDSC_FILE}"; continue; }

      # If LD score already has MAF, skip
      # (disable pipefail in the subshell to avoid SIGPIPE from head)
      if ( set +o pipefail; zcat "${LDSC_FILE}" | head -1 | grep -q 'MAF' ); then
        echo "[SKIP] ${POP} chr${CHR}: MAF already present"
        continue
      fi

      # Make (or reuse) PLINK frequency file
      FREQ_DIR="${RES_DIR}/freq_${POP}"
      mkdir -p "${FREQ_DIR}"
      FRQ="${FREQ_DIR}/${POP}_chr${CHR}.frq"
      if [[ ! -f "${FRQ}" ]]; then
        echo "[$(date +%H:%M:%S)] [FREQ] ${POP} chr${CHR}"
        plink --bfile "${BFILE}" --freq --out "${FREQ_DIR}/${POP}_chr${CHR}"
      fi

      # Unzip ldscore to a temp plain-text file
      TMP_TXT="${OUTROOT}/${POP}_w_hm3/chr${CHR}.l2.ldscore.txt"
      zcat "${LDSC_FILE}" > "${TMP_TXT}" || { echo "[WARN] cannot zcat ${LDSC_FILE}"; rm -f "${TMP_TXT}"; continue; }

      # Output gzip (with MAF appended)
      TMP_GZ="${OUTROOT}/${POP}_w_hm3/chr${CHR}.l2.ldscore.maf.gz"

      echo "[$(date +%H:%M:%S)] [MAF] join MAF into $(basename "${LDSC_FILE}")"
      LC_ALL=C awk -v OFS='\t' '
        # Pass 1: read .frq (whitespace-delimited); build SNP->MAF map (skip header)
        NR==FNR {
          if (FNR==1) next;
          # .frq default cols: CHR SNP A1 A2 MAF NCHROBS
          if ($2!="") maf[$2]=$5;
          next;
        }
        # Pass 2: ldscore table; detect header and SNP column
        FNR==1 {
          snpcol=0;
          for (i=1;i<=NF;i++) if ($i=="SNP") { snpcol=i; break; }
          if (snpcol==0) snpcol=2;     # fallback: 2nd column is SNP
          haveMAF=0;
          for (i=1;i<=NF;i++) if ($i=="MAF") haveMAF=1;
          if (haveMAF==0) print $0,"MAF"; else print $0;
          next;
        }
        {
          snp=$snpcol; m=maf[snp]; if (m=="") m="NA";
          if (haveMAF==0) print $0,m; else print $0;
        }
      ' "${FRQ}" "${TMP_TXT}" | gzip -c > "${TMP_GZ}"

      # Pipefail-safe header check on the tmp file
      header_has_maf=false
      if [[ -s "${TMP_GZ}" ]]; then
        if ( set +o pipefail; zcat "${TMP_GZ}" | head -1 | grep -q 'MAF' ); then
          header_has_maf=true
        fi
      fi

      if [[ "${header_has_maf}" == true ]]; then
        mv -f "${TMP_GZ}" "${LDSC_FILE}"
        rm -f "${TMP_TXT}"
      else
        echo "[WARN] failed to add MAF for ${POP} chr${CHR}"
        echo "       [DBG] header of original ldscore:"
        ( set +o pipefail; zcat "${LDSC_FILE}" | head -1 ) || true
        echo "       [DBG] header of joined (tmp) file:"
        ( set +o pipefail; zcat "${TMP_GZ}" 2>/dev/null | head -1 ) || echo "       (no tmp content)"
        echo "       [DBG] counts: ldscore=$( ( set +o pipefail; zcat "${LDSC_FILE}" | wc -l ) | awk "{print \$1}") tmp=$( ( set +o pipefail; zcat "${TMP_GZ}" 2>/dev/null | wc -l ) | awk "{print \$1}") frq=$(wc -l < "${FRQ}")"
        # keep tmp files for inspection
        mv -f "${TMP_TXT}" "${TMP_TXT}.keep" 2>/dev/null || true
        mv -f "${TMP_GZ}"  "${TMP_GZ}.keep"  2>/dev/null || true
      fi
    done
  done
}

# ---------- Merge per-chromosome outputs (MAF column preserved) ----------
merge_outputs() {
  shopt -s nullglob
  for POP in "${POPS[@]}"; do
    D="${OUTROOT}/${POP}_w_hm3"
    [[ -d "$D" ]] || { echo "[WARN] missing ${D}"; continue; }

    merged_ld="${D}/${POP}.l2.ldscore.gz"
    merged_M="${D}/${POP}.l2.M"
    merged_M55="${D}/${POP}.l2.M_5_50"

    if [[ "${FORCE_MERGE:-0}" -eq 0 && -f "${merged_ld}" ]]; then
      echo "[SKIP] ${POP}: merged exists (FORCE_MERGE=1 to overwrite)"
      continue
    fi

    # Collect per-chr ldscore files in numeric order
    files=()
    for C in "${CHRS[@]}"; do
      f="${D}/chr${C}.l2.ldscore.gz"
      if [[ -f "$f" ]] && gzip -t "$f" 2>/dev/null; then
        files+=("$f")
      fi
    done
    if [[ ${#files[@]} -eq 0 ]]; then
      echo "[WARN] ${POP}: no usable chr .l2.ldscore.gz files"; continue
    fi

    echo "[MERGE] ${POP} (${#files[@]} chromosomes)"

    # Merge ldscore tables (keep header once); be pipefail-safe for head
    {
      set +o pipefail
      zcat "${files[0]}" | head -n 1
      set -o pipefail
      for f in "${files[@]}"; do
        zcat "$f" | awk 'NR>1'
      done
    } | gzip -c > "${merged_ld}" || { echo "[ERROR] merge failed for ${POP}"; continue; }

    # Sum M and M_5_50 if present
    if compgen -G "${D}/chr*.l2.M" > /dev/null; then
      awk '{s+=$1} END{print s}' "${D}"/chr*.l2.M > "${merged_M}" || true
    fi
    if compgen -G "${D}/chr*.l2.M_5_50" > /dev/null; then
      awk '{s+=$1} END{print s}' "${D}"/chr*.l2.M_5_50 > "${merged_M55}" || true
    fi

    rows=$( ( set +o pipefail; zcat "${merged_ld}" | wc -l ) | awk '{print $1}')
    echo "  done ${POP}: rows=${rows}"
    [[ -f "${merged_M}"   ]] && echo "           M=$(cat "${merged_M}")"
    [[ -f "${merged_M55}" ]] && echo "      M_5_50=$(cat "${merged_M55}")"
  done
  echo "[DONE] merged outputs in ${OUTROOT}/<POP>_w_hm3/"
}

# ----------------- run selected stages -----------------
# ... define: prep_hm3_ids, do_hm3_filter_split, build_maps_3col, paint_cm,
#             run_ldsc_per_chr, add_maf_to_ldscore, merge_outputs

# ----------------- run selected stages -----------------
[[ "${DO_HM3}"        -eq 1 ]] && do_hm3_filter_split
[[ "${DO_BUILD_MAPS}" -eq 1 ]] && build_maps_3col
[[ "${DO_PAINT_CM}"   -eq 1 ]] && paint_cm
[[ "${DO_LDSC}"       -eq 1 ]] && run_ldsc_per_chr
# [[ "${DO_ADD_MAF}"    -eq 1 ]] && add_maf_to_ldscore
[[ "${DO_MERGE}"      -eq 1 ]] && merge_outputs

echo "[ALL DONE]"