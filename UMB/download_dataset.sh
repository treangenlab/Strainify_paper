#!/usr/bin/env bash
set -euo pipefail

BIOPROJECT="PRJNA400628"
PATTERN="UMB18"   # substring to match in isolate/attributes
OUTDIR="/dodo/rl152/UMB/PRJNA400628_UMB18"

THREADS=32
PREFETCH_TRIES=5
PREFETCH_SLEEP=20

# ---- Requirements ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing '$1' in PATH" >&2; exit 1; }; }
need esearch
need efetch
need prefetch
need fasterq-dump
need python3
need gzip

mkdir -p "$OUTDIR"
cd "$OUTDIR"

mkdir -p sra fastq

retry_prefetch () {
  local srr="$1"
  local tries="${2:-$PREFETCH_TRIES}"
  local sleep_s="${3:-$PREFETCH_SLEEP}"

  # clean partial lite objects that can block retries
  rm -f "sra/${srr}.sralite"* 2>/dev/null || true

  for i in $(seq 1 "$tries"); do
    echo "    prefetch ${srr} (attempt ${i}/${tries})"
    if prefetch -O sra "$srr"; then
      return 0
    fi
    echo "    prefetch failed for ${srr}; sleeping ${sleep_s}s..." >&2
    sleep "$sleep_s"
  done
  return 1
}

is_downloaded () {
  local sample="$1"
  # consider downloaded if both mates exist either gz or plain fastq and are non-empty
  if [[ -s "fastq/${sample}_r1.fastq.gz" && -s "fastq/${sample}_r2.fastq.gz" ]]; then
    return 0
  fi
  if [[ -s "fastq/${sample}_r1.fastq" && -s "fastq/${sample}_r2.fastq" ]]; then
    return 0
  fi
  return 1
}

echo "==> Fetching SRA RunInfo for ${BIOPROJECT}"
esearch -db sra -query "${BIOPROJECT}" | efetch -format runinfo > runinfo.csv

echo "==> Filtering for: isolate/attrs contain '${PATTERN}', Strategy=WGS, Source=METAGENOMIC, Layout=PAIRED; adding collection date to mapping"
python3 - <<'PY'
import csv, re
from pathlib import Path

PATTERN = "UMB18"
pattern = re.compile(re.escape(PATTERN), re.IGNORECASE)

REQ_STRATEGY = "WGS"
REQ_SOURCE   = "METAGENOMIC"
REQ_LAYOUT   = "PAIRED"

in_csv  = Path("runinfo.csv")
out_map = Path(f"{PATTERN}_mapping.tsv")
out_runs= Path(f"{PATTERN}_SRRs.txt")

with in_csv.open(newline="", encoding="utf-8") as f:
    rdr = csv.DictReader(f)
    rows = list(rdr)

if not rows:
    raise SystemExit("No rows found in runinfo.csv (is the BioProject public / reachable?)")

fieldnames = list(rows[0].keys())

def get_col_case_insensitive(target_lc: str):
    for c in fieldnames:
        if c.strip().lower() == target_lc:
            return c
    return None

def getv(row, col):
    return (row.get(col) or "").strip()

# --- isolate matching ---
iso_col = get_col_case_insensitive("isolate")

def row_matches_pattern(row):
    if iso_col is not None:
        return bool(pattern.search(getv(row, iso_col)))
    return any(v and pattern.search(str(v)) for v in row.values())

# --- library filters ---
layout_col   = get_col_case_insensitive("librarylayout")
strategy_col = get_col_case_insensitive("librarystrategy")
source_col   = get_col_case_insensitive("librarysource")

missing = [name for name, col in [("LibraryLayout",layout_col),
                                 ("LibraryStrategy",strategy_col),
                                 ("LibrarySource",source_col)] if col is None]
if missing:
    cols_preview = ", ".join(list(fieldnames)[:80])
    raise SystemExit(
        "ERROR: Missing required RunInfo columns: " + ", ".join(missing) + "\n"
        f"Columns seen (first ~80): {cols_preview}\n"
        "Tip: open runinfo.csv and confirm the header names.\n"
    )

def passes_library_filters(row):
    layout   = getv(row, layout_col).upper()
    strategy = getv(row, strategy_col).upper()
    source   = getv(row, source_col).upper()
    return (layout == REQ_LAYOUT and strategy == REQ_STRATEGY and source == REQ_SOURCE)

# --- collection date extraction ---
collection_candidates = [
    "collection_date",
    "collection date",
    "collectiondate",
    "date_collected",
    "sample collection date",
]
field_lc_map = {c.strip().lower(): c for c in fieldnames}
collection_col = None
for cand in collection_candidates:
    if cand in field_lc_map:
        collection_col = field_lc_map[cand]
        break
if collection_col is None:
    for c in fieldnames:
        cl = c.strip().lower()
        if "collection" in cl and "date" in cl:
            collection_col = c
            break

def get_collection_date(row):
    if collection_col is None:
        return ""
    return getv(row, collection_col)

# --- apply filters ---
hits = [r for r in rows if row_matches_pattern(r) and passes_library_filters(r)]

# fallback: if isolate col exists but nothing matched, scan-all + library filters
if not hits and iso_col is not None:
    hits = []
    for r in rows:
        if passes_library_filters(r) and any(v and pattern.search(str(v)) for v in r.values()):
            hits.append(r)

if not hits:
    raise SystemExit(
        f"No runs matched {PATTERN} with Strategy={REQ_STRATEGY}, Source={REQ_SOURCE}, Layout={REQ_LAYOUT}.\n"
        "Tip: open runinfo.csv and confirm values are exactly those strings.\n"
    )

hits.sort(key=lambda r: (r.get("BioSample",""), r.get("Run","")))

def fmt(i): return f"{PATTERN}_{i:02d}"

with out_map.open("w", encoding="utf-8", newline="") as g:
    g.write("\t".join([
        "sample_name",
        "BioSample(SAMN)",
        "Run(SRR)",
        "Experiment(SRX)",
        "LibraryStrategy",
        "LibrarySource",
        "LibraryLayout",
        "CollectionDate",
        "SampleName",
    ]) + "\n")

    for i, r in enumerate(hits, start=1):
        sample = fmt(i)
        g.write("\t".join([
            sample,
            getv(r, "BioSample"),
            getv(r, "Run"),
            getv(r, "Experiment"),
            getv(r, strategy_col),
            getv(r, source_col),
            getv(r, layout_col),
            get_collection_date(r),
            getv(r, "SampleName"),
        ]) + "\n")

with out_runs.open("w", encoding="utf-8") as g:
    for r in hits:
        srr = getv(r, "Run")
        if srr:
            g.write(srr + "\n")

print(f"Matched {len(hits)} runs. Wrote {out_map} and {out_runs}.")
if collection_col:
    print(f"CollectionDate source column: {collection_col}")
else:
    print("CollectionDate: no column found in RunInfo; left blank.")
PY

echo "==> Mapping written to: UMB18_mapping.tsv"
echo "==> SRR list written to: UMB18_SRRs.txt"

echo "==> Downloading SRA + converting to FASTQ (filtered; skip already-downloaded)"
# Columns: sample_name, SAMN, SRR, SRX, strategy, source, layout, collection_date, SampleName
tail -n +2 UMB18_mapping.tsv | while IFS=$'\t' read -r SAMPLE SAMN SRR SRX STRATEGY SOURCE LAYOUT COLLECTION_DATE SAMPLENAME; do
  [[ -z "${SRR:-}" ]] && continue

  # enforce filters again (extra safety)
  if [[ "${LAYOUT^^}" != "PAIRED" || "${STRATEGY^^}" != "WGS" || "${SOURCE^^}" != "METAGENOMIC" ]]; then
    echo "---- ${SAMPLE}: strategy=${STRATEGY} source=${SOURCE} layout=${LAYOUT} -> skipping ----"
    continue
  fi

  # skip if already downloaded (both mates exist)
  if is_downloaded "$SAMPLE"; then
    echo "---- ${SAMPLE} already downloaded -> skipping ----"
    continue
  fi

  echo "---- ${SAMPLE}  SRR=${SRR}  SAMN=${SAMN}  date=${COLLECTION_DATE} ----"

  # Download SRA (resume-safe + retries)
  if [[ ! -s "sra/${SRR}/${SRR}.sra" ]]; then
    rm -rf "sra/${SRR}" 2>/dev/null || true
    retry_prefetch "$SRR"
  fi

  sra_path="sra/${SRR}/${SRR}.sra"
  if [[ ! -s "$sra_path" ]]; then
    echo "ERROR: expected SRA file not found after prefetch: $sra_path" >&2
    exit 1
  fi

  TMPDIR="$(mktemp -d)"
  trap 'rm -rf "$TMPDIR"' EXIT

  # remove partial fastqs for this sample before re-creating
  rm -f "fastq/${SAMPLE}_r1.fastq"* "fastq/${SAMPLE}_r2.fastq"* 2>/dev/null || true

  fasterq-dump -O "$TMPDIR" --split-files --threads "$THREADS" "$sra_path"

  if [[ -s "$TMPDIR/${SRR}_1.fastq" && -s "$TMPDIR/${SRR}_2.fastq" ]]; then
    mv "$TMPDIR/${SRR}_1.fastq" "fastq/${SAMPLE}_r1.fastq"
    mv "$TMPDIR/${SRR}_2.fastq" "fastq/${SAMPLE}_r2.fastq"

    # gzip outputs (comment if you want plain .fastq)
    #gzip -f "fastq/${SAMPLE}_r1.fastq" "fastq/${SAMPLE}_r2.fastq"
  else
    echo "ERROR: Expected paired FASTQs for ${SRR} but did not find both *_1.fastq and *_2.fastq" >&2
    echo "TMPDIR contents:" >&2
    ls -lah "$TMPDIR" >&2 || true
    exit 1
  fi

  rm -rf "$TMPDIR"
  trap - EXIT
done

echo "==> Done."
echo "Outputs:"
echo "  - runinfo.csv"
echo "  - UMB18_mapping.tsv   (sample_name <-> SAMN <-> SRR + collection date + library fields)"
echo "  - fastq/              (renamed FASTQs; filtered; skips downloaded)"
echo "  - sra/                (downloaded .sra files)"
