#!/usr/bin/env bash
# Download the GTDB metadata and the NCBI proteomes or assemblies the plan needs.
#   ./download_genomes.sh check        # check the tools and the network first
#   ./download_genomes.sh bootstrap [archaea|bacteria]   # the GTDB files, once
#   ./download_genomes.sh archaea [N]  # archaeal proteomes; N = prefix of the plan
#   ./download_genomes.sh archaea-ncbi # archaeal proteomes without GTDB
#   ./download_genomes.sh pilot [500]  # bacterial pilot, measures the prevalence
#   ./download_genomes.sh bacteria N   # extend the same plan up to N
#   ./download_genomes.sh missing [N] [archaea|bacteria]  # plan rows with no proteome
#   ./download_genomes.sh predict [N] [archaea|bacteria]  # nucleotides, then Prodigal
#   ./download_genomes.sh nucleotides  # assemblies for the exporter negatives
#   ./download_genomes.sh reindex      # rebuild the journal from the files present
#   ./download_genomes.sh audit [Archaea|Bacteria]   # cost of a completeness filter
#   ./download_genomes.sh sizing [Archaea|Bacteria]  # how many genomes to draw
#   ./download_genomes.sh paired [N]   # same genome: NCBI annotation vs Prodigal
#
# Archaeal proteomes go to data/gtdb/proteomes, bacterial ones to
# data/gtdb/proteomes_bac. Any subcommand can be interrupted: re-running it
# resumes where it stopped.
set -uo pipefail

# Printed on an unknown subcommand. Spelled out here rather than sliced out of
# this file's own header by line number, which breaks silently when the header
# changes.
usage() {
    cat <<'USAGE'
Usage: ./download_genomes.sh <subcommand> [args]

  check                             check the tools and the network first
  bootstrap [archaea|bacteria]      the GTDB files, once
  archaea [N]                       archaeal proteomes; N = prefix of the plan
  archaea-ncbi                      archaeal proteomes without GTDB
  pilot [500]                       bacterial pilot, measures the prevalence
  bacteria N                        extend the same plan up to N
  missing [N] [archaea|bacteria]    plan rows with no proteome
  predict [N] [archaea|bacteria]    nucleotides, then Prodigal
  nucleotides                       assemblies for the exporter negatives
  reindex                           rebuild the journal from the files present
  audit [Archaea|Bacteria]          cost of a completeness filter
  sizing [Archaea|Bacteria]         how many genomes to draw
USAGE
}
D="data/gtdb"; mkdir -p "$D" "$D/proteomes" "$D/genomes" "$D/log"
# The link to GTDB (Australia) often fails before it succeeds, so the retry
# budget is large and -c resumes an interrupted transfer.
WG=(--tries=100 --waitretry=15 --retry-connrefused --timeout=120 --read-timeout=180 -c)
BATCH="${BATCH:-50}"
# Minimum completeness. 50 % is GTDB's own floor: completeness is not used as a
# hard filter but recorded in the plan and used as a covariate in the analysis.
# A 95 % threshold removes entire archaeal phyla (see `audit`).
MIN_COMP="${MIN_COMP:-50}"

# The plan must be rebuilt if its format has changed (completeness column).
plan_stale() {
  local f="$1"
  [ ! -s "$f" ] && return 0
  head -1 "$f" | grep -q "completeness" || return 0
  return 1
}
say() { echo "$@" >&2; }      # messages to stderr, never into a substitution

# Tools.
check_tools() {
  local ok=0
  for t in datasets wget unzip gzip awk; do
    if command -v "$t" >/dev/null; then say "  [ok]      $t"
    else say "  [MISSING] $t"; ok=1; fi
  done
  return $ok
}

# Release detection.
find_release() {
  # No network probe: on a link that keeps dropping, a --spider probe costs
  # about 2.5 min per failure and settles nothing. Take the default release and
  # let the download itself, with its large retry budget, settle the matter.
  local cached="$D/log/release.txt"
  if [ -n "${GTDB_RELEASE:-}" ]; then echo "$GTDB_RELEASE" > "$cached"; fi
  if [ -s "$cached" ]; then cat "$cached"; return 0; fi
  echo 226 > "$cached"; echo 226
}

# Metadata.
get_meta() {                       # $1 = ar53|bac120; prints the path, nothing else
  local kind="$1"
  local rel; rel=$(find_release)
  local f="$D/${kind}_metadata.tsv.gz"
  # Accept the name as it comes from the GTDB site: a file fetched by hand from
  # a browser (see NETWORK.md) is called bac120_metadata_r226.tsv.gz, and is
  # recognised under that name without being renamed.
  local c
  for c in "$f" "$D/${kind}_metadata_r"*.tsv.gz; do
    if [ -s "$c" ] && gzip -t "$c" 2>/dev/null; then echo "$c"; return 0; fi
  done

  local BASE="https://data.gtdb.ecogenomic.org/releases"
  local cand rr partial
  local candidates=()
  for rr in "$rel" 220 214 207; do
    candidates+=("release${rr}/${rr}.0/${kind}_metadata_r${rr}.tsv.gz")
  done
  candidates+=("release${rel}/${rel}.0/${kind}_metadata.tsv.gz")

  for cand in "${candidates[@]}"; do
    # One partial file per URL, so wget -c resumes exactly where it stopped.
    partial="$D/.part_$(echo "$cand" | tr '/.' '__')"
    say "  trying: ${cand}"
    say "         (unstable long-distance link: up to 100 attempts, automatic resume)"
    wget "${WG[@]}" -O "$partial" "${BASE}/${cand}" 2>&1 \
      | grep -Ei "saved|failed|error|refused|timed out" | tail -3 >&2
    if [ -s "$partial" ] && gzip -t "$partial" 2>/dev/null; then
      mv "$partial" "$f"; say "  [ok] retrieved: $(du -h "$f" | cut -f1)"
      echo "$f"; return 0
    fi
    if [ -s "$partial" ]; then
      say "  incomplete transfer ($(du -h "$partial" | cut -f1) received) -- RE-RUN the command,"
      say "  wget will resume exactly where it stopped."
      return 1
    fi
  done
  say "  [FAILED] no URL responded"
  return 1
}

get_tree() {                       # GTDB tree, small, needed for the phylogeny
  local kind="$1"
  local rel; rel=$(find_release) || return 1
  local f="$D/${kind}_r${rel}.tree.gz"
  [ -s "$f" ] && return 0
  wget "${WG[@]}" -q -O "$f" \
    "https://data.gtdb.ecogenomic.org/releases/release${rel}/${rel}.0/${kind}_r${rel}.tree.gz" \
    2>/dev/null && say "  ${kind} tree retrieved ($(du -h "$f" | cut -f1)) -- will serve for the phylogenetic correction" \
    || rm -f "$f"
}

# Extraction of the species representatives.
reps() {                           # $1 = ar53|bac120; writes $D/$1_reps.tsv
  local kind="$1"
  local m; m=$(get_meta "$kind") || return 1
  zcat "$m" | awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i
            if(!("gtdb_representative" in h) || !("accession" in h) || !("gtdb_taxonomy" in h)){
              print "MISSING_COLUMNS" > "/dev/stderr"; exit 1 }
            next }
    $(h["gtdb_representative"])=="t" {
      acc=$(h["accession"]); sub(/^(RS_|GB_)/,"",acc)
      n=split($(h["gtdb_taxonomy"]),t,";"); d=""; p=""
      for(j=1;j<=n;j++){ if(t[j]~/^d__/) d=substr(t[j],4); if(t[j]~/^p__/) p=substr(t[j],4) }
      print acc"\t"d"\t"p }' > "$D/${kind}_reps.tsv"
  local n; n=$(wc -l < "$D/${kind}_reps.tsv")
  [ "$n" -eq 0 ] && { say "  [FAILED] 0 representative extracted from $m"; return 1; }
  echo "$n"
}

# Download.
fetch() {           # $1 = accession file; $2 = protein|genome; $3 = suffix
  # The suffix separates the domains. Bacterial proteomes must not land in the
  # same directory as the archaeal ones, or a later scan mixes the two and the
  # resume journal becomes wrong.
  local list="$1"
  local what="$2"
  local suf="${3:-}"
  local done_f="$D/log/done_${what}${suf}.txt"
  local outdir; [ "$what" = protein ] && outdir="$D/proteomes${suf}" || outdir="$D/genomes${suf}"
  mkdir -p "$outdir"
  touch "$done_f"
  sort -u "$done_f" -o "$done_f"          # keep the journal deduplicated
  local todo; todo=$(mktemp)
  # grep -v exits 1 when it selects nothing, so do not chain a || onto it:
  # the resume would then ask for everything again.
  if [ -s "$done_f" ]; then
    grep -v '^accession' "$list" | cut -f1 | grep '^GC' | sort -u | grep -vxF -f "$done_f" > "$todo"
  else
    grep -v '^accession' "$list" | cut -f1 | grep '^GC' | sort -u > "$todo"
  fi
  local n; n=$(wc -l < "$todo")
  local n0; n0=$(wc -l < "$done_f")
  say ">> $n genomes to retrieve ($what) ; $n0 already done"
  if [ "$n" -eq 0 ]; then rm -f "$todo"; return 0; fi
  rm -f "$D"/log/batch_*
  split -l "$BATCH" "$todo" "$D/log/batch_"
  local i=0 tot; tot=$(ls "$D"/log/batch_* | wc -l)
  for b in "$D"/log/batch_*; do
    i=$((i+1))
    say "   batch $i/$tot ($(wc -l < "$b") genomes)  --  this run: $(( $(wc -l < "$done_f") - n0 ))/$n"
    if datasets download genome accession --inputfile "$b" --include "$what" \
         --no-progressbar --filename "$D/tmp.zip" >>"$D/log/errors.txt" 2>&1; then
      rm -rf "$D/tmp"; unzip -qo "$D/tmp.zip" -d "$D/tmp" 2>/dev/null
      local pat; [ "$what" = protein ] && pat="protein.faa" || pat="*_genomic.fna"
      local ext; [ "$what" = protein ] && ext="faa" || ext="fna"
      while IFS= read -r p; do
        [ -s "$p" ] || continue
        local g; g=$(basename "$(dirname "$p")")
        if [ -s "$outdir/${g}.${ext}.gz" ] && gzip -t "$outdir/${g}.${ext}.gz" 2>/dev/null; then
          grep -qxF "$g" "$done_f" || echo "$g" >> "$done_f"; continue; fi
        gzip -c "$p" > "$outdir/${g}.${ext}.gz" && echo "$g" >> "$done_f"
      done < <(find "$D/tmp" -name "$pat" 2>/dev/null)
      rm -rf "$D/tmp" "$D/tmp.zip"
    else
      say "   [batch $i failed -- will be resumed on the next run]"
    fi
    rm -f "$b"
  done
  rm -f "$todo"
  say ">> total retrieved: $(wc -l < "$done_f") ; volume: $(du -sh "$outdir" | cut -f1)"
}

# Main.
case "${1:-}" in
  check)
    say "=== Check of the environment ==="
    check_tools || { say "  install what is missing, then re-run"; exit 1; }
    say ""
    say "  testing the connection to GTDB (90 KB file, limited patience)..."
    rel=$(find_release)
    if wget -q --tries=8 --waitretry=10 --retry-connrefused --timeout=90 -O /dev/null \
        "https://data.gtdb.ecogenomic.org/releases/release${rel}/${rel}.0/ar53_r${rel}.tree.gz" 2>/dev/null
    then say "  [ok]      GTDB responds (r${rel})"
    else
      say "  [SLOW]    GTDB (Australia) did not respond in 8 attempts."
      say "            This is not blocking: run"
      say "                ./scripts/download_genomes.sh bootstrap"
      say "            which insists up to 100 times and resumes interrupted transfers."
      say "            Only 5.5 MB has to come from GTDB, once."
    fi
    say "  testing the connection to NCBI datasets..."
    if datasets summary genome accession GCF_000005845.2 --as-json-lines >/dev/null 2>&1; then
      say "  [ok]      NCBI datasets responds"
    else
      say "  [MISSING] NCBI datasets does not respond -- check the network/proxy"; exit 1
    fi
    say ""; say "  Everything is ready. Run:  ./scripts/download_genomes.sh archaea"
    ;;
  bootstrap)
    # Retrieve only the GTDB files, the ones that go over the long-distance
    # link. Everything else comes afterwards from NCBI, which is fast.
    DOM="${2:-archaea}"
    if [ "$DOM" = bacteria ] || [ "$DOM" = bacteries ]; then KIND=bac120; else KIND=ar53; fi
    say "=== GTDB bootstrap ($KIND) -- the only step that uses the GTDB link ==="
    if [ "$KIND" = bac120 ]; then
      say "    bac120_metadata_r226.tsv.gz   ~100 MB   (far bigger than the archaeal one)"
      say "    bac120_r226.tree.gz            ~20 MB   (tree, for the phylogenetic correction)"
      say "    It is slow. wget retries up to 100 times and resumes an interrupted"
      say "    transfer, so re-running the same command never restarts from zero."
    else
      say "    ar53_metadata_r226.tsv.gz  ~5.4 MB   /   ar53_r226.tree.gz  ~90 KB"
    fi
    say ""
    check_tools || exit 1
    m=$(get_meta "$KIND") || { say ""
      say "  Not yet. Two options:"
      say "    1) re-run this command -- the transfer resumes where it had got to"
      say "    2) download the file by hand (see NETWORK.md) and drop it"
      say "       into $D/ : the script recognises the name as it stands"; exit 1; }
    say "  [ok] metadata: $m ($(du -h "$m" | cut -f1))"
    get_tree "$KIND"
    say ""
    if [ "$KIND" = bac120 ]; then
      say "  GTDB is no longer needed. Check the sizing, then run the pilot:"
      say "      ./scripts/download_genomes.sh sizing Bacteria"
      say "      ./scripts/download_genomes.sh pilot 500"
    else
      n=$(reps ar53) || exit 1
      say "  [ok] $n archaeal species representatives extracted"
      say "  GTDB is no longer needed. What follows goes through NCBI:"
      say "      ./scripts/download_genomes.sh archaea"
    fi
    ;;
  archaea)
    say "=== All the archaea ==="
    check_tools || exit 1
    get_tree ar53
    m=$(get_meta ar53) || exit 1
    plan="$D/ar53_plan.tsv"
    if plan_stale "$plan"; then
      [ -s "$plan" ] && say "  plan in the old format (without completeness column) -- rebuilding"
      python3 scripts/sample_design.py --meta "$m" --domain Archaea \
        --stratum family --min-completeness "$MIN_COMP" --out "$plan" || exit 1
    else
      say "  plan already built: $plan"
    fi
    N="${2:-0}"
    if [ "$N" -gt 0 ]; then
      head -n $((N+1)) "$plan" | tail -n +2 | cut -f1,6 > "$D/ar53_reps.tsv"
      say "  requested prefix: $N genomes (the plan stays complete on disk)"
    else
      cut -f1,6 "$plan" | tail -n +2 > "$D/ar53_reps.tsv"
    fi
    n=$(wc -l < "$D/ar53_reps.tsv")
    say "  genomes retained: $n across $(cut -f2 "$D/ar53_reps.tsv" | sort -u | wc -l) families"
    say "  estimated volume (proteins only): ~$((n*436/1000)) MB"
    fetch "$D/ar53_reps.tsv" protein
    ;;
  pilot|bacteria)
    # Balanced sampling plan: one genome per family, families in a reproducible
    # pseudo-random order. See scripts/power.py:
    #   proportional draw   -> N_eff saturates around 450, however many are drawn
    #   1 genome per family -> N_eff about 2400 at N = 4000
    # Every prefix of the list is a valid sample, so the pilot of 500 is not
    # thrown away but forms part of the final sample.
    if [ "$1" = pilot ]; then N="${2:-500}"; else N="${2:-3000}"; fi
    say "=== Balanced bacterial sample (1 per family), N=$N ==="
    check_tools || exit 1
    get_tree bac120
    m=$(get_meta bac120) || exit 1
    plan="$D/bac120_plan.tsv"
    if plan_stale "$plan"; then
      [ -s "$plan" ] && say "  plan in the old format -- rebuilding"
      python3 scripts/sample_design.py --meta "$m" --domain Bacteria \
        --stratum family --min-completeness "$MIN_COMP" --out "$plan" || exit 1
    else
      say "  plan already built: $plan ($(( $(wc -l < "$plan") - 1 )) ordered genomes)"
    fi
    head -n $((N+1)) "$plan" | tail -n +2 | cut -f1,6 > "$D/bac120_sample.tsv"
    say "  sample: $(wc -l < "$D/bac120_sample.tsv") genomes across \
$(cut -f2 "$D/bac120_sample.tsv" | sort -u | wc -l) families"
    if [ "$1" = pilot ]; then
      say ""
      say "  Pilot run. The FAcD prevalence measured here fixes the final N,"
      say "  and results/gtdb/pilot_*.json records it. Extend the sample with"
      say "    ./scripts/download_genomes.sh bacteria <N>"
      say "  The genomes already retrieved are not downloaded again."
    fi
    fetch "$D/bac120_sample.tsv" protein _bac
    say ""
    say "  Bacterial proteomes in data/gtdb/proteomes_bac/ -- separate from the archaea."
    ;;
  paired)
    # Paired test on the provenance of the annotation. Take N genomes that
    # already have a downloaded proteome, retrieve their nucleotide assembly,
    # predict their genes with Prodigal into a separate directory, and scan
    # both. Each genome is then its own control, so any difference comes from
    # the gene calling rather than from the clade, the genome size or the
    # completeness.
    N="${2:-100}"; P="${3:-1259}"
    plan="$D/ar53_plan.tsv"
    [ -s "$plan" ] || { say "  missing $plan"; exit 1; }
    check_tools || exit 1
    mkdir -p "$D/genomes_test" "$D/proteomes_test" results/gtdb
    pred="$D/log/proteomes_predicted.tsv"
    head -n $((P+1)) "$plan" | tail -n +2 | cut -f1 | sort -u > "$D/log/.pref.txt"
    ls "$D/proteomes"/*.faa.gz 2>/dev/null | sed 's|.*/||; s|\.faa\.gz$||' | sort -u > "$D/log/.present.txt"
    if [ -s "$pred" ]; then tail -n +2 "$pred" | cut -f1 | sort -u > "$D/log/.pred.txt"
    else : > "$D/log/.pred.txt"; fi
    # Candidates: in the prefix, proteome present, and not yet predicted.
    comm -12 "$D/log/.pref.txt" "$D/log/.present.txt" | comm -23 - "$D/log/.pred.txt" \
      | head -n "$N" > "$D/paired.txt"
    m=$(wc -l < "$D/paired.txt")
    say "=== Paired test: $m genomes that already have a downloaded proteome ==="
    say "    (~$((m*700/1000)) MB of nucleotides)"
    awk '{print $1"\t-"}' "$D/paired.txt" > "$D/paired_fetch.tsv"
    D_SAVE="$D"
    fetch_dir="$D/genomes_test"
    # Download the nucleotides into a dedicated directory.
    ( D="$D_SAVE"; outdir_override="$fetch_dir"; : )
    say "  downloading the assemblies..."
    rm -f "$D/log/done_genome_test.txt"; touch "$D/log/done_genome_test.txt"
    split -l "$BATCH" "$D/paired.txt" "$D/log/tbatch_"
    i=0; tot=$(ls "$D"/log/tbatch_* | wc -l)
    for b in "$D"/log/tbatch_*; do
      i=$((i+1)); say "   batch $i/$tot"
      if datasets download genome accession --inputfile "$b" --include genome \
           --no-progressbar --filename "$D/tmp_t.zip" >>"$D/log/errors.txt" 2>&1; then
        rm -rf "$D/tmp_t"; unzip -qo "$D/tmp_t.zip" -d "$D/tmp_t" 2>/dev/null
        while IFS= read -r f; do
          [ -s "$f" ] || continue
          g=$(basename "$(dirname "$f")")
          gzip -c "$f" > "$D/genomes_test/${g}.fna.gz" && echo "$g" >> "$D/log/done_genome_test.txt"
        done < <(find "$D/tmp_t" -name "*_genomic.fna" 2>/dev/null)
        rm -rf "$D/tmp_t" "$D/tmp_t.zip"
      else say "   [batch $i failed -- re-run the command]"; fi
      rm -f "$b"
    done
    say "  $(wc -l < "$D/log/done_genome_test.txt") assemblies retrieved"
    say ""
    say "=== Gene prediction, in a separate directory ==="
    python3 scripts/predict_genes.py --genomes "$D/genomes_test" \
      --out "$D/proteomes_test" --manifest "$D/log/predits_test.tsv" || exit 1
    say ""
    say "=== Scan of the predicted proteomes ==="
    python3 scripts/scan_proteomes.py --proteomes "$D/proteomes_test" \
      --hmm results/hmm/Fluc_CrcB.hmm --hmm results/hmm/CLC_F.hmm \
      --out results/gtdb/scan_paired.tsv --jobs 6 --cpu 2 || exit 1
    say ""
    say "  In results/gtdb/scan_paired.tsv each genome is its own control."
    ;;
  missing)
    # Which genomes of the plan have no proteome at all, that is which GenBank
    # assemblies are unannotated. $2 = prefix of the plan;
    # $3 = archaea (default) or bacteria.
    N="${2:-1259}"; DOM="${3:-archaea}"
    if [ "$DOM" = bacteria ] || [ "$DOM" = bacteries ]; then
      plan="$D/bac120_plan.tsv"; SUF="_bac"; BOOTSTRAP="pilot"
    else
      plan="$D/ar53_plan.tsv"; SUF=""; BOOTSTRAP="archaea"
    fi
    [ -s "$plan" ] || { say "  missing $plan -- run '$BOOTSTRAP' first"; exit 1; }
    mkdir -p results/gtdb "$D/proteomes$SUF"
    head -n $((N+1)) "$plan" | tail -n +2 | cut -f1 > "$D/log/.targeted.txt"
    ls "$D/proteomes$SUF"/*.faa.gz 2>/dev/null | sed 's|.*/||; s|\.faa\.gz$||' | sort -u > "$D/log/.present.txt"
    sort -u "$D/log/.targeted.txt" -o "$D/log/.targeted.txt"
    comm -23 "$D/log/.targeted.txt" "$D/log/.present.txt" > "$D/missing.txt"
    m=$(wc -l < "$D/missing.txt")
    say "=== $m genomes of the plan (prefix $N) have no proteome ==="
    say "    they will be retrieved as nucleotides then annotated by Prodigal:"
    say "        ./scripts/download_genomes.sh predict $N $DOM"
    awk 'NR==FNR{m[$1];next} FNR>1 && ($1 in m){c[$3]++} END{for(k in c) print "    "c[k]"\t"k}' \
      "$D/missing.txt" "$plan" | sort -rn | head -20 >&2
    ;;
  predict)
    N="${2:-1259}"; DOM="${3:-archaea}"
    if [ "$DOM" = bacteria ] || [ "$DOM" = bacteries ]; then SUF="_bac"; else SUF=""; fi
    check_tools || exit 1
    "$0" missing "$N" "$DOM" || exit 1
    m=$(wc -l < "$D/missing.txt")
    [ "$m" -eq 0 ] && { say "  nothing to predict."; exit 0; }
    say ""
    say "=== Nucleotides for the $m genomes without a proteome (~$((m*700/1000)) MB) ==="
    awk '{print $1"\t-"}' "$D/missing.txt" > "$D/missing_fetch.tsv"
    fetch "$D/missing_fetch.tsv" genome "$SUF"
    say ""
    say "=== Gene prediction (Prodigal) ==="
    python3 scripts/predict_genes.py --genomes "$D/genomes$SUF" --out "$D/proteomes$SUF" \
      --manifest "$D/log/proteomes_predicted$SUF.tsv" || exit 1
    say ""
    say "  Then re-run the scan: it will only redo the new genomes."
    ;;
  reindex)
    # Rebuild the download journal from the files present. Run this after
    # copying an existing data/ directory in: without it, `fetch` sees an empty
    # journal and downloads everything again.
    say "=== Reindexing the journal from the files present ==="
    for w in protein:faa: genome:fna: protein:faa:_bac genome:fna:_bac; do
      what="${w%%:*}"; rest="${w#*:}"; ext="${rest%%:*}"; suf="${rest#*:}"
      [ "$what" = protein ] && dir="$D/proteomes${suf}" || dir="$D/genomes${suf}"
      [ -d "$dir" ] || continue
      j="$D/log/done_${what}${suf}.txt"; mkdir -p "$D/log"
      n_before=0; [ -s "$j" ] && n_before=$(wc -l < "$j")
      find "$dir" -name "*.${ext}.gz" -size +0 -printf '%f\n' 2>/dev/null \
        | sed "s/\.${ext}\.gz$//" | sort -u > "$j"
      say "  ${what}${suf:+ $suf} : $(wc -l < "$j") genomes indexed (previous journal: $n_before)"
    done
    say "  Now re-run the wanted subcommand: nothing will be re-downloaded."
    ;;
  sizing)
    dom="${2:-Archaea}"
    say "=== Sizing ($dom) -- no download ==="
    mkdir -p results/gtdb
    python3 scripts/plan_sizing.py --meta-dir "$D" --domain "$dom" \
      --out "results/gtdb/sizing_$(echo "$dom" | tr 'A-Z' 'a-z').json" || exit 1
    ;;
  audit)
    # Downloads nothing: measures what a completeness threshold would cost,
    # and which phyla would pay it.
    dom="${2:-Archaea}"
    say "=== Audit of the completeness filter ($dom) -- no download ==="
    mkdir -p results/gtdb
    o="results/gtdb/completeness_$(echo "$dom" | tr 'A-Z' 'a-z').json"
    python3 scripts/completeness_audit.py --meta-dir "$D" --domain "$dom" --out "$o" || exit 1
    ;;
  archaea-ncbi)
    say "=== Fallback route: archaea via NCBI alone (GTDB unreachable) ==="
    check_tools || exit 1
    mkdir -p data/ncbi
    python3 scripts/ncbi_reps.py --taxon archaea --out data/ncbi/ar_reps.tsv || exit 1
    fetch data/ncbi/ar_reps.tsv protein
    ;;
  nucleotides)
    say "=== Nucleotides, for the exporter negatives only ==="
    check_tools || exit 1
    f=results/gtdb/negatives.txt
    [ -s "$f" ] || { say "  missing $f (produced by the pipeline)"; exit 1; }
    awk 'NF{print $1"\t-"}' "$f" > "$D/neg.tsv"
    say "  $(wc -l < "$D/neg.tsv") genomes"
    fetch "$D/neg.tsv" genome
    ;;
  *) usage >&2; exit 1 ;;
esac
