#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the seed of the FAcD equivalog from UniProt, auditably.

  python3 scripts/facd_seed.py --out seeds
Requires network access to UniProt, plus mafft and hmmer.

"Everything that carries EC 3.8.1.3" cannot serve as the positive set, for two
reasons:

  1. Q01399 (DehH2 from Moraxella sp. B) carries EC 3.8.1.3 in the manually
     reviewed part of UniProt, while its own entry states that "dehH2 acts on
     chloro-, bromo- and iodoacetate, but not on fluoroacetate", and it
     belongs to the HAD superfamily (PF00702) rather than to the alpha/beta
     hydrolases (PF00561).
  2. Among the EC 3.8.1.3 entries that are in PF00561, about half lack the
     HDRG nucleophile elbow shared by the three experimentally characterised
     enzymes. Those are epoxide hydrolases and other subfamilies of the same
     fold.

The seed is therefore built the same way as for CLC-F, where the same problem
arose:

  a. three experimentally characterised anchors:
       Q6NAM1  RPA1163      Rhodopseudomonas palustris
       Q1JU72  FAc-DEX FA1  Burkholderia sp.
       Q01398  DehH1        Moraxella sp. B
  b. a preliminary model built on those three anchors alone;
  c. the whole candidate pool scored with that model;
  d. the largest discontinuity in the score distribution located. The
     threshold is not fixed a priori but read from the data, and the three
     largest discontinuities are reported so the choice can be checked;
  e. above the cut, the canonical group, which becomes the seed and the
     validation set; below it, unlabelled sequences, which are neither
     positives nor negatives. Treating them as negatives is what the CLC-F
     model had to avoid.

Q01399 is the negative control and must fall below the cut. If it does not,
the model fails to separate what it has to separate, and the report says so.
"""

import argparse, json, os, re, subprocess, sys, urllib.parse, urllib.request
from collections import defaultdict

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
ANCHORS = ["Q6NAM1", "Q1JU72", "Q01398"]
NEGATIVE_CONTROL = "Q01399"        # DehH2: EC 3.8.1.3, but not on fluoroacetate


def say(*a):
    print(*a, file=sys.stderr)


CACHE = None


def fetch(query, fields=None, fmt="fasta", size=500):
    """UniProt REST query with cursor pagination (the Link header).
    With a cache configured, a pool already downloaded is not fetched again, so
    an interrupted run resumes immediately."""
    if CACHE:
        import hashlib
        f = os.path.join(CACHE, hashlib.md5(query.encode()).hexdigest() + "." + fmt)
        if os.path.exists(f) and os.path.getsize(f) > 0:
            say(f"      (cache) {query[:56]}...")
            return open(f, encoding="utf-8").read()
    p = {"query": query, "format": fmt, "size": str(size)}
    if fields:
        p["fields"] = fields
    url = f"{UNIPROT}?{urllib.parse.urlencode(p)}"
    out = []
    while url:
        req = urllib.request.Request(url, headers={"User-Agent": "facd-seed/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            out.append(r.read().decode("utf-8", "replace"))
            link = r.headers.get("Link", "")
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None
    txt = "".join(out)
    if CACHE:
        import hashlib
        os.makedirs(CACHE, exist_ok=True)
        open(os.path.join(CACHE, hashlib.md5(query.encode()).hexdigest()
                          + "." + fmt), "w", encoding="utf-8").write(txt)
    return txt


def read_fasta(txt):
    """Index on the first word of the header, which is what hmmsearch puts in
    column 0 of --tblout. The full header is kept alongside, because the genus
    can only be read there (OS=...)."""
    d, tok, hdr, buf = {}, None, None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if tok:
                d[tok] = (hdr, "".join(buf))
            hdr = line[1:]
            tok = hdr.split()[0] if hdr.split() else hdr
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if tok:
        d[tok] = (hdr, "".join(buf))
    return d


def acc(tok):
    p = tok.split("|")
    return p[1] if len(p) > 2 else tok


def genus(entry):
    """entry = (full header, sequence)"""
    m = re.search(r"OS=(\S+)", entry[0] or "")
    return m.group(1) if m else "?"


def write_fasta(path, d):
    with open(path, "w") as fh:
        for tok, (hdr, seq) in d.items():
            fh.write(f">{hdr}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    return len(d)


def run(cmd, **kw):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"[FAIL] {cmd}\n{r.stderr[:1500]}")
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="seeds")
    ap.add_argument("--tmp", default="seeds/.facd_tmp")
    ap.add_argument("--cache", default="seeds/.facd_cache",
                    help="cache directory for the UniProt pools; "
                         "emptying this directory forces a fresh download")
    ap.add_argument("--pool-query",
                    default='(ec:3.8.1.3 OR protein_name:"fluoroacetate dehalogenase") '
                            'AND (xref:pfam-PF00561)')
    ap.add_argument("--controls", default="A9BLX5",
                    help="additional accessions whose position is reported "
                         "without using them as anchors (comma-separated list)")
    ap.add_argument("--threshold", type=float, default=None,
                    help="force the partition threshold (default: the largest discontinuity)")
    a = ap.parse_args()
    global CACHE
    CACHE = a.cache
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.tmp, exist_ok=True)

    # 1. Download the UniProt pools.
    say("[1/5] downloading the UniProt pools...")
    pool = read_fasta(fetch(a.pool_query))
    say(f"      candidate pool         : {len(pool)} sequences")
    anchors = read_fasta(fetch(" OR ".join(f"accession:{x}" for x in ANCHORS)))
    say(f"      characterised anchors  : {len(anchors)}")
    neg_control = read_fasta(fetch(f"accession:{NEGATIVE_CONTROL}"))
    say(f"      DehH2 negative control : {len(neg_control)}")
    # Observed controls are never used as anchors; the report only says where
    # they fall. A9BLX5 is an alpha/beta hydrolase from Delftia acidovorans
    # SPH-1, a candidate for the DEF1 used by Dodge et al. 2024 and by Dodge,
    # O'Connor & Wackett 2026. It must be checked against the published
    # sequence before being promoted to an anchor.
    obs = [x.strip() for x in a.controls.split(",") if x.strip()]
    observed = read_fasta(fetch(" OR ".join(f"accession:{x}" for x in obs))) if obs else {}
    say(f"      observed controls      : {len(observed)} ({', '.join(obs)})")
    neg = read_fasta(fetch("(ec:3.8.1.5 OR ec:3.3.2.9 OR ec:3.3.2.10) "
                           "AND (xref:pfam-PF00561)"))
    say(f"      close negatives        : {len(neg)} "
        f"(haloalkane dehalogenases, epoxide hydrolases)")
    broad = read_fasta(fetch("(xref:pfam-PF00561) AND (reviewed:true) "
                             "NOT (ec:3.8.1.3) NOT (ec:3.8.1.5)"))
    say(f"      broad negatives        : {len(broad)} (PF00561 reviewed, outside 3.8.1.*)")
    if len(anchors) < 3:
        sys.exit("[FAIL] the three anchors were not retrieved")

    # 2. Preliminary model, from the three anchors alone.
    say("[2/5] preliminary model from the three anchors...")
    fa = f"{a.tmp}/anchors.faa"; write_fasta(fa, anchors)
    run(f"mafft --auto --quiet {fa} > {a.tmp}/anchors.aln")
    run(f"hmmbuild --amino -n FAcD_prelim {a.tmp}/prelim.hmm {a.tmp}/anchors.aln")

    # 3. Score the candidate pool with that model.
    say("[3/5] scoring of the candidate pool...")
    all_seqs = dict(pool); all_seqs.update(anchors); all_seqs.update(neg_control)
    all_seqs.update(observed)
    fp = f"{a.tmp}/pool.faa"; write_fasta(fp, all_seqs)
    run(f"hmmsearch --max -E 1000 --tblout {a.tmp}/pool.tbl -o /dev/null "
        f"{a.tmp}/prelim.hmm {fp}")
    score = {}
    for line in open(f"{a.tmp}/pool.tbl"):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) > 5:
            score[f[0]] = max(score.get(f[0], -1e9), float(f[5]))
    # A sequence the model does not reach even at E=1000 is not an anomaly: it
    # is the clearest available signal that the sequence belongs to another
    # family. It is scored -inf so that it sorts to the bottom.
    silent = [e for e in all_seqs if e not in score]
    for e in silent:
        score[e] = float("-inf")
    say(f"      {len(all_seqs) - len(silent)} sequences scored; "
        f"{len(silent)} without any hit (score -inf)")
    for e in silent[:15]:
        say(f"        no hit: {acc(e)}")
    if len(silent) > 15:
        say(f"        ... and {len(silent) - 15} others (full list in the report)")
    if len(silent) > 0.5 * len(all_seqs):
        say("      [ALERT] more than half of the pool obtains no hit at all; "
            "check the query and the tools before reading the result.")

    # 4. Read the threshold off the score distribution.
    say("[4/5] search for the discontinuity...")
    ranked = sorted(((s, e) for e, s in score.items()), reverse=True)
    n = len(ranked)
    # An equivalog has to contain its own characterised members, so the cut is
    # searched for strictly below the lowest anchor score. Without that
    # constraint the largest discontinuity could fall among the anchors.
    s_anchors = [score.get(k, float("-inf")) for k in all_seqs if acc(k) in ANCHORS]
    floor = min(s_anchors) if s_anchors else float("inf")
    say(f"      anchor scores : "
        + ", ".join(f"{round(x,1)}" for x in sorted(s_anchors, reverse=True))
        + f"  -> the cut will be looked for below {floor:.1f} bits")
    bounds = [(ranked[i][0] - ranked[i + 1][0], ranked[i][0], ranked[i + 1][0], i + 1)
              for i in range(n - 1)
              if ranked[i][0] < floor and ranked[i + 1][0] > float("-inf")]
    if not bounds:
        sys.exit("[FAIL] no usable discontinuity below the anchors")
    bounds.sort(reverse=True)
    say("      three largest discontinuities:")
    for d, above, below, i in bounds[:3]:
        say(f"        gap {d:7.1f} bits between {above:7.1f} and {below:7.1f} "
            f"-> {i} sequences above")
    threshold = a.threshold if a.threshold is not None else (bounds[0][1] + bounds[0][2]) / 2
    say(f"      threshold retained : {threshold:.1f} bits"
        + ("  (imposed on the command line)" if a.threshold is not None else
           "  (midpoint of the largest discontinuity)"))

    above = {e: all_seqs[e] for s, e in ranked if s >= threshold and e in all_seqs}
    below = {e: all_seqs[e] for s, e in ranked if s < threshold and e in all_seqs}

    # 5. Check the controls and write the output.
    say("[5/5] controls...")
    ctrl = {}
    for name, lst in (("anchors", ANCHORS), ("negative_control", [NEGATIVE_CONTROL]),
                      ("observed_controls", obs)):
        for x in lst:
            e = next((k for k in all_seqs if acc(k) == x), None)
            s = score.get(e) if e else None
            side = "ABOVE" if (s is not None and s >= threshold) else "below"
            shown = "absent" if s is None else ("no hit" if s == float("-inf")
                                                else f"{s:.1f}")
            ctrl[x] = dict(score=None if s in (None, float("-inf")) else s,
                           hit=s is not None and s != float("-inf"), side=side)
            say(f"      {x:<8} score {shown:>10}  -> {side}")
    ok = all(ctrl[x]["side"] == "ABOVE" for x in ANCHORS) and ctrl[NEGATIVE_CONTROL]["side"] == "below"
    say("      controls: " + ("the three anchors fall above the cut and DehH2 "
                              "below it, so the split separates what it must"
                              if ok else
                              "the controls do not fall on the expected sides, "
                              "so the split is not usable"))

    # One sequence per genus goes into the seed; the rest is held for validation.
    by_genus = defaultdict(list)
    for e, v in above.items():
        by_genus[genus(v)].append(e)
    seed, valid = {}, {}
    for g, es in sorted(by_genus.items()):
        es.sort()
        seed[es[0]] = above[es[0]]
        for e in es[1:]:
            valid[e] = above[e]

    negatives = dict(neg); negatives.update(broad); negatives.update(neg_control)
    for e in list(observed):            # an observed control is not a negative
        negatives.pop(e, None)
    for e in below:                    # unlabelled sequences are not negatives
        negatives.pop(e, None)

    r = {
        "seed": write_fasta(f"{a.out}/facd_pos.faa", seed),
        "validation": write_fasta(f"{a.out}/facd_truth.faa", valid),
        "unlabelled": write_fasta(f"{a.out}/facd_unlabelled.faa", below),
        "negatives": write_fasta(f"{a.out}/facd_neg.faa", negatives),
    }
    say("")
    say(f"      seed            : {r['seed']} sequences, {len(by_genus)} genera")
    say(f"      validation      : {r['validation']} sequences")
    say(f"      unlabelled      : {r['unlabelled']} (neither positives nor negatives)")
    say(f"      negatives       : {r['negatives']}")

    json.dump(dict(threshold=threshold, discontinuities=[dict(gap=d, above=h, below=b, n_above=i)
                                                         for d, h, b, i in bounds[:5]],
                   controls=ctrl, split_valid=ok, counts=r,
                   genera=sorted(by_genus), query=a.pool_query,
                   no_hit=[acc(e) for e in silent],
                   scores={acc(e): round(s, 1) for s, e in ranked}),
              open(f"{a.out}/facd_seed_report.json", "w"), indent=2, ensure_ascii=False)
    say(f"\n[ok] {a.out}/facd_seed_report.json")
    if not ok:
        say("[WARNING] the split is invalid, but the report was written all the same.")


if __name__ == "__main__":
    main()
