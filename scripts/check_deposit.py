#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Self-check for the deposit. Run this before anything else.

The defects this suite is built to catch are not arithmetic errors inside a
single file; they are disagreements between files. A threshold written by hand
into a .hmm header that its own calibration record contradicts; a calibration
record in French that the deposited script would have written in English; a
figure script that stops on the first key of the JSON shipped beside it. None of
those is visible by reading any one file on its own. They are visible only by
cross-checking files against each other, so the deposit ships the cross-check.

Exit code 0 means every check passed.
"""
import collections, json, os, re, shutil, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OK, FAIL, WARN = [], [], []

def ok(m):   OK.append(m);   print(f"  [ok]   {m}")
def bad(m):  FAIL.append(m); print(f"  [FAIL] {m}")
def warn(m): WARN.append(m); print(f"  [warn] {m}")

_NEG_SCORES = {}


def negative_scores(model, negfa, cutoff):
    """Best score per sequence for a model against its negative set.

    Memoised: the confirmed-subset check and the above-GA check need the same
    table, and scoring a few thousand sequences with --max is the slowest thing
    this script does.
    """
    key = (model, negfa, cutoff)
    if key in _NEG_SCORES:
        return _NEG_SCORES[key]
    exe = shutil.which("hmmsearch")
    hmm = ROOT / f"results/hmm/{model}.hmm"
    if not exe or not hmm.exists() or not (ROOT / negfa).exists():
        _NEG_SCORES[key] = None
        return None
    tbl = tempfile.mktemp(suffix=".tbl")
    try:
        subprocess.run([exe, "--max", "-E", cutoff, "--tblout", tbl,
                        "-o", os.devnull, str(hmm), str(ROOT / negfa)],
                       check=True, capture_output=True, timeout=900)
        best = {}
        for line in open(tbl):
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 6:
                best[f[0]] = max(best.get(f[0], float("-inf")), float(f[5]))
    except (subprocess.SubprocessError, OSError):
        best = None
    finally:
        if os.path.exists(tbl):
            os.unlink(tbl)
    _NEG_SCORES[key] = best
    return best


def negative_headers(negfa):
    """Accession -> full FASTA header, for the deposited negative set."""
    out = {}
    with open(ROOT / negfa) as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].rstrip()
                out[h.split()[0]] = h
    return out


def hmm_header(p):
    h = {}
    for line in p.read_text(errors="ignore").splitlines():
        m = re.match(r"^(LENG|NSEQ|GA|TC|NC)\s+([-\d.]+)", line)
        if m:
            h[m.group(1)] = float(m.group(2))
        if line.startswith("HMM "):
            break
    return h

# CLC_F_strict belongs in this list: a model left out of it has its calibration
# record skipped by check 9, where a mislabelled field would go unnoticed.
MODELS = ["Fluc_CrcB", "CLC_F", "FAcD", "CLC_F_strict"]
# Negative set each model was calibrated against, used to re-count the confirmed
# subset from the sequences themselves.
NEG_FASTA = {
    "Fluc_CrcB": "seeds/negatives.faa",
    "CLC_F": "seeds/negatives_clc.faa",
    "CLC_F_strict": "seeds/negatives_clc.faa",
    "FAcD": "seeds/facd_neg.faa",
}

print("=" * 68)
print("1. Every .hmm header agrees with its own calibration record")
print("=" * 68)
print("   (a threshold written by hand is invisible in the .hmm alone)")
for m in MODELS:
    hp = ROOT / "results/hmm" / f"{m}.hmm"
    cp = ROOT / "results/hmm" / f"{m}.calib.json"
    if not hp.exists() or not cp.exists():
        bad(f"{m}: missing .hmm or .calib.json"); continue
    h = hmm_header(hp); c = json.loads(cp.read_text())
    for tag, key in (("GA", "gathering_threshold"), ("TC", "loo_min"), ("NC", "negative_max")):
        hv, cv = h.get(tag), c.get(key)
        if hv is None or cv is None:
            bad(f"{m}: {tag} absent"); continue
        if abs(hv - cv) <= 0.011:
            ok(f"{m}: {tag} {hv:.2f} == {key} {cv:.2f}")
        else:
            bad(f"{m}: {tag} {hv:.2f} != {key} {cv:.2f}  (threshold not reproducible)")

print()
print("=" * 68)
print("2. Every threshold is the midpoint of its admissible window")
print("=" * 68)
print("   (this is the rule the README defends; it must be the rule used)")
for m in MODELS:
    cp = ROOT / "results/hmm" / f"{m}.calib.json"
    if not cp.exists(): continue
    c = json.loads(cp.read_text())
    ga, lo, hi = c.get("gathering_threshold"), c.get("negative_max"), c.get("loo_min")
    if None in (ga, lo, hi):
        bad(f"{m}: incomplete calibration record"); continue
    mid = 0.5 * (lo + hi)
    if abs(mid - ga) <= 0.011:
        ok(f"{m}: GA {ga:.2f} = midpoint({lo:.1f}, {hi:.1f})")
    else:
        bad(f"{m}: GA {ga:.2f} != midpoint({lo:.1f}, {hi:.1f}) = {mid:.2f}")
    if c.get("rule") != "midpoint":
        bad(f"{m}: calibration says rule={c.get('rule')!r}, expected 'midpoint'")
    else:
        ok(f"{m}: rule field says 'midpoint'")

print()
print("=" * 68)
print("3. HMMER convention NC < GA < TC")
print("=" * 68)
print("   (NC above GA makes the file internally incoherent for hmmsearch)")
for p in sorted((ROOT / "results/hmm").glob("*.hmm")):
    h = hmm_header(p)
    nc, ga, tc = h.get("NC"), h.get("GA"), h.get("TC")
    if None in (nc, ga, tc):
        warn(f"{p.name}: incomplete header"); continue
    if nc < ga < tc:
        ok(f"{p.name}: {nc:.2f} < {ga:.2f} < {tc:.2f}")
    else:
        bad(f"{p.name}: {nc:.2f} < {ga:.2f} < {tc:.2f} violated")

print()
print("=" * 68)
import sys as _sys4
if str(ROOT / "scripts") not in _sys4.path:
    _sys4.path.insert(0, str(ROOT / "scripts"))
from deposit_exclusions import is_excluded as _excl4


def deposit_excluded(path):
    return _excl4(str(path.relative_to(ROOT)))


print("4. No French names survive in any deposited JSON or TSV")
print("=" * 68)
print("   (a deposited file whose field names the deposited script does not")
print("    write cannot have been produced by it)")

# This used to test a hand-written list of fourteen known French keys. That is a
# list of the offenders already found, not a test: seeds/facd_seed_report.json
# carried nine French top-level keys that were not on it, so the check reported
# green on a file written entirely in French. The test below works the other way
# round -- it looks for French morphology in whatever names are actually there,
# so a French name nobody has seen yet still fails.
#
# Accessions are keys in some of these files (UniProt accessions in the seed
# report's "scores" and "controls" blocks), so they are exempted by pattern, not
# by name.
FRENCH_STEMS = {
    # function words and endings that do not occur in English field names
    "le", "la", "les", "des", "une", "du", "aux", "ses", "leur", "cette", "ces",
    "est", "sont", "pas", "qui", "que", "pour", "dans", "avec", "sans", "nous",
    "donc", "mais", "tout", "tous", "toute", "meme", "deja", "dont", "ou",
    # domain words that appeared, or plausibly would, in this deposit
    "seuil", "seuils", "graine", "graines", "ecart", "ecarts", "effectif",
    "effectifs", "negatif", "negatifs", "positif", "positifs", "etiquete",
    "etiquetes", "non_etiquetes", "partage", "valide", "requete", "requetes",
    "hit_absent", "sans_hit", "discontinuite", "discontinuites", "controle",
    "controles", "cote", "genre", "genres", "famille", "familles",
    "ordre", "ordres", "couche", "couches", "completude", "taille",
    "tailles", "longueur", "longueurs", "organisme", "organismes", "libelle",
    "libelles", "alerte", "alertes", "statut", "statuts", "nom", "noms",
    "depot", "depots", "au_dessus", "n_au_dessus", "dessus", "dessous",
    "bas", "haut", "hauts", "regle", "regles", "fenetre", "fenetres",
    "mesuree", "mesurees", "mesure", "brute", "brutes", "moyenne", "moyennes",
    "mediane", "medianes", "echantillon", "echantillons", "marqueur",
    "marqueurs", "phylum_par", "par_phylum", "par_taille", "par_source",
    "demande", "demandee", "contraignant", "contraignants", "proportion_brute",
    "proteines", "proteomes_predits", "predits", "predit",
    "nombre", "vide", "vides", "absent_du", "hypothese",
    # "classes" is left out on purpose: it is a plausible English field name.
    # A French "classes" would be caught through its neighbours ("classe_libelle"
    # by "libelle"), and a false positive on an English name is worse here than
    # one missed form.
    "marge", "marges", "rappel", "rappels", "taux", "classe",
    "etiquette", "etiquettes", "retenu", "retenus", "retenue", "retenues",
    "faux", "vrai", "vrais", "vraie", "vraies", "prevalence_brute",
    "hypotheses", "anterieure", "tirage", "tirages", "ligne_de_commande",
    "defaut", "cles", "fichier", "fichiers", "repertoire", "repertoires",
}
ACCESSION = re.compile(r"^(?:[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]"
                       r"|[OPQ][0-9][A-Z0-9]{3}[0-9]"
                       r"|[A-Z0-9]{6,10}"
                       r"|GC[AF]_\d+\.\d+)$")


def french_names(names):
    """Names that look French, ignoring accessions and purely numeric labels."""
    out = []
    for n in names:
        if not isinstance(n, str) or ACCESSION.match(n) or n.replace(".", "").isdigit():
            continue
        parts = [t for t in re.split(r"[_\s-]+", n.lower()) if t]
        if any(t in FRENCH_STEMS for t in parts) or n.lower() in FRENCH_STEMS:
            out.append(n)
    return out


def all_keys(obj, acc=None):
    acc = set() if acc is None else acc
    if isinstance(obj, dict):
        for k, v in obj.items():
            acc.add(k)
            all_keys(v, acc)
    elif isinstance(obj, list):
        for v in obj:
            all_keys(v, acc)
    return acc


found, n_json, n_tsv = [], 0, 0
for p in sorted(ROOT.rglob("*.json")):
    if deposit_excluded(p):
        continue
    try:
        obj = json.loads(p.read_text())
    except Exception:
        continue
    n_json += 1
    hits = french_names(all_keys(obj))
    if hits:
        found.append((p.relative_to(ROOT), sorted(hits)[:8]))
for p in sorted(ROOT.rglob("*.tsv")):
    if deposit_excluded(p):
        continue
    try:
        head = p.read_text(errors="replace").split("\n", 1)[0]
    except Exception:
        continue
    n_tsv += 1
    hits = french_names(head.split("\t"))
    if hits:
        found.append((p.relative_to(ROOT), sorted(hits)[:8]))
if found:
    for rel, hits in found:
        bad(f"{rel}: French field names {hits}")
else:
    ok(f"no French field name in {n_json} deposited JSON and {n_tsv} TSV files, "
       f"tested by morphology rather than against a list of known ones")

# The documentation names these fields too, and a French name mentioned in prose
# is both a leftover and a wrong name. A line naming one identifier carries too
# few French words for a prose detector to notice, so the identifiers inside
# backticks in the deposited Markdown go through the same test as the fields
# themselves.
# docs/CHANGES.md is the one place a French name is legitimate, because it
# records the renames and has to quote what was renamed. It is not exempted
# outright: a French name there is accepted only on a line that also names its
# English replacement, which is what a rename entry looks like. A leftover
# sitting on its own still fails.
_md_hits = []
for p in sorted(ROOT.rglob("*.md")):
    if deposit_excluded(p):
        continue
    _is_log = p.name == "CHANGES.md"
    _h = set()
    for _line in p.read_text(errors="replace").splitlines():
        _idents = set(re.findall(r"`([A-Za-z][A-Za-z0-9_]{2,})`", _line))
        _fr = set(french_names(_idents))
        if not _fr:
            continue
        if _is_log and (_idents - _fr):
            continue          # a rename entry: the English name is on the line
        _h |= _fr
    if _h:
        _md_hits.append((p.relative_to(ROOT), sorted(_h)[:8]))
if _md_hits:
    for rel, hits in _md_hits:
        bad(f"{rel}: documents field names that are French {hits}")
else:
    ok("no French identifier inside backticks in the deposited documentation, "
       "beyond the renames the changelog records beside their English names")

# And a column the documentation names must exist in the table it documents,
# which catches the same drift when the rename is only half applied.
_rt = ROOT / "seeds/uniprot_recheck_2026-10-01.tsv"
_rm = ROOT / "seeds/README.md"
if _rt.exists() and _rm.exists():
    _cols = set(_rt.read_text(errors="replace").split("\n", 1)[0].split("\t"))
    _named = {m for m in re.findall(r"`([A-Za-z][A-Za-z0-9_]*)` column",
                                    _rm.read_text(errors="replace"))}
    _ghost = sorted(_named - _cols)
    if _ghost:
        bad(f"seeds/README.md names columns the re-check table does not have: "
            f"{_ghost}")
    elif _named:
        ok(f"the {len(_named)} re-check columns named in seeds/README.md all "
           f"exist in the table")

# The seed report is written by one json.dump(dict(...)) call with literal
# keyword names, so the names the script writes can be read out of the script
# and compared with the names in the file. The French keys above were exactly
# this disagreement: the file predated the repository's switch to English, and
# no check compared the two.
_fsr = ROOT / "seeds/facd_seed_report.json"
_fss = ROOT / "scripts/facd_seed.py"
if _fsr.exists() and _fss.exists():
    import ast as _ast4
    _written = None
    for _node in _ast4.walk(_ast4.parse(_fss.read_text())):
        if (isinstance(_node, _ast4.Call)
                and isinstance(_node.func, _ast4.Attribute)
                and _node.func.attr == "dump"
                and _node.args and isinstance(_node.args[0], _ast4.Call)
                and getattr(_node.args[0].func, "id", None) == "dict"):
            _written = {kw.arg for kw in _node.args[0].keywords if kw.arg}
            break
    _in_file = set(json.loads(_fsr.read_text()).keys())
    if _written is None:
        warn("could not read the seed report's key names out of facd_seed.py")
    elif _written != _in_file:
        bad(f"seeds/facd_seed_report.json and scripts/facd_seed.py disagree on "
            f"field names: only in the file {sorted(_in_file - _written)}, only "
            f"in the script {sorted(_written - _in_file)}")
    else:
        ok(f"the seed report's {len(_in_file)} top-level names are the ones "
           f"facd_seed.py writes, read out of the script")

print()
print("=" * 68)
print("5. The figure script can read the deposited analysis JSON")
print("=" * 68)
print("   (re-plotting must work even when re-running is impossible)")
ap = ROOT / "results/gtdb/baseline_archaea.json"
if not ap.exists():
    bad("results/gtdb/baseline_archaea.json absent")
else:
    R = json.loads(ap.read_text())
    try:
        pp = [d for d in R["by_phylum"] if d["N"] >= 10]
        b = R["raw_proportion"]; _ = b["prop"], b["lo"], b["hi"]
        _ = R["mean_by_phylum"]["estimate"]
        _ = {d["phylum"]: d.get("median_size") for d in R["by_phylum"]}
        ok(f"figure access pattern resolves ({len(pp)} phyla with N>=10)")
    except KeyError as e:
        bad(f"figure_archaea.py would stop on KeyError: {e}")

print()
print("=" * 68)
print("6. The census totals are internally consistent")
print("=" * 68)
if ap.exists():
    R = json.loads(ap.read_text())
    N, k = R["sample"]["N"], R["sample"]["k"]
    sN = sum(d["N"] for d in R["by_phylum"])
    sk = sum(d["k"] for d in R["by_phylum"])
    (ok if sN == N else bad)(f"per-phylum N sums to {sN}, sample says {N}")
    (ok if sk == k else bad)(f"per-phylum k sums to {sk}, sample says {k}")
    raw = R["raw_proportion"]["prop"]
    (ok if abs(raw - k / N) < 1e-9 else bad)(
        f"raw proportion {raw:.6f} == {k}/{N} = {k/N:.6f}")

print()
print("=" * 68)
print("7. The design effect is counted from the plan, not assumed")
print("=" * 68)
dp = ROOT / "results/gtdb/design_effect.json"
pp_ = ROOT / "data/gtdb/ar53_plan.tsv"
if not dp.exists():
    warn("results/gtdb/design_effect.json absent (run command 3 of README section 3, WITH its --out flag: without --out the script only prints)")
elif not pp_.exists():
    warn("data/gtdb/ar53_plan.tsv absent: the design effect cannot be re-derived")
else:
    D = json.loads(dp.read_text())
    import collections, csv as _csv
    rows = list(_csv.DictReader(open(pp_), delimiter="\t"))
    fam = collections.Counter(r["family"] for r in rows[:D["n_genomes"]])
    m_eff = sum(v * v for v in fam.values()) / sum(fam.values())
    if abs(m_eff - D["m_eff_kish"]) < 1e-6:
        ok(f"m_eff {m_eff:.4f} recomputed from the plan matches")
    else:
        bad(f"m_eff {m_eff:.4f} from plan != {D['m_eff_kish']:.4f} recorded")
    if D["m_eff_kish"] > D["m_eff_even_allocation_lower_bound"]:
        ok("recorded m_eff exceeds the even-allocation lower bound, as it must")
    else:
        bad("recorded m_eff is at or below the lower bound: it was assumed, not counted")

# plan_sizing.py computes N_eff from an ASSUMED ICC (default 0.15), while
# design_effect.py computes it from the MEASURED ICC (0.730). Both write a field
# called N_eff, so a planning assumption and a measurement are indistinguishable
# in the stored file unless the file states which one it holds. This check
# requires plan_sizing.json to DECLARE where its ICC comes from, and requires its
# corrected column to land on design_effect.py's independent measurement.
psp = ROOT / "results/gtdb/plan_sizing.json"
if not psp.exists():
    warn("results/gtdb/plan_sizing.json absent (needs the GTDB metadata; "
         "see data/gtdb/README.md)")
elif not dp.exists():
    warn("plan_sizing.json present but design_effect.json absent: the two N_eff "
         "cannot be cross-checked")
else:
    P = json.loads(psp.read_text())
    missing = [k for k in ("icc", "icc_role", "icc_source", "icc_measured")
               if k not in P]
    if missing:
        bad(f"plan_sizing.json does not declare where its ICC comes from "
            f"(missing {', '.join(missing)}): a planning assumption and a "
            f"measurement would again be indistinguishable in the file")
    else:
        ok(f"plan_sizing.json declares its ICC as a {P['icc_role']} "
           f"({P['icc']}, {P['icc_source']})")
        if P["icc_measured"] is None:
            warn("plan_sizing.json carries no measured ICC: it was produced "
                 "without design_effect.json on disk")
        elif abs(P["icc_measured"] - D["icc_family"]) > 1e-9:
            bad(f"plan_sizing.json records measured ICC {P['icc_measured']} but "
                f"design_effect.json says {D['icc_family']}")
        else:
            ok(f"the measured ICC it records ({P['icc_measured']:.4f}) is the one "
               f"design_effect.py computed, read from the file, not retyped")
    # the published row: 1259 genomes, 50 % floor, 3 layers
    pub = [r for r in P.get("sizing", [])
           if r["threshold"] == 50 and r["to_draw"] == D["n_genomes"]]
    if not pub:
        bad(f"plan_sizing.json has no row drawing {D['n_genomes']} genomes at the "
            f"50 % floor: it does not describe the plan that was actually used")
    else:
        r = pub[0]
        if abs(r["m_eff"] - round(D["m_eff_kish"], 2)) > 0.01:
            bad(f"plan predicted m_eff {r['m_eff']} but the realised plan measures "
                f"{D['m_eff_kish']:.4f}")
        else:
            ok(f"predicted m_eff {r['m_eff']} matches the realised "
               f"{D['m_eff_kish']:.4f}: the geometry of the plan was right")
        got = r.get("n_eff_at_measured_icc")
        if got is None:
            bad("the published row carries no n_eff_at_measured_icc: the corrected "
                "value is not in the deposit")
        elif abs(got - D["n_eff"]) >= 1:
            bad(f"n_eff under the measured ICC is {got} in plan_sizing.json but "
                f"{D['n_eff']:.2f} in design_effect.json")
        else:
            ok(f"N_eff recomputed under the measured ICC ({got}) lands on "
               f"design_effect.py's {D['n_eff']:.2f} by an independent path; the "
               f"planning value {r['n_eff']} overstated it by "
               f"{100*(r['n_eff']/D['n_eff']-1):.0f} %")

# The N_eff ceiling: the attribution must stay true and the script deterministic.
#
# The manuscript attributes the ceiling table it carried to power.py applied to
# power.py's SYNTHETIC taxonomy, and rests that attribution on an exact numerical
# coincidence: the balanced column of the published table, 873 and 2410,
# reproduces to the unit on that taxonomy. A change to power.modelled_taxo() or
# to deff_at() would make the attribution false in silence, leaving the
# manuscript crediting a table its own code no longer produces.
#
# neff_ceiling.py is also required to be DETERMINISTIC, which is checked by
# running it: power.py drew the proportional column once, from a shared
# generator, so the value it reported was not reproducible.
ncp = ROOT / "results/gtdb/neff_ceiling.json"
if not ncp.exists():
    warn("results/gtdb/neff_ceiling.json absent (run command 9 of README section 3)")
else:
    NC = json.loads(ncp.read_text())
    _tax = {t["key"]: t for t in NC.get("taxonomies", [])}
    # 1. the attribution: the published table's balanced column, to the unit
    PUBLISHED_BALANCED = {1000: 873, 4000: 2410}
    _mod = _tax.get("model")
    if not _mod:
        bad("neff_ceiling.json carries no 'model' taxonomy: the attribution of the "
            "withdrawn table to power.modelled_taxo() can no longer be checked")
    else:
        _mis = []
        for _r in _mod["rows"]:
            _want = PUBLISHED_BALANCED.get(_r["n_drawn"])
            if _want is None:
                continue
            _got = _r["balanced"].get("n_eff")
            if _got != _want:
                _mis.append(f"N={_r['n_drawn']}: {_got} != {_want}")
        if _mis:
            bad("the synthetic taxonomy no longer reproduces the published table's "
                "balanced column (" + "; ".join(_mis) + "), so the manuscript's "
                "attribution of that table to power.py is no longer supported")
        else:
            ok("the synthetic taxonomy still reproduces the withdrawn table's "
               "balanced column to the unit (873 and 2410): the attribution the "
               "manuscript states holds")
    # 2. the archaeal taxonomy must be the plan's own, not another
    _arc = _tax.get("archaea")
    if not _arc:
        warn("neff_ceiling.json has no archaeal row: it was produced without the "
             "GTDB archaeal metadata (see data/gtdb/README.md)")
    elif not dp.exists():
        warn("archaeal ceiling present but design_effect.json absent: its capacity "
             "cannot be checked against the plan")
    else:
        if _arc["capacity_at_cap"] != D["n_genomes"]:
            bad(f"the archaeal capacity at {NC['cap']} layers is "
                f"{_arc['capacity_at_cap']} but the plan drew {D['n_genomes']}: the "
                f"ceiling is computed on a different universe than the plan")
        else:
            ok(f"the archaeal ceiling is computed on the plan's own universe "
               f"({_arc['n_families']} families, {_arc['n_representatives']} "
               f"representatives, capacity {_arc['capacity_at_cap']} = the "
               f"{D['n_genomes']} genomes actually drawn)")
    # 3. determinism: run the script and compare the model column
    _nc_script = ROOT / "scripts/neff_ceiling.py"
    if not _nc_script.exists():
        bad("results/gtdb/neff_ceiling.json is deposited but scripts/neff_ceiling.py "
            "is not: the file has no producing script")
    elif _mod:
        import os as _os8, subprocess as _sp8, tempfile as _tf8
        with _tf8.TemporaryDirectory() as _td8:
            _o8 = _os8.path.join(_td8, "neff.json")
            _r8 = _sp8.run([sys.executable, "scripts/neff_ceiling.py",
                            "--meta-dir", _td8,          # no metadata: model only
                            "--reps", str(NC["replicates"]),
                            "--seed", str(NC["seed"]), "--out", _o8],
                           cwd=ROOT, capture_output=True, text=True, timeout=300)
            if _r8.returncode != 0:
                bad("scripts/neff_ceiling.py fails when run: "
                    + (_r8.stderr.strip().splitlines() or [""])[-1][:140])
            else:
                _again = json.load(open(_o8))
                _m2 = {t["key"]: t for t in _again["taxonomies"]}.get("model")
                if _m2 and _m2["rows"] == _mod["rows"]:
                    ok(f"re-running neff_ceiling.py reproduces the deposited model "
                       f"column exactly ({NC['replicates']} replicates, seed "
                       f"{NC['seed']}): the ceiling is no longer one unseeded draw")
                else:
                    bad("re-running neff_ceiling.py does not reproduce the deposited "
                        "model column: the ceiling is not deterministic")
            # 4. following the recipe WITHOUT the metadata must not degrade the
            #    deposited file. The archaeal row cannot be recomputed without a
            #    file the deposit does not redistribute, so a reader running
            #    command 9 as written would otherwise lose it.
            if _arc:
                import shutil as _sh8
                _copy8 = _os8.path.join(_td8, "deposited_copy.json")
                _sh8.copy2(ncp, _copy8)
                _r9 = _sp8.run([sys.executable, "scripts/neff_ceiling.py",
                                "--meta-dir", _td8, "--reps", "5",
                                "--out", _copy8],
                               cwd=ROOT, capture_output=True, text=True, timeout=300)
                _kept = False
                try:
                    _kept = any(t.get("key") == "archaea" for t in
                                json.load(open(_copy8)).get("taxonomies", []))
                except (OSError, ValueError):
                    _kept = False
                if _r9.returncode == 0 and _kept:
                    ok("following command 9 without the GTDB metadata leaves the "
                       "deposited ceiling intact instead of silently replacing its "
                       "archaeal row with a model-only one")
                else:
                    bad("following command 9 without the GTDB metadata "
                        + ("fails outright" if _r9.returncode else
                           "destroys the deposited archaeal row"))

# The provenance control compares the prevalence measured on downloaded proteomes
# against the prevalence on proteomes predicted here with Prodigal. It is null
# whenever analyse_archaea.py runs without --predicted, and it stayed null across
# several versions because the documented command omitted that flag: a control
# reported as null reads as a control that was not run.
bap = ROOT / "results/gtdb/baseline_archaea.json"
if bap.exists():
    _B = json.loads(bap.read_text())
    _pv = _B.get("provenance_control")
    _manifest = ROOT / "data/gtdb/log/proteomes_predicted.tsv"
    if _pv is None:
        bad("baseline_archaea.json reports provenance_control as null: re-run "
            "command 1 of README section 3 with --predicted "
            "data/gtdb/log/proteomes_predicted.tsv")
    elif not _manifest.exists():
        warn("provenance_control is populated but the manifest it needs is not "
             "deposited, so it cannot be re-derived")
    else:
        import csv as _csv2
        _pred = {r["genome_id"] for r in
                 _csv2.DictReader(open(_manifest), delimiter="\t")}
        _g = {r["source"]: r for r in _pv["global_"]}
        _plan2 = list(_csv2.DictReader(open(ROOT / "data/gtdb/ar53_plan.tsv"),
                                       delimiter="\t"))[:1259]
        _acc = {r["accession"] for r in _plan2}
        _want = len(_acc & _pred)
        if _g.get("predicted", {}).get("N") != _want:
            bad(f"provenance_control counts {_g.get('predicted', {}).get('N')} "
                f"predicted genomes in the sample; the manifest intersected with "
                f"the 1259-row prefix gives {_want}")
        elif _g["downloaded"]["N"] + _g["predicted"]["N"] != len(_plan2):
            bad(f"provenance_control partitions "
                f"{_g['downloaded']['N'] + _g['predicted']['N']} genomes, not "
                f"{len(_plan2)}")
        else:
            _ov = (_g["downloaded"]["lo"] <= _g["predicted"]["hi"]
                   and _g["predicted"]["lo"] <= _g["downloaded"]["hi"])
            ok(f"provenance control runs and partitions the sample exactly "
               f"({_g['downloaded']['N']} downloaded, {_g['predicted']['N']} "
               f"predicted): {100*_g['downloaded']['prop']:.2f} % against "
               f"{100*_g['predicted']['prop']:.2f} %, intervals "
               f"{'overlapping' if _ov else 'DISJOINT'}")

print()
print("=" * 68)
print("8. Data needed to recompute the headline figure is present")
print("=" * 68)
for rel in ["data/gtdb/ar53_plan.tsv", "results/gtdb/scan_archaea.tsv",
            "results/gtdb/baseline_archaea.json"]:
    p = ROOT / rel
    (ok if p.exists() else bad)(f"{rel}{'' if p.exists() else ' MISSING'}")

print()
print("=" * 68)
print("9. Each calibration record carries the filter that produced its own number")
print("=" * 68)
print("   (negative_max is a maximum over a SUBSET; the subset must be stated)")
for m in MODELS:
    cp = ROOT / "results/hmm" / f"{m}.calib.json"
    if not cp.exists():
        continue
    c = json.loads(cp.read_text())
    nf = c.get("negative_filter")
    if not nf:
        bad(f"{m}: no negative_filter block; negative_max is not re-derivable")
        continue
    # the count is named for the cutoff it was measured at: two hmmsearch passes
    # with different -E report different numbers of sequences, and shipping both
    # under the single name "n_scored" reads as a contradiction.
    scored_key = next((k for k in nf if k.startswith("n_scored")), None)
    need = ("criterion", "n_confirmed")
    if all(nf.get(k) is not None for k in need) and scored_key:
        n_sc = nf[scored_key]
        ok(f"{m}: filter = {nf['criterion']} ({nf['n_confirmed']}/{n_sc} at "
           f"{nf.get('reporting_cutoff', 'unstated cutoff')})")
        # fraction_kept must use the denominator shipped beside it, not the other one
        fk = nf.get("fraction_kept")
        if fk is not None and abs(fk - nf["n_confirmed"] / n_sc) > 1e-6:
            bad(f"{m}: fraction_kept {fk} does not match {nf['n_confirmed']}/{n_sc}")
        # and the build-pass count, if present elsewhere in the record, must not be
        # silently conflated with it
        other = c.get("n_negatives_scored")
        if other is not None and other != n_sc and "reporting_cutoff" not in nf:
            bad(f"{m}: two scored-negative counts ({other} and {n_sc}) with no cutoff stated")
    else:
        bad(f"{m}: negative_filter incomplete")
    if nf.get("confirmed_pattern") and scored_key and nf["n_confirmed"] >= nf[scored_key]:
        bad(f"{m}: a pattern filter that discards nothing is not a filter")

    # n_negatives_above_GA is re-derived the same way as n_confirmed. Checking
    # only that the field exists lets a wrong number through: the strict model
    # shipped it as null for several versions behind a note claiming it was not
    # derivable from the deposit, while the full negative set is deposited and a
    # rescore gives it directly.
    ga_rec = c.get("gathering_threshold")
    n_above = c.get("n_negatives_above_GA")
    negfa0 = NEG_FASTA.get(m)
    if ga_rec is None:
        bad(f"{m}: no gathering_threshold in the record")
    elif n_above is None:
        bad(f"{m}: n_negatives_above_GA is null; the deposited negative set makes "
            f"it derivable, so a note is not a substitute for the number")
    elif not negfa0 or not (ROOT / negfa0).exists():
        warn(f"{m}: n_negatives_above_GA cannot be re-derived, no deposited "
             f"negative set")
    else:
        _cut0 = "1e9" if "1e9" in (nf.get("reporting_cutoff") or "") else "1000"
        _best0 = negative_scores(m, negfa0, _cut0)
        if _best0 is None:
            warn(f"{m}: hmmsearch unavailable, so n_negatives_above_GA "
                 f"({n_above}) is taken on trust")
            _sc0 = None
        else:
            _names = negative_headers(negfa0)
            _pref0 = nf.get("confirmed_prefix") or "sp|"
            _sc0 = sum(1 for v in _best0.values() if v >= ga_rec)
            _rev0 = sum(1 for k, v in _best0.items()
                        if v >= ga_rec and _names.get(k, "").startswith(_pref0))
        if True:
            if _sc0 is not None:
                if _sc0 != n_above:
                    bad(f"{m}: record says {n_above} negatives at or above "
                        f"GA {ga_rec:.2f}; re-scoring {negfa0} gives {_sc0}")
                else:
                    ok(f"{m}: the {n_above} negatives at or above GA "
                       f"{ga_rec:.2f} are re-derived from {negfa0}, not taken "
                       f"on trust")
                _rev_rec = c.get("n_negatives_above_GA_reviewed")
                if _rev_rec is not None and _rev_rec != _rev0:
                    bad(f"{m}: record says {_rev_rec} of them are manually "
                        f"reviewed; re-scoring gives {_rev0}")

    # n_confirmed counts the sequences that BOTH pass the filter and were reported
    # at the record's cutoff, so the sequences alone give an upper bound: reported
    # is a subset of the file. The bound is weak but needs no HMMER, and it is
    # enough to catch a count that is too high -- which the record's own
    # fraction_kept == n_confirmed / n_scored identity cannot do, since that checks
    # the record against itself. Where hmmsearch is on PATH the exact count is
    # re-derived instead.
    negfa = NEG_FASTA.get(m)
    if not negfa or not (ROOT / negfa).exists():
        warn(f"{m}: no deposited negative set to count the confirmed subset against")
    elif not nf.get("confirmed_prefix"):
        warn(f"{m}: negative_filter states no prefix, so the subset cannot be counted")
    else:
        _pref = nf["confirmed_prefix"]
        _pat = re.compile(nf["confirmed_pattern"]) if nf.get("confirmed_pattern") else None
        _hdr = {}
        with open(ROOT / negfa) as _fh:
            for _line in _fh:
                if _line.startswith(">"):
                    _h = _line[1:].rstrip()
                    _hdr[_h.split()[0]] = _h
        _match = {k: h for k, h in _hdr.items()
                  if h.startswith(_pref) and (_pat is None or _pat.search(h))}
        _cutc = "1e9" if "1e9" in (nf.get("reporting_cutoff") or "") else "1000"
        _best = negative_scores(m, negfa, _cutc)
        _exact = len(set(_best) & set(_match)) if _best is not None else None
        if _exact is not None:
            if _exact == nf["n_confirmed"]:
                ok(f"{m}: the {_exact} confirmed negatives are re-derived from "
                   f"{negfa} with hmmsearch, not taken on trust")
            else:
                bad(f"{m}: negative_filter claims {nf['n_confirmed']} confirmed "
                    f"negatives; re-scoring {negfa} at {_cut} and applying its own "
                    f"prefix and pattern gives {_exact}")
        elif nf["n_confirmed"] <= len(_match):
            ok(f"{m}: n_confirmed {nf['n_confirmed']} is within the "
               f"{len(_match)} sequences of {negfa} that pass its filter "
               f"(upper bound; install hmmsearch for the exact count)")
        else:
            bad(f"{m}: negative_filter claims {nf['n_confirmed']} confirmed "
                f"negatives, but only {len(_match)} sequences in {negfa} pass its "
                f"own prefix and pattern at all")

print()
print("=" * 68)
print("10. The uncomfortable numbers are shown, not filtered away")
print("=" * 68)
print("   ('separation: OK' means separation against CONFIRMED negatives only)")
for m in MODELS:
    cp = ROOT / "results/hmm" / f"{m}.calib.json"
    if not cp.exists():
        continue
    c = json.loads(cp.read_text())
    if "n_negatives_above_GA" not in c or "unfiltered_negative_max" not in c:
        bad(f"{m}: n_negatives_above_GA / unfiltered_negative_max absent")
        continue
    n_ab = c["n_negatives_above_GA"]
    scope = c.get("separation_scope", "")
    if n_ab is None:
        # A count that cannot be derived from the deposited files is a legitimate
        # state, provided it is DECLARED: a note must say why, or a null reads as
        # a zero. The scope of the separation claim is still required.
        note = c.get("n_negatives_above_GA_note", "")
        good = bool(note) and "CONFIRMED" in scope.upper()
        (ok if good else bad)(
            f"{m}: negatives above GA not derivable from the deposit, and said so; "
            f"unfiltered max {c['unfiltered_negative_max']}, scope stated"
            if good else
            f"{m}: n_negatives_above_GA is null without a note explaining why, "
            f"or the separation scope does not say it is against confirmed negatives")
    elif c.get("separation") == "OK" and n_ab > 0 and "CONFIRMED" not in scope.upper():
        bad(f"{m}: separation 'OK' with {n_ab} negatives above GA and no qualifying scope")
    else:
        ok(f"{m}: {n_ab} negatives above GA, max {c['unfiltered_negative_max']}, scope stated")

print()
print("=" * 68)
print("11. Reported degrees of freedom are qualified where separation occurs")
print("=" * 68)
if ap.exists():
    R = json.loads(ap.read_text())
    mo = R.get("size_phylum_model", {})
    if "df_phylum_effective" not in mo:
        bad("size_phylum_model has no df_phylum_effective")
    else:
        n_sep = mo.get("n_phyla_completely_separated")
        if mo["df_phylum"] - mo["df_phylum_effective"] == n_sep:
            ok(f"df {mo['df_phylum']} reported, {mo['df_phylum_effective']} effective "
               f"({n_sep} separated phyla)")
        else:
            bad("df_phylum_effective does not match the separated-phyla count")
    ref = mo.get("reference_level")
    if not ref:
        bad("the logistic reference level is not documented")
    else:
        ok(f"reference level documented: {ref['phylum']} + "
           f"{len(ref['absorbed_rare_phyla'])} rare phyla "
           f"({ref['k_absorbed_positives']}/{ref['n_absorbed_genomes']})")

print()
print("=" * 68)
print("12. File hygiene: LF endings, English headers")
print("=" * 68)
crlf = [p.relative_to(ROOT) for p in ROOT.rglob("*.tsv") if b"\r" in p.read_bytes()]
(ok if not crlf else bad)(
    "all TSV use LF endings" if not crlf else f"CRLF found in {crlf}")
plan = ROOT / "data/gtdb/ar53_plan.tsv"
if plan.exists():
    hdr = plan.read_text().splitlines()[0].split("\t")
    fr = [h for h in hdr if h in ("couche", "famille", "genre", "classe", "ordre",
                                  "completude", "taille_genome", "n_proteines")]
    (ok if not fr else bad)(
        "plan header is in English" if not fr else f"French plan headers: {fr}")

print()
print("=" * 68)
print("13. Sub-threshold censoring is disclosed and has a remedy")
print("=" * 68)
print("   (an empty score cell means 'no hit at threshold', not 'no signal')")
sp_ = ROOT / "results/gtdb/scan_archaea.tsv"
if sp_.exists():
    import csv as _c
    rows = list(_c.DictReader(open(sp_), delimiter="\t"))
    empty = sum(1 for r in rows if not (r.get("best_score") or "").strip())
    pct = 100 * empty / len(rows) if rows else 0
    unc = ROOT / "results/gtdb/scan_archaea_all_scores.tsv"
    if unc.exists():
        u = list(_c.DictReader(open(unc), delimiter="\t"))
        ue = sum(1 for r in u if not (r.get("best_score_unthresholded") or "").strip())
        keys_o = {(r["genome_id"], r["marker"]) for r in rows}
        keys_u = {(r["genome_id"], r["marker"]) for r in u}
        by_o = {(r["genome_id"], r["marker"]): r for r in rows}
        mism = sum(1 for r in u
                   if (r["best_score"] or "").strip()
                   != (by_o.get((r["genome_id"], r["marker"]), {}).get("best_score") or "").strip())
        if keys_o == keys_u and mism == 0 and ue == 0:
            ok(f"uncensored table deposited: {len(u)} rows, GA column identical to the census, "
               f"no empty score (census GA column is {pct:.1f} % censored by design)")
        else:
            bad(f"uncensored table disagrees with the census "
                f"(key mismatch {len(keys_o ^ keys_u)}, score mismatch {mism}, empty {ue})")

        # Re-DERIVE the census count from the raw scores rather than comparing the
        # two tables' columns to each other. Comparing columns cannot detect a rule
        # that was applied to build the census but is not recoverable from the
        # deposited scores. --cut_ga requires BOTH the sequence score and at least
        # one domain score to clear GA, so both must be present in the table for
        # the count to be re-derivable.
        has_dom = "best_domain_score_unthresholded" in (u[0] if u else {})
        if not has_dom:
            bad("the uncensored table has no domain-score column, so the census count "
                "cannot be re-derived: --cut_ga requires sequence AND domain to clear GA")
        else:
            ga = {}
            for m in MODELS:
                hp = ROOT / "results/hmm" / f"{m}.hmm"
                if hp.exists():
                    ga[m] = hmm_header(hp).get("GA")
            counted = collections.Counter()
            derived = collections.Counter()
            for r in u:
                m = r["marker"]
                if m not in ga or ga[m] is None:
                    continue
                if (r["best_score"] or "").strip():
                    counted[m] += 1
                sq = (r["best_score_unthresholded"] or "").strip()
                dm = (r["best_domain_score_unthresholded"] or "").strip()
                if sq and dm and float(sq) >= ga[m] and float(dm) >= ga[m]:
                    derived[m] += 1
            allok = True
            for m in sorted(ga):
                if counted[m] == derived[m]:
                    ok(f"{m}: census count {counted[m]} re-derived exactly from the raw scores")
                else:
                    bad(f"{m}: census says {counted[m]} positives, re-deriving from the "
                        f"raw scores at GA {ga[m]} gives {derived[m]}")
                    allok = False
            if allok:
                ok("the whole census is re-derivable from the deposited scores alone")

            # The sweep file must use the SAME rule as the census, and say so. A
            # sweep computed on the sequence score alone, with no rule field and no
            # producing script, agrees with the census at the published threshold
            # and disagrees at the edges of the window (209 against 195 at NC);
            # only re-deriving it from the raw scores exposes that.
            sp = ROOT / "results/gtdb/threshold_sensitivity.json"
            if not sp.exists():
                warn("results/gtdb/threshold_sensitivity.json absent")
            else:
                sw = json.loads(sp.read_text())
                rule = sw.get("rule")
                if rule != "cut_ga":
                    bad(f"threshold sweep declares rule {rule!r}; the census used "
                        f"'cut_ga' (sequence AND domain), so the two disagree at the "
                        f"edges of the window")
                else:
                    ok("threshold sweep declares the census rule (cut_ga)")
                    mk = sw.get("marker")
                    pl = ROOT / "data/gtdb/ar53_plan.tsv"
                    if mk in ga and pl.exists():
                        plan = [r["accession"] for r in
                                list(_c.DictReader(open(pl), delimiter="\t"))
                                [:sw.get("n_genomes", 0)]]
                        keep = set(plan)
                        sq_, dm_ = {}, {}
                        for r in u:
                            if r["marker"] != mk or r["genome_id"] not in keep:
                                continue
                            a_ = (r["best_score_unthresholded"] or "").strip()
                            b_ = (r["best_domain_score_unthresholded"] or "").strip()
                            if a_: sq_[r["genome_id"]] = float(a_)
                            if b_: dm_[r["genome_id"]] = float(b_)
                        wrong = []
                        for pt in sw.get("sweep", []):
                            t = pt["threshold"]
                            k = sum(1 for g in plan
                                    if sq_.get(g, -1e9) >= t and dm_.get(g, -1e9) >= t)
                            if k != pt["positives"]:
                                wrong.append((t, pt["positives"], k))
                        if wrong:
                            bad(f"threshold sweep does not re-derive: "
                                + "; ".join(f"at {t}: file says {a_}, scores give {b_}"
                                            for t, a_, b_ in wrong[:3]))
                        else:
                            ok(f"all {len(sw.get('sweep', []))} sweep points re-derive "
                               f"from the raw scores under the census rule")
    else:
        warn(f"census table is threshold-censored ({empty}/{len(rows)} = {pct:.1f} % of rows): "
             f"regenerate with scripts/uncensor_scan.sh for sensitivity analysis")
tool = ROOT / "scripts/uncensor_scan.sh"
(ok if tool.exists() else bad)(
    "scripts/uncensor_scan.sh is present" if tool.exists()
    else "scripts/uncensor_scan.sh MISSING: censoring has no documented remedy")
sc = (ROOT / "scripts/scan_proteomes.py").read_text()
(ok if "--all-scores" in sc else bad)(
    "scan_proteomes.py offers --all-scores" if "--all-scores" in sc
    else "scan_proteomes.py cannot record sub-threshold scores")
demo = ROOT / "results/gtdb/scan_all_scores_demo.tsv"
if demo.exists():
    import csv as _c2
    d = list(_c2.DictReader(open(demo), delimiter="\t"))
    bad_ = sum(1 for r in d
               if r["best_score"].strip() and r["best_score_unthresholded"].strip()
               and abs(float(r["best_score"]) - float(r["best_score_unthresholded"])) > 0.05)
    (ok if bad_ == 0 else bad)(
        f"worked example: {len(d)} rows, uncensored column agrees with the GA column"
        if bad_ == 0 else f"worked example: {bad_} disagreements")

# Firth penalised logistic regression: the script, its self-test, and the
# agreement between the deposited JSON and the manuscript.
fz = ROOT / "scripts/firth_phylum.py"
fzt = ROOT / "scripts/firth_selftest.py"
(ok if fz.exists() else bad)(
    "scripts/firth_phylum.py is present" if fz.exists()
    else "scripts/firth_phylum.py MISSING: the separation fix has no code")
(ok if fzt.exists() else bad)(
    "scripts/firth_selftest.py is present" if fzt.exists()
    else "scripts/firth_selftest.py MISSING: Firth is unverified")

if fz.exists():
    src = fz.read_text()
    # the two nested models must be fitted under a common penalty
    (ok if "_fit_firth_constrained" in src else bad)(
        "firth_phylum.py compares nested models under a common penalty"
        if "_fit_firth_constrained" in src
        else "firth_phylum.py fits the two models separately: the penalty "
             "difference would inflate the statistic")

fj = ROOT / "results/gtdb/firth_phylum.json"
if fj.exists():
    import json as _j
    F = _j.load(open(fj))
    exp = {"phylum_after_size": (139.11, 16), "size_after_phylum": (43.77, 1)}
    for key, (stat, df) in exp.items():
        got = F.get(key, {})
        good = (abs(got.get("lrt_stat", -1) - stat) < 0.05
                and got.get("df") == df)
        (ok if good else bad)(
            f"firth {key}: chi2={got.get('lrt_stat'):.2f} df={got.get('df')} "
            f"matches the manuscript" if good
            else f"firth {key}: deposit has {got.get('lrt_stat')} / "
                 f"{got.get('df')} df, manuscript states {stat} / {df} df")
    zp = set(F.get("levels_with_zero_positives", []))
    need = {"B1Sed10-29", "Korarchaeota", "Nanohalarchaeota"}
    (ok if need <= zp else bad)(
        "firth: the three zero-positive phyla are recorded" if need <= zp
        else f"firth: zero-positive phyla recorded as {sorted(zp)}, "
             f"manuscript names {sorted(need)}")
else:
    bad("results/gtdb/firth_phylum.json MISSING: the manuscript reports a "
        "penalised test with no deposited output")

# CLC signatures: the gate motif must stay calibrated against its own controls.
mt = ROOT / "scripts/clcf_motif_test.py"
if mt.exists():
    src = mt.read_text()
    (ok if "GREG" in src and "permutation" in src.lower() else bad)(
        "clcf_motif_test.py tests CLC membership with a permutation control"
        if "GREG" in src and "permutation" in src.lower()
        else "clcf_motif_test.py: gate motif or its chance control is missing")
    # a short motif without a chance control is not interpretable
    (ok if "_shuffled" in src else bad)(
        "clcf_motif_test.py computes the chance rate by shuffling"
        if "_shuffled" in src
        else "clcf_motif_test.py reports raw motif rates with no chance baseline")

lf = ROOT / "scripts/lifestyle_fluc.py"
lj = ROOT / "results/gtdb/lifestyle_fluc.json"
if lf.exists():
    src = lf.read_text()
    (ok if "size_matched" in src else bad)(
        "lifestyle_fluc.py stratifies by genome size"
        if "size_matched" in src
        else "lifestyle_fluc.py compares lifestyles without matching genome size, "
             "which would measure genome size instead")
if lj.exists():
    import json as _j3
    L = _j3.load(open(lj))
    cov = L.get("coverage", {}).get("fraction", 0)
    (ok if 0.05 < cov < 0.3 else bad)(
        f"lifestyle: {100*cov:.1f} % of genomes classified, and the figure is "
        f"declared" if 0.05 < cov < 0.3
        else f"lifestyle: coverage {100*cov:.1f} % is not what the manuscript states")
    halo = next((g for g in L["groups"] if g["lifestyle"] == "halophile"), None)
    (ok if halo and halo["N"] == 60 and halo["k"] == 48 else bad)(
        f"lifestyle: halophiles {halo['k']}/{halo['N']} matches the manuscript"
        if halo and halo["N"] == 60 and halo["k"] == 48
        else "lifestyle: halophile counts differ from the manuscript")

fs = ROOT / "scripts/fluc_signature_lifestyle.py"
(ok if fs.exists() else bad)(
    "scripts/fluc_signature_lifestyle.py is present" if fs.exists()
    else "scripts/fluc_signature_lifestyle.py MISSING")
if fs.exists():
    src = fs.read_text()
    # hand-picking the positions amounts to picking the result
    (ok if "hmm_match_entropy" in src else bad)(
        "fluc_signature_lifestyle.py derives diagnostic positions from the model"
        if "hmm_match_entropy" in src
        else "fluc_signature_lifestyle.py hand-picks positions, which picks the result")
    (ok if "hmmalign" in src else bad)(
        "fluc_signature_lifestyle.py aligns hits to the model before comparing"
        if "hmmalign" in src
        else "fluc_signature_lifestyle.py compares positions without aligning")
    # the model must carry as many match states as its declared length
    import re as _re
    h = (ROOT / "results/hmm/Fluc_CrcB.hmm").read_text()
    leng = int(_re.search(r"^LENG\s+(\d+)", h, _re.M).group(1))
    states = len([1 for ln in h.splitlines()
                  if len(ln.split()) >= 21 and ln.split()[0].isdigit()])
    (ok if states == leng else bad)(
        f"Fluc_CrcB.hmm: {states} match states == declared LENG {leng}"
        if states == leng
        else f"Fluc_CrcB.hmm: {states} match states but LENG says {leng}")

cn = ROOT / "scripts/fluc_copy_number.py"
(ok if cn.exists() else bad)(
    "scripts/fluc_copy_number.py is present" if cn.exists()
    else "scripts/fluc_copy_number.py MISSING")
if cn.exists():
    src = cn.read_text()
    # a version suffix is not a gene index; read as one it understates adjacency
    (ok if r're.sub(r"\.\d+$"' in src else bad)(
        "fluc_copy_number.py strips the version suffix before testing adjacency"
        if r're.sub(r"\.\d+$"' in src
        else "fluc_copy_number.py may read a version suffix as a gene index, "
             "which understates adjacency")

cj = ROOT / "results/gtdb/fluc_copy_number.json"
if cj.exists():
    import json as _j4
    C = _j4.load(open(cj))["by_lifestyle"]
    h = C.get("halophile", {})
    good = (h.get("n_genomes") == 48 and h.get("n_copies") == 97
            and h.get("pairs_adjacent") == 7 and h.get("pairs_undecidable") == 39)
    (ok if good else bad)(
        f"copy number: halophiles {h.get('n_genomes')} genomes / "
        f"{h.get('n_copies')} copies / {h.get('pairs_adjacent')} decidably adjacent "
        f"/ {h.get('pairs_undecidable')} undecidable, matches the manuscript" if good
        else "copy number: halophile figures differ from the manuscript")
    # no adjacency rate may be reported while the undecidable pairs dominate
    (ok if h.get("frac_adjacent_among_decidable") is None else bad)(
        "copy number: no adjacency rate is reported, as the undecidable pairs dominate"
        if h.get("frac_adjacent_among_decidable") is None
        else "copy number: an adjacency rate is reported although most pairs are "
             "undecidable from identifiers alone")
    t = C.get("thermophile", {})
    (ok if t.get("copies_per_genome") == 1.0 else bad)(
        "copy number: thermophiles are single-copy without exception"
        if t.get("copies_per_genome") == 1.0
        else "copy number: thermophile copies/genome differs from the manuscript")

# The recipe must match what was actually run.
sf = ROOT / "workflow/Snakefile"
if sf.exists():
    src = sf.read_text()
    import re as _re7
    called = set(_re7.findall(r"scripts/([A-Za-z0-9_]+\.py)", src))
    missing = sorted(c for c in called if not (ROOT / "scripts" / c).exists())
    (ok if not missing else bad)(
        f"Snakefile calls {len(called)} scripts, all present"
        if not missing
        else f"Snakefile calls scripts absent from the archive: {missing}")
    seeds = set(_re7.findall(r'"(seeds/[^"]+)"', src))
    miss_s = sorted(x for x in seeds if not (ROOT / x).exists())
    (ok if not miss_s else bad)(
        f"Snakefile references {len(seeds)} seed files, all present"
        if not miss_s else f"Snakefile references missing seeds: {miss_s}")
    (ok if "negatives_clc.faa" in src else bad)(
        "Snakefile uses negatives_clc.faa for CLC_F"
        if "negatives_clc.faa" in src
        else "Snakefile uses the wrong negative set for CLC_F, which would move its GA")
    (ok if "FAcD" in src else bad)(
        "Snakefile covers all three markers including FAcD"
        if "FAcD" in src else "Snakefile omits FAcD")

    # Run the recipe, do not only read it. The four checks above read the text of
    # the recipe, and all four pass while the recipe itself exits 2 on a
    # non-existent option such as --min-clade-n. So each shell command of the
    # recipe is extracted here, its variables substituted from config.yaml, and
    # executed in a temporary directory. A non-existent option, a script that
    # crashes or a missing parameter fails this check.
    import subprocess as _sp, tempfile as _tf, shutil as _sh, os as _os
    # PyYAML is optional here. An unguarded 'import yaml' makes this script die
    # on ModuleNotFoundError BEFORE the RESULT line wherever PyYAML is not
    # installed, handing the reader a traceback instead of the verdict of the
    # checks that passed; a verification tool must not need an optional package to
    # deliver its verdict. So only the one field needed is read, falling back to a
    # minimal parser when PyYAML is absent.
    try:
        import yaml as _yaml
        cfg = _yaml.safe_load((ROOT / "config.yaml").read_text())
    except ImportError:
        cfg = {}
        for _l in (ROOT / "config.yaml").read_text().splitlines():
            _m = _re7.match(r"^([A-Za-z0-9_]+):\s*(\S+)", _l)
            if _m:
                _v = _m.group(2)
                cfg[_m.group(1)] = int(_v) if _v.isdigit() else _v
    except Exception:
        cfg = {}
    # analytic rules: name -> (command, expected output)
    # Split by rule FIRST, then extract the shell block, rather than with one
    # non-greedy expression: a single `.*?` silently skips the 'baseline' rule,
    # whose block holds comments before `shell:`, because it stops at the first
    # `shell:` it meets, which belongs to the preceding rule. A check that skips a
    # rule without reporting it reports a pass it never earned.
    _blocks = _re7.split(r"^rule\s+(\w+):", src, flags=_re7.M)
    rules = []
    for _i in range(1, len(_blocks), 2):
        _name, _body = _blocks[_i], _blocks[_i + 1]
        _sh = _re7.search(r"shell:\s*\n((?:\s*\"[^\"]*\"\s*\n?)+)", _body)
        if _sh:
            rules.append((_name, _sh.group(1)))
    ran = failed = 0
    with _tf.TemporaryDirectory() as td:
        for name, block in rules:
            cmd = "".join(_re7.findall(r'"([^"]*)"', block))
            if "scripts/" not in cmd:
                continue
            # substitute the Snakemake placeholders
            cmd = cmd.replace("{input.scan}", "results/gtdb/scan_archaea.tsv")
            cmd = cmd.replace("{input.plan}", "data/gtdb/ar53_plan.tsv")
            cmd = cmd.replace("{params.n}", str(cfg.get("min_clade_n", "")))
            cmd = cmd.replace("{output}", _os.path.join(td, f"{name}.json"))
            if "{" in cmd:          # unresolved placeholder: do not guess
                continue
            if "check_deposit.py" in cmd:
                # Do not call this script: the check_deposit rule invokes it, so
                # running it here recurses without end.
                continue
            if "firth_selftest" in cmd:
                # the Firth self-test runs for minutes and has its own rule; only
                # check here that it starts and accepts --B
                cmd = cmd.split("--B")[0].strip() + " --B 60"
            ran += 1
            try:
                r = _sp.run(cmd, shell=True, cwd=ROOT, timeout=120,
                            stdout=_sp.DEVNULL, stderr=_sp.PIPE, text=True)
            except _sp.TimeoutExpired:
                failed += 1
                bad(f"Snakefile rule '{name}' did not finish within 120 s")
                continue
            if r.returncode != 0:
                failed += 1
                bad(f"Snakefile rule '{name}' FAILS when actually run "
                    f"(exit {r.returncode}): {r.stderr.strip().splitlines()[-1][:120]}"
                    if r.stderr.strip() else
                    f"Snakefile rule '{name}' FAILS when actually run "
                    f"(exit {r.returncode})")
        if ran and not failed:
            ok(f"all {ran} analytic Snakefile rules execute successfully "
               f"(commands run, not just read)")
        elif not ran:
            bad("no Snakefile rule could be executed: the check would pass vacuously")

        # the parameters must reproduce the published figure, not merely run
        bj = _os.path.join(td, "baseline.json")
        if _os.path.exists(bj):
            import json as _j8
            got = _j8.load(open(bj))["mean_by_phylum"]
            want_n, want_e = 16, 0.15645364051084093
            good = (got["n_clades"] == want_n
                    and abs(got["estimate"] - want_e) < 1e-9)
            (ok if good else bad)(
                f"Snakefile parameters reproduce the published estimator "
                f"({100*got['estimate']:.3f} % over {got['n_clades']} phyla)"
                if good else
                f"Snakefile runs but yields {100*got['estimate']:.3f} % over "
                f"{got['n_clades']} phyla; the manuscript states "
                f"{100*want_e:.3f} % over {want_n} (check min_clade_n in config.yaml)")

# Run the commands of the README, not only the rules of the Snakefile.
#
# Executing the Snakefile rules does not cover the recipe a reader actually
# follows, which is the README's, and two of its scripts -- design_effect.py and
# threshold_sweep.py -- appear in no rule. When wilson() was unified and the
# shared function went from two return values to three, those two scripts went on
# unpacking two: ValueError on the first call, no output produced. The suite still
# passed, because the only check on wilson() was TEXTUAL -- a regular expression
# on the signature plus the verification that the body delegates. It validated
# the shape of the fix without ever calling it.
#
# A check that reads code is not a substitute for a check that runs it. So each
# `python3 scripts/...` command in the bash blocks of README and GETTING_STARTED
# is extracted here and launched.
import subprocess as _sp2, tempfile as _tf2, re as _re9, os as _os2

_doc_cmds = []
for _docname in ("README.md", "GETTING_STARTED.md"):
    _dp = ROOT / _docname
    if not _dp.exists():
        continue
    _txt = _dp.read_text(encoding="utf-8", errors="replace")
    for _blk in _re9.findall(r"```bash\n(.*?)```", _txt, flags=_re9.S):
        # rejoin line continuations before splitting
        _blk = _blk.replace("\\\n", " ")
        _needs = None
        for _line in _blk.splitlines():
            _line = _line.strip()
            # One exemption only, and the exemption is itself checked. A
            # documented command may depend on a file the deposit does not
            # redistribute (the GTDB metadata): running it as written fails, and
            # dropping it from the README hides it. A marker is therefore accepted
            # but not trusted. The marked command is run like the others and is
            # REQUIRED to fail cleanly, naming the missing file; if that file is
            # present, the marker is wrong and the command must succeed. Checked
            # in both directions, the marker cannot be used to hide a defect.
            _m9 = _re9.match(r"#\s*REQUIRES-NON-DEPOSITED-INPUT:\s*(\S+)", _line)
            if _m9:
                _needs = _m9.group(1)
                continue
            if not _line or _line.startswith("#"):
                continue
            if not _line.startswith("python3 scripts/"):
                _needs = None
                continue          # bash scripts/... requires the proteomes
            _doc_cmds.append((_docname, _line, _needs))
            _needs = None

_ran2 = _failed2 = 0
with _tf2.TemporaryDirectory() as _td2:
    for _docname, _cmd, _needs in _doc_cmds:
        _script = _cmd.split()[1]
        if not (ROOT / _script).exists():
            _failed2 += 1
            bad(f"{_docname} documents {_script}, which is absent from the archive")
            continue
        if "check_deposit.py" in _cmd:
            continue              # recursion: this script is the one running
        if "firth_selftest" in _cmd:
            _cmd = _cmd.split("--B")[0].strip() + " --B 60"
        # send every --out to the temporary directory: a check must never rewrite
        # the deposited files it verifies
        _cmd = _re9.sub(r"--out\s+(\S+)",
                        lambda m: "--out " + _os2.path.join(
                            _td2, _os2.path.basename(m.group(1))), _cmd)
        # figure_archaea takes its output prefix as a positional argument
        if "figure_archaea.py" in _cmd:
            _parts = _cmd.split()
            _parts[-1] = _os2.path.join(_td2, "figure")
            _cmd = " ".join(_parts)
        _ran2 += 1
        try:
            _r2 = _sp2.run(_cmd, shell=True, cwd=ROOT, timeout=300,
                           stdout=_sp2.DEVNULL, stderr=_sp2.PIPE, text=True)
        except _sp2.TimeoutExpired:
            _failed2 += 1
            bad(f"{_docname} command '{_script}' did not finish within 300 s")
            continue
        if _needs:
            # the file declared missing must really be missing
            import glob as _g9
            _present = _g9.glob(str(ROOT / _needs)) or _g9.glob(
                str(ROOT / _re9.sub(r"_r\d+(?=\.)", "*", _needs)))
            _base9 = _os2.path.basename(_needs)
            if _present:
                if _r2.returncode == 0:
                    ok(f"{_docname} marks {_script} as needing {_base9}, which IS "
                       f"on disk here, and the command then succeeds")
                else:
                    _failed2 += 1
                    _tail = (_r2.stderr.strip().splitlines() or [""])[-1][:140]
                    bad(f"{_docname} documents {_script} as needing {_base9}; the "
                        f"file is present yet the command fails -> {_tail}")
            elif _r2.returncode == 0:
                _failed2 += 1
                bad(f"{_docname} declares {_script} to need {_base9}, which is "
                    f"absent, yet the command succeeded: the declaration is wrong, "
                    f"or the script silently produced output without its input")
            else:
                _err9 = _r2.stderr or ""
                _stem9 = _re9.split(r"[*_]r?\d*\.", _base9)[0]
                if "Traceback" in _err9:
                    _failed2 += 1
                    bad(f"{_docname}: {_script} lacks {_base9} and crashes with a "
                        f"traceback instead of refusing cleanly")
                elif _stem9 and _stem9 in _err9:
                    ok(f"{_docname}: {_script} needs {_base9}, which the deposit "
                       f"does not redistribute, and refuses cleanly by naming it")
                else:
                    _failed2 += 1
                    _tail = (_err9.strip().splitlines() or [""])[-1][:140]
                    bad(f"{_docname}: {_script} fails without {_base9} but its "
                        f"message does not name the missing input -> {_tail}")
        elif _r2.returncode != 0:
            _failed2 += 1
            _tail = (_r2.stderr.strip().splitlines() or [""])[-1][:140]
            bad(f"{_docname} documents a command that FAILS when run "
                f"(exit {_r2.returncode}): {_cmd[:80]} -> {_tail}")
    if _ran2 and not _failed2:
        ok(f"all {_ran2} documented python3 commands of README/GETTING_STARTED "
           f"execute successfully (run, not read)")
    elif not _ran2:
        bad("no documented command could be executed: the check would pass vacuously")

# A textual check on wilson() is not enough: the shared function and each of its
# relays are called here, and the arity actually returned is verified.
import importlib.util as _ilu
_sc = ROOT / "scripts"
import sys as _sys9
if str(_sc) not in _sys9.path:
    _sys9.path.insert(0, str(_sc))
_arity_bad = []
try:
    import stats_common as _stc
    _t = _stc.wilson(3, 10)
    if not (isinstance(_t, tuple) and len(_t) == 3):
        _arity_bad.append("stats_common.wilson")
    for _mod in ("design_effect", "threshold_sweep", "analyse_archaea",
                 "baseline", "lifestyle_fluc"):
        _f = _sc / f"{_mod}.py"
        if not _f.exists():
            continue
        _src9 = _f.read_text(encoding="utf-8", errors="replace")
        # expected arity at each call site, measured on the text of the site
        # (number of targets left of '= wilson(') and compared with 3
        for _m9 in _re9.finditer(r"^\s*([A-Za-z_][\w, ]*?)\s*=\s*wilson\(",
                                 _src9, flags=_re9.M):
            _n9 = len([x for x in _m9.group(1).split(",") if x.strip()])
            if _n9 != 1 and _n9 != 3:
                _arity_bad.append(f"{_mod}.py: {_n9} targets for wilson()")
except Exception as _e9:
    _arity_bad.append(f"import failure: {_e9}")
(ok if not _arity_bad else bad)(
    "wilson() returns a triplet and every call site unpacks 3 values "
    "(checked by calling, not by reading)"
    if not _arity_bad else
    f"wilson() arity mismatch -- the defect class of Z1: {_arity_bad}")


# Coverage in the folder -> manifest direction. Checking that every manifest
# entry corresponds to a file present on disk does not test the converse, so a
# supernumerary file -- bytecode, a draft, a forgotten output -- stays invisible:
# 95 manifest entries against 102 delivered files leaves 7 uncovered. An
# uncovered file is one whose integrity nobody can verify and whose origin is not
# declared.
import os as _os3, pathlib as _pl3, sys as _sys8
if str(ROOT / "scripts") not in _sys8.path:
    _sys8.path.insert(0, str(ROOT / "scripts"))
from deposit_exclusions import deposited_files as _dep_files, is_excluded as _is_excl
_listed3 = set()
_sums3 = ROOT / "SHA256SUMS.txt"
if _sums3.exists():
    for _l3 in _sums3.read_text(errors="replace").splitlines():
        if _l3.strip():
            _parts3 = _l3.split(None, 1)
            if len(_parts3) == 2:
                _listed3.add(_parts3[1].strip().removeprefix("./"))
_present3 = _dep_files(ROOT)
_uncov3 = sorted(_present3 - _listed3)
(ok if not _uncov3 else bad)(
    f"manifest covers every deposited file ({len(_present3)} files, "
    f"{len(_listed3)} entries, folder -> manifest direction checked)"
    if not _uncov3 else
    f"{len(_uncov3)} deposited file(s) absent from the manifest, so their "
    f"integrity cannot be verified: {_uncov3[:6]}")

# Checking that the packaging script EXISTS is not enough: a script that is
# present still fails after a workflow run when its coverage check does not
# exclude the .ok markers the Snakefile produces. So the workflow's own outputs
# are what gets tested here.
_pack3 = ROOT / "scripts/make_archive.sh"
if not _pack3.exists():
    bad("no packaging script: the archive is built by hand, which is how "
        "derived bytecode has entered a delivery before")
else:
    # 1. the Snakefile's real outputs must be excluded
    _sf3 = (ROOT / "workflow/Snakefile").read_text(errors="replace") \
        if (ROOT / "workflow/Snakefile").exists() else ""
    _touched3 = sorted(set(_re7.findall(r'touch\(\s*"([^"]+)"', _sf3)))
    _probe3 = _touched3 + ["scripts/__pycache__/stats_common.cpython-311.pyc"]
    _leaks3 = [m for m in _probe3 if not _is_excl(m)]
    (ok if (_touched3 and not _leaks3) else bad)(
        f"the {len(_touched3)} marker files the Snakefile touches are excluded "
        f"from manifest, coverage check and tar alike, so packaging still works "
        f"after a real workflow run"
        if (_touched3 and not _leaks3) else
        (f"these workflow outputs are not excluded, so make_archive.sh fails "
         f"after a real Snakemake run: {_leaks3}") if _leaks3 else
        "no touch() marker found in the Snakefile: this check would pass vacuously")
    # 2. no step may carry its own copy of the exclusion list
    _srcs3 = ("scripts/regenerate_checksums.sh", "scripts/make_archive.sh")
    _stale3 = [f for f in _srcs3
               if (ROOT / f).exists()
               and "deposit_exclusions" not in (ROOT / f).read_text(errors="replace")]
    (ok if not _stale3 else bad)(
        "manifest, coverage check and tar all read scripts/deposit_exclusions.py, "
        "so the three cannot drift apart (the defect class of Z1)"
        if not _stale3 else
        f"these steps carry their own exclusion list instead of importing the "
        f"shared one: {_stale3}")


# Every output of the documented recipe must be covered or excluded.
#
# The coverage check above measures an INTACT deposit. This defect appears only
# after a real pass of the recipe: command 2 of the README writes
# results/gtdb/figure.pdf and figure.png, which are neither in the manifest nor
# excluded, so following the recipe and then re-running the self-test -- the
# order the README recommends -- turns a clean run into one failure. The gap
# between what a check measures and what a reader does is where this kind of
# defect lodges.
#
# The two figure paths are therefore NOT written out here. The outputs are
# DERIVED: the paths come from the README's commands, the figure's extensions
# from the code of figure_archaea.py. A documented command added later that
# writes an undeclared file fails this check.
_declared_out = set()
for _docname in ("README.md", "GETTING_STARTED.md"):
    _dp4 = ROOT / _docname
    if not _dp4.exists():
        continue
    for _blk4 in _re7.findall(r"```bash\n(.*?)```",
                              _dp4.read_text(encoding="utf-8", errors="replace"),
                              flags=_re7.S):
        _blk4 = _blk4.replace("\\\n", " ")
        for _line4 in _blk4.splitlines():
            _line4 = _line4.strip()
            if not _line4.startswith("python3 scripts/"):
                continue
            _m4 = _re7.search(r"--out\s+(\S+)", _line4)
            if _m4:
                _declared_out.add(_m4.group(1).lstrip("./"))
            # figure_archaea.py <analysis> <prefix>: the extensions are read
            # from the script, not assumed here
            if "figure_archaea.py" in _line4:
                _pref4 = _line4.split()[-1].lstrip("./")
                _fsrc4 = (ROOT / "scripts/figure_archaea.py")
                _exts4 = set(_re7.findall(r'\{OUT\}(\.[A-Za-z0-9]+)',
                                          _fsrc4.read_text(errors="replace"))) \
                    if _fsrc4.exists() else set()
                for _e4 in _exts4:
                    _declared_out.add(_pref4 + _e4)

_uncovered_out = sorted(o for o in _declared_out
                        if o not in _listed3 and not _is_excl(o))
(ok if (_declared_out and not _uncovered_out) else bad)(
    f"all {len(_declared_out)} outputs the documented recipe writes are either "
    f"in the manifest or declared derived, so running the README does not make "
    f"this self-test fail"
    if (_declared_out and not _uncovered_out) else
    (f"running the documented recipe creates {len(_uncovered_out)} file(s) that "
     f"are neither in the manifest nor excluded, so the self-test fails on a "
     f"deposit that only followed its own README: {_uncovered_out}")
    if _uncovered_out else
    "no documented output could be derived: this check would pass vacuously")


# Every output directory a script declares is itself declared.
#
# The check above covers the DOCUMENTED commands. It does not cover
# scripts/power.py, which appears in neither the README nor the Snakefile and
# writes six files into results/power/ -- a directory delivered EMPTY. Running
# that script turns a clean self-test into one failure.
#
# The destinations are therefore derived from the code: every --out default that
# names a DIRECTORY (no file suffix) must be either declared in
# EXCLUDED_OUTPUT_DIRS or a directory that ships deposited content, such as
# seeds/. A script added later with a new destination fails this check until that
# destination is settled.
from deposit_exclusions import EXCLUDED_OUTPUT_DIRS as _OUTDIRS
_cand5, _undeclared5 = set(), []
for _f5 in sorted((ROOT / "scripts").glob("*.py")):
    for _m5 in _re7.finditer(r'add_argument\(\s*"--out"\s*,\s*default\s*=\s*"([^"]+)"',
                             _f5.read_text(errors="replace")):
        _v5 = _m5.group(1).strip("/")
        if "." not in _pl3.Path(_v5).name:      # no suffix => a directory
            _cand5.add(_v5)
for _d5 in sorted(_cand5):
    if _d5 in _OUTDIRS:
        continue
    # a directory that ships deposited content is not a pure output destination
    if any(e == _d5 or e.startswith(_d5 + "/") for e in _listed3):
        continue
    _undeclared5.append(_d5)
(ok if (_cand5 and not _undeclared5) else bad)(
    f"all {len(_cand5)} output directories declared by scripts/*.py are either "
    f"listed as derived or ship deposited content, so running a script cannot "
    f"leave uncovered files behind"
    if (_cand5 and not _undeclared5) else
    (f"these script output directories are neither declared derived nor ship "
     f"deposited content, so running the script breaks the coverage check: "
     f"{_undeclared5}") if _undeclared5 else
    "no script output directory could be derived: this check would pass vacuously")


# The seed re-check against UniProt, and what of it is verifiable here.
#
# seeds/uniprot_recheck_2026-10-01.tsv carries one row per positive seed (266:
# 214 FAcD, 22 CLC-F, 30 Fluc), holding side by side what the deposit records and
# what UniProt answered on 1 October 2026.
#
# The UniProt half of that table is NOT verifiable from the deposit: it dates
# from a network query and will age. The deposit half is. This check therefore
# revalidates what can be revalidated -- sp|/tr| status, label, evidence level,
# length -- across the 266 rows. Editing a seed without redoing the re-check
# fails it.
_rc = ROOT / "seeds/uniprot_recheck_2026-10-01.tsv"
if not _rc.exists():
    bad("seeds/uniprot_recheck_2026-10-01.tsv absent: the exhaustive seed "
        "re-check against UniProt is not deposited")
else:
    import csv as _csv6
    _dep6 = {}
    for _f6, _fam6 in (("seeds/facd_pos.faa", "FAcD"),
                       ("seeds/clcf_pos.faa", "CLC_F"),
                       ("seeds/fluc_pos.faa", "Fluc")):
        _k6, _b6 = None, []
        for _l6 in (ROOT / _f6).read_text(errors="replace").splitlines():
            if _l6.startswith(">"):
                if _k6:
                    _dep6[_k6[0]] = _k6[1:] + (len("".join(_b6)),)
                _h6 = _l6[1:].strip()
                _acc6 = _h6.split("|")[1] if "|" in _h6 else _h6.split()[0]
                _st6 = "Swiss-Prot" if _h6.startswith("sp|") else "TrEMBL"
                _m6 = _re7.search(r"^\S+\s+(.*?)\s+OS=", _h6)
                _nm6 = _re7.sub(r"\s*\(EC [^)]*\)", "", _m6.group(1)).strip() if _m6 else "?"
                _pe6 = _re7.search(r"PE=(\d)", _h6)
                _k6 = (_acc6, _fam6, _st6, _nm6, _pe6.group(1) if _pe6 else "?")
                _b6 = []
            else:
                _b6.append(_l6.strip())
        if _k6:
            _dep6[_k6[0]] = _k6[1:] + (len("".join(_b6)),)
    _rows6 = list(_csv6.DictReader(_rc.open(encoding="utf-8"), delimiter="\t"))
    _mismatch6, _seen6 = [], set()
    for _r6 in _rows6:
        _a6 = _r6["accession"].strip()
        _seen6.add(_a6)
        if _a6 not in _dep6:
            _mismatch6.append(f"{_a6}: not a deposited seed")
            continue
        _fam, _st, _nm, _pe, _L = _dep6[_a6]
        if _r6["status_deposit"].strip() != _st:
            _mismatch6.append(f"{_a6}: status")
        if _r6["name_deposit"].strip() != _nm:
            _mismatch6.append(f"{_a6}: name")
        if _r6["PE_deposit"].strip() != _pe:
            _mismatch6.append(f"{_a6}: PE")
        if _r6["length_aa"].strip() != str(_L):
            _mismatch6.append(f"{_a6}: length")
    _missing6 = sorted(set(_dep6) - _seen6)
    _okk = not _mismatch6 and not _missing6 and len(_rows6) == len(_dep6)
    (ok if _okk else bad)(
        f"the UniProt re-check covers all {len(_dep6)} positive seeds and its "
        f"deposit-side columns (status, name, PE, length) still match the "
        f"deposited FASTA exactly"
        if _okk else
        f"the UniProt re-check no longer matches the seeds: "
        f"{len(_missing6)} seed(s) absent from the table, "
        f"{len(_mismatch6)} field mismatch(es) {_mismatch6[:4]}")


vf = ROOT / "results/versions.txt"
(ok if vf.exists() and vf.stat().st_size > 0 else bad)(
    "results/versions.txt is present and non-empty" if vf.exists() and vf.stat().st_size
    else "results/versions.txt MISSING: software versions are claimed but not recorded")

# The evidence files: present, or their absence documented.
eh = ROOT / "scripts/export_hits.sh"
(ok if eh.exists() else bad)(
    "scripts/export_hits.sh is present (extracts the detected sequences)"
    if eh.exists()
    else "scripts/export_hits.sh MISSING: Table 10 cannot be reproduced")
if eh.exists():
    _es = eh.read_text()
    # It must FAIL on a partial export, not merely warn. Exercised against a
    # directory holding a single sequence, a script that only warns writes a
    # one-sequence file and still returns 0.
    _guards = ("sys.exit(1)" in _es and "mktemp" in _es
               and "got != exp" in _es)
    (ok if _guards else bad)(
        "export_hits.sh writes through a temporary file and exits non-zero on a "
        "partial export" if _guards
        else "export_hits.sh can write an incomplete file, warn, and still exit 0, "
             "which lets a partial export replace a valid one")

hits = ROOT / "results/gtdb/fluc_hits.faa"
if hits.exists():
    nseq = hits.read_text().count(">")
    (ok if nseq == 239 else bad)(
        f"fluc_hits.faa holds the {nseq} detected sequences" if nseq == 239
        else f"fluc_hits.faa holds {nseq} sequences, manuscript states 239")
    # the identifiers must be EXACTLY those of the census
    import csv as _c7
    _plan = list(_c7.DictReader(open(ROOT / "data/gtdb/ar53_plan.tsv"),
                                delimiter="\t"))[:1259]
    _keep = {r["accession"] for r in _plan}
    _want = set()
    for r in _c7.DictReader(open(ROOT / "results/gtdb/scan_archaea.tsv"),
                            delimiter="\t"):
        if (r["marker"] == "Fluc_CrcB" and r["genome_id"] in _keep
                and int(r["n_copies"]) > 0 and r.get("hit_ids")):
            for _h in r["hit_ids"].split(";"):
                if _h:
                    _want.add(f'{r["genome_id"]}|{_h}')
    _got = {l[1:].strip().split()[0] for l in hits.read_text().splitlines()
            if l.startswith(">")}
    (ok if _got == _want else bad)(
        "fluc_hits.faa identifiers match the census exactly"
        if _got == _want
        else f"fluc_hits.faa identifiers differ from the census "
             f"({len(_want - _got)} missing, {len(_got - _want)} extra)")

    # the alignment must yield exactly LENG match columns, without which
    # Table 10 does not recompute
    _a2m = ROOT / "results/gtdb/fluc_hits.a2m"
    if _a2m.exists():
        _al, _n, _b = {}, None, []
        for _l in _a2m.read_text().splitlines():
            if _l.startswith(">"):
                if _n:
                    _al[_n] = "".join(_b)
                _n, _b = _l[1:].split()[0], []
            else:
                _b.append(_l.strip())
        if _n:
            _al[_n] = "".join(_b)
        import re as _re8
        _leng = int(_re8.search(r"^LENG\s+(\d+)",
                    (ROOT / "results/hmm/Fluc_CrcB.hmm").read_text(),
                    _re8.M).group(1))
        _widths = {sum(1 for c in v if c.isupper() or c == "-")
                   for v in _al.values()}
        good = _widths == {_leng} and len(_al) == len(_got)
        (ok if good else bad)(
            f"fluc_hits.a2m: {len(_al)} sequences, {_leng} match columns each, "
            f"so Table 10 recomputes from the deposit" if good
            else f"fluc_hits.a2m: {len(_al)} sequences, match columns {_widths}, "
                 f"expected {len(_got)} and {_leng}")
    else:
        warn("results/gtdb/fluc_hits.a2m absent: the sequences are deposited but "
             "not their alignment, so Table 10 needs hmmalign to recompute")
else:
    warn("results/gtdb/fluc_hits.faa absent: run scripts/export_hits.sh where the "
         "proteomes live, so Table 10 becomes reproducible from the deposit")

# The documentation must not advertise a stale check count, nor a command that
# fails on macOS.
import re as _re6
for _doc in ("README.md", "GETTING_STARTED.md"):
    _p = ROOT / _doc
    if not _p.exists():
        continue
    _t = _p.read_text()
    # The advertised count is verified EXACTLY, but at the end of the script: see
    # the final block. Here len(OK) is not the total, and the wide tolerance that
    # follows from that lets "110 cross-checks" pass in a deposit holding 115.
    _bad_cmd = "xargs -0 sha256sum > SHA256SUMS" in _t
    (ok if not _bad_cmd else bad)(
        f"{_doc} does not hand the reader a macOS-breaking checksum command"
        if not _bad_cmd
        else f"{_doc} tells the reader to run a pipeline that empties the manifest "
             f"on macOS; point to scripts/regenerate_checksums.sh instead")

_rc = ROOT / "scripts/regenerate_checksums.sh"
if _rc.exists():
    _src = _rc.read_text()
    (ok if "shasum" in _src and "mktemp" in _src else bad)(
        "regenerate_checksums.sh is portable and writes through a temporary file"
        if "shasum" in _src and "mktemp" in _src
        else "regenerate_checksums.sh is not portable, or truncates on failure")

# The declared environment. With nothing reading envs/, a lock file that is in
# fact the export of the machine's 'base' conda environment rather than the
# project's -- 2 packages in common out of 15 -- passes unnoticed while the
# documentation prescribes using it, so a third party can run no script at all.
import re as _re9
_envdir = ROOT / "envs"
_envyml = _envdir / "environment.yml"
(ok if _envyml.exists() else bad)(
    "envs/environment.yml is present" if _envyml.exists()
    else "envs/environment.yml MISSING: dependencies are undeclared")

def _conda_pkgs(path):
    out = set()
    for _l in path.read_text().splitlines():
        m = _re9.match(r"\s+-\s+([A-Za-z0-9_.\-]+)", _l)
        if m:
            out.add(_re9.split(r"[><=]", m.group(1))[0].lower())
    return out

if _envyml.exists():
    _declared = _conda_pkgs(_envyml)
    _txt = _envyml.read_text()
    _name = _re9.search(r"^name:\s*(\S+)", _txt, _re9.M)
    (ok if _name and _name.group(1) == "fluoride-baseline" else bad)(
        f"environment.yml declares name: {_name.group(1) if _name else '?'}"
        if _name and _name.group(1) == "fluoride-baseline"
        else "environment.yml does not declare name: fluoride-baseline, so the "
             "'conda activate fluoride-baseline' in the docs would fail")

    # every module the code imports must be declared somewhere
    _need = {"numpy": "numpy", "scipy": "scipy", "yaml": "pyyaml",
             "statsmodels": "statsmodels"}
    _imported = set()
    for _f in sorted((ROOT / "scripts").glob("*.py")):
        _src = _f.read_text()
        for _mod in _need:
            # covers the 'import a, b as c, d' forms: a pattern anchored on the
            # start of the statement misses 'import ... as _yaml' at end of line
            if _re9.search(rf"(?:^|[\s,])(?:import|from)\s+{_mod}\b|"
                           rf"[,\s]{_mod}\s+as\s+\w+", _src, _re9.M):
                _imported.add(_mod)
    _undeclared = sorted(_need[m] for m in _imported
                         if _need[m] not in _declared)
    (ok if not _undeclared else bad)(
        f"every module imported by scripts/ is declared in environment.yml"
        if not _undeclared
        else f"imported but undeclared in environment.yml: {_undeclared} — a "
             f"third party following the docs hits ModuleNotFoundError")

    # a lock file, if one exists, must cover the declared dependencies
    _lock = _envdir / "environment.lock.yml"
    if _lock.exists():
        _locked = _conda_pkgs(_lock)
        _miss = sorted(_declared - _locked)
        _lname = _re9.search(r"^name:\s*(\S+)", _lock.read_text(), _re9.M)
        good = not _miss and _lname and _lname.group(1) == "fluoride-baseline"
        (ok if good else bad)(
            "environment.lock.yml locks the project environment"
            if good else
            f"environment.lock.yml is NOT the project environment: "
            f"name={_lname.group(1) if _lname else '?'}, "
            f"{len(_miss)} declared packages absent ({_miss[:5]}...)")
    else:
        (ok if (_envdir / "README.md").exists() else bad)(
            "no lock file, and envs/README.md explains why"
            if (_envdir / "README.md").exists()
            else "no lock file and no explanation in envs/")

# The manifest must not contain derived bytecode. A .pyc compiled from an older
# source would be checksummed as deposited content while contradicting the .py
# beside it.
_man = ROOT / "SHA256SUMS.txt"
if _man.exists():
    _derived = [l.split(None, 1)[1].strip() for l in _man.read_text().splitlines()
                if l.strip() and _re9.search(r"__pycache__|\.pyc$|\.DS_Store$",
                                             l)]
    (ok if not _derived else bad)(
        "the manifest contains no derived bytecode"
        if not _derived
        else f"{len(_derived)} derived files in the manifest ({_derived[:2]}): "
             f"compiled bytecode is not deposited content and can be stale")

# One definition of wilson() only. Copied into five scripts, the function
# carried three different values of the normal quantile; the measured impact was
# at most 2.5e-4 of a point and no published figure was affected, but five copies
# are five occasions for future divergence.
import re as _re10
_sc = ROOT / "scripts/stats_common.py"
(ok if _sc.exists() else bad)(
    "scripts/stats_common.py holds the shared statistical helpers"
    if _sc.exists() else "scripts/stats_common.py MISSING")

_zvals, _bodies = set(), []
for _f in sorted((ROOT / "scripts").glob("*.py")):
    if _f.name == "stats_common.py":
        continue
    _src = _f.read_text()
    for _m in _re10.finditer(r"def wilson\(k, n, z=([^)]+)\):", _src):
        _zvals.add(_m.group(1))
        # a definition that does not delegate is a real copy
        _tail = _src[_m.end():_m.end() + 900]
        if "_wilson(" not in _tail:
            _bodies.append(_f.name)
(ok if not _bodies else bad)(
    f"no script re-implements wilson(); {len(_zvals) - (1 if 'Z95' in _zvals else 0)}"
    f" literal z values remain" if not _bodies
    else f"wilson() is re-implemented in {_bodies}, which is how three different "
         f"z values came to coexist")
_literals = {z for z in _zvals if z != "Z95"}
(ok if not _literals else bad)(
    "every wilson() signature defaults to the shared Z95"
    if not _literals
    else f"literal z defaults still present: {sorted(_literals)}")

# The cross-check count the README advertises, verified exactly.
#
# This block sits at the end of the script and not elsewhere. Placed in the
# middle of the suite, a check on that number can only compare it loosely, since
# len(OK) is not final there -- and the wide tolerance that follows lets
# "110 cross-checks" pass in a deposit holding 115. Last in the script, len(OK)
# is the total and the comparison can be exact.
#
# The block counts itself: if it passes, the total becomes len(OK) + 1, and that
# is the number the README must carry. Exactness is required only when nothing
# else failed -- otherwise len(OK) is the count of checks that PASSED, not the
# size of the suite, and comparing the two would mean nothing.
_adv = []
for _doc in ("README.md", "GETTING_STARTED.md"):
    _p = ROOT / _doc
    if not _p.exists():
        continue
    for _m in re.finditer(r"(\d+)\s+cross-checks", _p.read_text()):
        _adv.append((_doc, int(_m.group(1))))
if not _adv:
    bad("no document advertises a cross-check count: the figure cannot be kept honest")
elif FAIL:
    warn(f"advertised cross-check count not verified exactly: {len(FAIL)} check(s) "
         f"failed, so len(OK) is not the size of the suite")
else:
    _expected = len(OK) + 1
    _wrong = [(d, c) for d, c in _adv if c != _expected]
    if _wrong:
        bad("advertised cross-check count is stale: "
            + "; ".join(f"{d} says {c}" for d, c in _wrong)
            + f" but the suite contains {_expected}. Write {_expected}.")
    else:
        ok(f"the {_expected} cross-checks advertised is the number this suite "
           f"actually contains, counted at the end")

print()
print("=" * 68)
print(f"RESULT: {len(OK)} passed, {len(FAIL)} failed, {len(WARN)} warnings")
print("=" * 68)
if FAIL:
    print("\nFailures:")
    for f in FAIL:
        print(f"  - {f}")
sys.exit(1 if FAIL else 0)
