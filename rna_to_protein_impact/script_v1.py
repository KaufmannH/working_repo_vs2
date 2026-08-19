#!/usr/bin/env python3
import os
import re
import math
import time
import argparse
import concurrent.futures as cf
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict

import pandas as pd
import requests
import matplotlib.pyplot as plt


# --- UniProt endpoints (stable) ---
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_JSON = "https://rest.uniprot.org/uniprotkb/{acc}.json"

TAXON_MOUSE = "10090"

# --- heuristics ---
PHOSPHO_IN_DESC = re.compile(r"\bphospho", re.IGNORECASE)
HALF_LIFE_REGEX = re.compile(
    r"Half-life:\s*([0-9]+(?:\.[0-9]+)?)\s*(h|hour|hours|min|mins|minutes|d|day|days)\b",
    re.IGNORECASE,
)
MONOMER_REGEX = re.compile(r"\bmonomer\b", re.IGNORECASE)
COMPLEX_CUES = re.compile(r"\b(complex|hetero|homo|subunit|component of|part of|forms a)\b", re.IGNORECASE)


@dataclass
class Weights:
    w_mrna_hl: float = 0.55
    w_prot_hl: float = 0.45
    w_complex_penalty: float = 0.25
    w_ptm_penalty: float = 0.15
    activity_emphasis: float = 1.0
    missing_penalty_per_feature: float = 0.08


def logistic(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1 / (1 + z)
    z = math.exp(x)
    return z / (1 + z)


def zscore(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    mu = s.mean(skipna=True)
    sd = s.std(skipna=True)
    if sd == 0 or pd.isna(sd):
        return (s - mu) * 0.0
    return (s - mu) / sd


def parse_half_life_hours_from_biophys(text: str) -> Optional[float]:
    m = HALF_LIFE_REGEX.search(text or "")
    if not m:
        return None
    value = float(m.group(1))
    unit = m.group(2).lower()
    if unit.startswith("min"):
        return value / 60.0
    if unit.startswith("d") or unit.startswith("day"):
        return value * 24.0
    return value


def extract_uniprot_fields(entry: dict) -> Dict:
    # accession
    acc = entry.get("primaryAccession", "")

    # protein name
    protein_name = ""
    try:
        pdsc = entry.get("proteinDescription") or {}
        rec = pdsc.get("recommendedName") or {}
        full = (rec.get("fullName") or {}).get("value")
        if isinstance(full, str):
            protein_name = full
    except Exception:
        pass

    # gene primary
    gene_primary = ""
    genes = entry.get("genes") or []
    if genes and isinstance(genes, list):
        g0 = genes[0] or {}
        gname = (g0.get("geneName") or {}).get("value")
        if isinstance(gname, str):
            gene_primary = gname

    # reviewed
    reviewed = entry.get("entryType", "")
    is_reviewed = True if reviewed == "UniProtKB reviewed (Swiss-Prot)" else False

    # SUBUNIT comment -> in_complex heuristic
    subunit_text = ""
    comments = entry.get("comments") or []
    if isinstance(comments, list):
        for c in comments:
            if (c or {}).get("commentType") == "SUBUNIT":
                texts = (c or {}).get("texts") or []
                vals = []
                for t in texts:
                    v = (t or {}).get("value")
                    if isinstance(v, str):
                        vals.append(v.strip())
                subunit_text = " ".join(vals).strip()
                break

    in_complex = None
    if subunit_text:
        if MONOMER_REGEX.search(subunit_text):
            in_complex = False
        elif COMPLEX_CUES.search(subunit_text):
            in_complex = True

    # phospho sites: count "Modified residue" features with phospho in description
    phospho_sites = 0
    feats = entry.get("features") or []
    if isinstance(feats, list):
        for f in feats:
            if (f or {}).get("type") == "Modified residue":
                desc = (f or {}).get("description", "")
                if isinstance(desc, str) and PHOSPHO_IN_DESC.search(desc):
                    phospho_sites += 1

    # protein half-life: from BIOPHYSICOCHEMICAL_PROPERTIES comment text (rare)
    prot_half_life_h = None
    biophys_text = ""
    if isinstance(comments, list):
        for c in comments:
            if (c or {}).get("commentType") == "BIOPHYSICOCHEMICAL_PROPERTIES":
                texts = (c or {}).get("texts") or []
                vals = []
                for t in texts:
                    v = (t or {}).get("value")
                    if isinstance(v, str):
                        vals.append(v.strip())
                biophys_text = " ".join(vals).strip()
                break
    if biophys_text:
        prot_half_life_h = parse_half_life_hours_from_biophys(biophys_text)

    return {
        "accession": acc,
        "gene_symbol": gene_primary,
        "protein_name": protein_name,
        "reviewed": is_reviewed,
        "subunit_comment": subunit_text,
        "in_complex": in_complex,
        "phospho_site_count": phospho_sites,
        "protein_half_life_hours": prot_half_life_h,
    }


def fetch_all_mouse_accessions(session: requests.Session, reviewed_only: bool) -> List[str]:
    # pull ONLY accessions via TSV (stable + tiny)
    query = f"organism_id:{TAXON_MOUSE}"
    if reviewed_only:
        query += " AND reviewed:true"

    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession",
        "size": 500,
    }

    url = UNIPROT_SEARCH
    accs: List[str] = []

    while True:
        r = session.get(url, params=params, timeout=60)
        r.raise_for_status()
        lines = r.text.splitlines()
        # first line header: "Accession"
        for line in lines[1:]:
            a = line.strip()
            if a:
                accs.append(a)

        link = r.headers.get("Link", "")
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        if not m:
            break
        url = m.group(1)
        params = None

    return sorted(set(accs))


def fetch_entry_json(session: requests.Session, acc: str, retries: int = 3) -> Optional[dict]:
    url = UNIPROT_ENTRY_JSON.format(acc=acc)
    for i in range(retries):
        try:
            r = session.get(url, timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except requests.HTTPError as e:
            # back off on 429/5xx
            code = getattr(e.response, "status_code", None)
            if code in (429, 500, 502, 503, 504):
                time.sleep(1.5 * (i + 1))
                continue
            raise
        except Exception:
            time.sleep(1.0 * (i + 1))
    return None


def build_scores(df: pd.DataFrame, weights: Weights) -> pd.DataFrame:
    df = df.copy()

    # placeholders unless you merge real half-life tables later
    if "mrna_half_life_hours" not in df.columns:
        df["mrna_half_life_hours"] = pd.NA

    df["mrna_hl_z"] = zscore(df["mrna_half_life_hours"])
    df["prot_hl_z"] = zscore(df["protein_half_life_hours"])

    df["complex_penalty"] = df["in_complex"].fillna(False).astype(bool).astype(float) * weights.w_complex_penalty
    df["ptm_penalty"] = df["phospho_site_count"].fillna(0).astype(float).apply(lambda x: math.log1p(x)) * weights.w_ptm_penalty

    miss = (
        df["mrna_half_life_hours"].isna().astype(int)
        + df["protein_half_life_hours"].isna().astype(int)
        + df["in_complex"].isna().astype(int)
    )
    df["missing_count"] = miss
    df["confidence"] = (1.0 - weights.missing_penalty_per_feature * miss).clip(0.0, 1.0)

    abundance_logit = (
        weights.w_mrna_hl * df["mrna_hl_z"].fillna(0.0)
        + weights.w_prot_hl * df["prot_hl_z"].fillna(0.0)
        - df["complex_penalty"]
        - df["ptm_penalty"]
    )
    df["abundance_score_raw"] = abundance_logit.apply(logistic)
    df["abundance_score"] = df["abundance_score_raw"] * df["confidence"]

    # activity penalty: PTM burden + complex membership
    activity_reg_index = (
        (df["phospho_site_count"].fillna(0).astype(float).apply(lambda x: math.log1p(x)) / math.log(1 + 20))
        + df["in_complex"].fillna(False).astype(bool).astype(float) * 0.35
    ).clip(0.0, 1.0)

    df["activity_reg_index"] = activity_reg_index
    df["activity_score"] = (df["abundance_score_raw"] * (1.0 - weights.activity_emphasis * df["activity_reg_index"])) * df["confidence"]

    return df


def make_plots(df: pd.DataFrame) -> None:
    # missingness
    missing = {
        "mRNA half-life": df["mrna_half_life_hours"].isna().mean(),
        "Protein half-life": df["protein_half_life_hours"].isna().mean(),
        "Complex": df["in_complex"].isna().mean(),
        "Phospho sites": df["phospho_site_count"].isna().mean(),
    }
    plt.figure()
    plt.bar(list(missing.keys()), list(missing.values()))
    plt.ylabel("Fraction missing")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig("missingness.png", dpi=200)
    plt.close()

    plt.figure()
    plt.hist(df["abundance_score"].dropna(), bins=50)
    plt.xlabel("abundance_score")
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig("abundance_score_hist.png", dpi=200)
    plt.close()

    plt.figure()
    plt.hist(df["activity_score"].dropna(), bins=50)
    plt.xlabel("activity_score")
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig("activity_score_hist.png", dpi=200)
    plt.close()


def write_report(df: pd.DataFrame, weights: Weights, reviewed_only: bool) -> None:
    lines = []
    lines.append("# Mouse RNA → protein impact table\n\n")
    lines.append(f"- rows: {len(df)}\n")
    lines.append(f"- reviewed_only: {reviewed_only}\n")
    lines.append(f"- complex annotated (heuristic): {int(df['in_complex'].notna().sum())}\n")
    lines.append(f"- phospho sites >0: {int((df['phospho_site_count'] > 0).sum())}\n")
    lines.append(f"- protein half-life present: {int(df['protein_half_life_hours'].notna().sum())}\n")
    lines.append(f"- mRNA half-life present: {int(df['mrna_half_life_hours'].notna().sum())}\n\n")

    lines.append("## Weights\n")
    for k, v in weights.__dict__.items():
        lines.append(f"- {k}: {v}\n")

    lines.append("\n## Files\n")
    lines.append("- mouse_rna2protein_impact.tsv\n")
    lines.append("- missingness.png\n")
    lines.append("- abundance_score_hist.png\n")
    lines.append("- activity_score_hist.png\n")

    with open("report.md", "w", encoding="utf-8") as f:
        f.write("".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reviewed-only", action="store_true", help="Swiss-Prot only (smaller, cleaner)")
    ap.add_argument("--threads", type=int, default=8, help="Parallel entry fetch")
    ap.add_argument("--activity-emphasis", type=float, default=1.0, help="1.0 activity focus, 0.0 abundance-only")
    ap.add_argument("--sleep", type=float, default=0.0, help="Extra delay between requests (seconds)")
    args = ap.parse_args()

    weights = Weights(activity_emphasis=args.activity_emphasis)

    session = requests.Session()
    session.headers.update({"User-Agent": "mouse_rna2protein_impact/2.0"})

    print("Fetching mouse accessions from UniProt...")
    accs = fetch_all_mouse_accessions(session, reviewed_only=args.reviewed_only)
    print(f"Accessions: {len(accs)}")

    # Fetch entry JSONs
    print("Downloading entry JSONs (UniProt)...")
    rows = []
    with cf.ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futs = {ex.submit(fetch_entry_json, session, a): a for a in accs}
        for i, fut in enumerate(cf.as_completed(futs), start=1):
            a = futs[fut]
            try:
                entry = fut.result()
            except Exception as e:
                print(f"Failed {a}: {e}")
                entry = None
            if entry:
                rows.append(extract_uniprot_fields(entry))
            if args.sleep and args.sleep > 0:
                time.sleep(args.sleep)
            if i % 1000 == 0:
                print(f"... {i}/{len(accs)}")

    df = pd.DataFrame(rows)

    # Score
    df = build_scores(df, weights)

    # Write next to script (current working directory)
    df.to_csv("mouse_rna2protein_impact.tsv", sep="\t", index=False)

    make_plots(df)
    write_report(df, weights, reviewed_only=args.reviewed_only)

    print("Done. Wrote files next to script.")


if __name__ == "__main__":
    main()