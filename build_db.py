#!/usr/bin/env python3
"""
Build the protein-protein interaction database from public sources.

Downloads and merges:
  - STRING v12.0 physical links (human, score >= 400)
  - BioGRID (human, physical interactions)
  - HuRI (Human Reference Interactome)

Outputs: data/interactions.csv.gz
"""

import csv
import gzip
import io
import os
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from urllib.request import urlretrieve

DATA_DIR = Path(__file__).parent / "data"
CACHE_DIR = DATA_DIR / "raw"

# --- Download URLs ---
STRING_PHYSICAL = "https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/9606.protein.physical.links.detailed.v12.0.txt.gz"
STRING_ALIASES = "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
STRING_INFO = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
BIOGRID_URL = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip"
HURI_URL = "http://www.interactome-atlas.org/data/HuRI.tsv"

STRING_MIN_SCORE = 400  # medium confidence


def download(url: str, dest: Path) -> Path:
    """Download a file if not already cached."""
    if dest.exists():
        print(f"  [cached] {dest.name}")
        return dest
    print(f"  Downloading {url} ...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    urlretrieve(url, dest)
    print(f"  -> {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def build_id_mapping() -> tuple[dict, dict]:
    """
    Build mappings from STRING Ensembl protein IDs to gene symbols and UniProt IDs.

    Returns:
        gene_map: {ensp_id: preferred_gene_symbol}
        uniprot_map: {ensp_id: uniprot_accession}
    """
    print("\n=== Building ID mappings from STRING aliases ===")
    aliases_path = download(STRING_ALIASES, CACHE_DIR / "9606.protein.aliases.txt.gz")
    info_path = download(STRING_INFO, CACHE_DIR / "9606.protein.info.txt.gz")

    # Get preferred names from protein.info (most reliable for gene symbols)
    gene_map: dict[str, str] = {}
    with gzip.open(info_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ensp = row["#string_protein_id"].replace("9606.", "")
            preferred = row["preferred_name"]
            if preferred:
                gene_map[ensp] = preferred

    # Get UniProt mappings from aliases
    uniprot_map: dict[str, str] = {}
    with gzip.open(aliases_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ensp = row["#string_protein_id"].replace("9606.", "")
            alias = row["alias"]
            source = row["source"]
            # Prefer Swiss-Prot (reviewed) UniProt IDs
            if "UniProt_AC" in source and alias and not alias.startswith("UPI"):
                # Keep Swiss-Prot over TrEMBL if we already have one
                if ensp not in uniprot_map or "Swiss-Prot" in source:
                    uniprot_map[ensp] = alias

    print(f"  Gene symbols mapped: {len(gene_map)}")
    print(f"  UniProt IDs mapped: {len(uniprot_map)}")
    return gene_map, uniprot_map


def load_string(gene_map: dict, uniprot_map: dict) -> dict:
    """
    Load STRING physical interactions.

    Returns:
        interactions: {(gene_a, gene_b): {"string_score": int, "string_experimental": int,
                        "string_database": int, "uniprot_a": str, "uniprot_b": str}}
    """
    print("\n=== Loading STRING physical interactions ===")
    path = download(STRING_PHYSICAL, CACHE_DIR / "9606.protein.physical.links.detailed.txt.gz")

    interactions: dict[tuple, dict] = {}
    skipped_no_gene = 0
    skipped_low_score = 0
    total = 0

    with gzip.open(path, "rt") as f:
        reader = csv.DictReader(f, delimiter=" ")
        for row in reader:
            total += 1
            combined = int(row["combined_score"])
            if combined < STRING_MIN_SCORE:
                skipped_low_score += 1
                continue

            ensp1 = row["protein1"].replace("9606.", "")
            ensp2 = row["protein2"].replace("9606.", "")
            gene1 = gene_map.get(ensp1)
            gene2 = gene_map.get(ensp2)
            if not gene1 or not gene2:
                skipped_no_gene += 1
                continue

            # Canonical ordering (alphabetical) to deduplicate
            if gene1 > gene2:
                gene1, gene2 = gene2, gene1
                ensp1, ensp2 = ensp2, ensp1

            key = (gene1, gene2)
            if key in interactions:
                # Keep the higher score
                if combined > interactions[key]["string_score"]:
                    interactions[key]["string_score"] = combined
                continue

            interactions[key] = {
                "string_score": combined,
                "string_experimental": int(row.get("experimental", 0)),
                "string_database": int(row.get("database", 0)),
                "uniprot_a": uniprot_map.get(ensp1, ""),
                "uniprot_b": uniprot_map.get(ensp2, ""),
            }

    print(f"  Total rows: {total}")
    print(f"  Skipped (low score): {skipped_low_score}")
    print(f"  Skipped (no gene map): {skipped_no_gene}")
    print(f"  Loaded interactions: {len(interactions)}")
    return interactions


def load_biogrid(gene_map: dict) -> dict:
    """
    Load BioGRID human physical interactions.

    Returns:
        interactions: {(gene_a, gene_b): {"biogrid_systems": set_of_experimental_systems}}
    """
    print("\n=== Loading BioGRID interactions ===")
    path = download(BIOGRID_URL, CACHE_DIR / "BIOGRID-ORGANISM-LATEST.tab3.zip")

    # Find the human file inside the zip
    interactions: dict[tuple, dict] = defaultdict(lambda: {"biogrid_systems": set()})
    total = 0
    kept = 0

    with zipfile.ZipFile(path) as zf:
        human_files = [n for n in zf.namelist() if "Homo_sapiens" in n and n.endswith(".tab3.txt")]
        if not human_files:
            print("  ERROR: No human file found in BioGRID zip")
            print(f"  Available files: {zf.namelist()[:10]}")
            return {}

        fname = human_files[0]
        print(f"  Using: {fname}")

        with zf.open(fname) as raw:
            reader = csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"), delimiter="\t")
            for row in reader:
                total += 1
                # Only physical interactions
                if row.get("Experimental System Type") != "physical":
                    continue

                gene1 = row.get("Official Symbol Interactor A", "").strip()
                gene2 = row.get("Official Symbol Interactor B", "").strip()
                if not gene1 or not gene2 or gene1 == "-" or gene2 == "-":
                    continue

                # Canonical ordering
                if gene1 > gene2:
                    gene1, gene2 = gene2, gene1

                system = row.get("Experimental System", "")
                interactions[(gene1, gene2)]["biogrid_systems"].add(system)
                kept += 1

    print(f"  Total rows: {total}")
    print(f"  Physical interaction rows: {kept}")
    print(f"  Unique pairs: {len(interactions)}")
    return dict(interactions)


def load_huri() -> set:
    """
    Load HuRI interactions.

    Returns:
        set of (gene_a, gene_b) tuples (using Ensembl gene IDs, need mapping)
    """
    print("\n=== Loading HuRI interactions ===")
    path = download(HURI_URL, CACHE_DIR / "HuRI.tsv")

    pairs = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                pairs.add((parts[0].strip(), parts[1].strip()))

    print(f"  Loaded {len(pairs)} pairs")

    # Check what ID format they use
    sample = list(pairs)[:3]
    print(f"  Sample pairs: {sample}")
    return pairs


def build_ensembl_gene_to_symbol() -> dict:
    """Build mapping from Ensembl gene IDs to gene symbols using STRING aliases."""
    print("\n=== Building Ensembl Gene ID -> Symbol mapping ===")
    aliases_path = CACHE_DIR / "9606.protein.aliases.txt.gz"

    ensg_to_symbol: dict[str, str] = {}
    with gzip.open(aliases_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            alias = row["alias"]
            source = row["source"]
            ensp = row["#string_protein_id"].replace("9606.", "")

            # Collect Ensembl gene ID mappings
            if alias.startswith("ENSG") and "Ensembl" in source:
                # Map ENSG to the preferred gene symbol via protein info
                pass  # We'll use a different approach

    # Actually, let's build ENSG -> gene symbol from the info + aliases files
    # First: ENSP -> gene symbol (from info file)
    info_path = CACHE_DIR / "9606.protein.info.txt.gz"
    ensp_to_gene: dict[str, str] = {}
    with gzip.open(info_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ensp = row["#string_protein_id"].replace("9606.", "")
            ensp_to_gene[ensp] = row["preferred_name"]

    # Then: ENSG -> ENSP (from aliases)
    ensg_to_gene: dict[str, str] = {}
    with gzip.open(aliases_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            alias = row["alias"]
            source = row["source"]
            if alias.startswith("ENSG") and "Ensembl" in source:
                ensp = row["#string_protein_id"].replace("9606.", "")
                gene = ensp_to_gene.get(ensp)
                if gene:
                    ensg_to_gene[alias] = gene

    print(f"  Mapped {len(ensg_to_gene)} Ensembl gene IDs to symbols")
    return ensg_to_gene


def resolve_huri_ids(huri_pairs: set, ensg_to_gene: dict) -> dict:
    """
    Resolve HuRI Ensembl gene IDs to gene symbols.

    Returns:
        {(gene_a, gene_b): True} for all resolved pairs
    """
    resolved = {}
    unresolved = 0

    for id1, id2 in huri_pairs:
        gene1 = ensg_to_gene.get(id1, id1)
        gene2 = ensg_to_gene.get(id2, id2)

        # If still looks like an Ensembl ID, skip
        if gene1.startswith("ENSG") or gene2.startswith("ENSG"):
            unresolved += 1
            continue

        if gene1 > gene2:
            gene1, gene2 = gene2, gene1
        resolved[(gene1, gene2)] = True

    print(f"  Resolved HuRI pairs: {len(resolved)}")
    print(f"  Unresolved (skipped): {unresolved}")
    return resolved


def build_gene_to_uniprot(uniprot_map: dict, gene_map: dict) -> dict:
    """Build gene_symbol -> UniProt mapping."""
    result: dict[str, str] = {}
    for ensp, gene in gene_map.items():
        up = uniprot_map.get(ensp, "")
        if up and (gene not in result or len(up) < len(result[gene])):
            result[gene] = up
    return result


def merge_and_write(
    string_data: dict,
    biogrid_data: dict,
    huri_data: dict,
    gene_to_uniprot: dict,
):
    """Merge all sources and write the output file."""
    print("\n=== Merging all sources ===")

    # Collect all unique pairs
    all_pairs = set(string_data.keys()) | set(biogrid_data.keys()) | set(huri_data.keys())
    print(f"  Total unique pairs: {len(all_pairs)}")

    # Count overlaps
    s_keys = set(string_data.keys())
    b_keys = set(biogrid_data.keys())
    h_keys = set(huri_data.keys())
    print(f"  STRING only: {len(s_keys - b_keys - h_keys)}")
    print(f"  BioGRID only: {len(b_keys - s_keys - h_keys)}")
    print(f"  HuRI only: {len(h_keys - s_keys - b_keys)}")
    print(f"  STRING & BioGRID: {len(s_keys & b_keys)}")
    print(f"  STRING & HuRI: {len(s_keys & h_keys)}")
    print(f"  BioGRID & HuRI: {len(b_keys & h_keys)}")
    print(f"  All three: {len(s_keys & b_keys & h_keys)}")

    # Build output rows
    rows = []
    for gene_a, gene_b in sorted(all_pairs):
        s = string_data.get((gene_a, gene_b), {})
        b = biogrid_data.get((gene_a, gene_b), {})
        h = (gene_a, gene_b) in huri_data

        sources = []
        if s:
            sources.append("STRING")
        if b:
            sources.append("BioGRID")
        if h:
            sources.append("HuRI")

        uniprot_a = s.get("uniprot_a", "") or gene_to_uniprot.get(gene_a, "")
        uniprot_b = s.get("uniprot_b", "") or gene_to_uniprot.get(gene_b, "")

        biogrid_systems = "; ".join(sorted(b.get("biogrid_systems", set()))) if b else ""

        rows.append({
            "gene_a": gene_a,
            "gene_b": gene_b,
            "uniprot_a": uniprot_a,
            "uniprot_b": uniprot_b,
            "string_score": s.get("string_score", ""),
            "string_experimental": s.get("string_experimental", ""),
            "string_database": s.get("string_database", ""),
            "biogrid_systems": biogrid_systems,
            "huri": 1 if h else 0,
            "sources": "|".join(sources),
            "num_sources": len(sources),
        })

    # Write compressed CSV
    out_path = DATA_DIR / "interactions.csv.gz"
    print(f"\n=== Writing {out_path} ===")
    with gzip.open(out_path, "wt", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    size_mb = out_path.stat().st_size / 1e6
    print(f"  Wrote {len(rows)} interactions ({size_mb:.1f} MB compressed)")

    # Also write the gene->uniprot mapping for the app
    mapping_path = DATA_DIR / "gene_uniprot_map.csv.gz"
    print(f"\n=== Writing {mapping_path} ===")
    with gzip.open(mapping_path, "wt", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_symbol", "uniprot_id"])
        for gene, uniprot in sorted(gene_to_uniprot.items()):
            writer.writerow([gene, uniprot])

    size_mb = mapping_path.stat().st_size / 1e6
    print(f"  Wrote {len(gene_to_uniprot)} mappings ({size_mb:.1f} MB compressed)")

    # Write feather files to public/ for WASM deployment
    import pandas as pd

    public_dir = Path(__file__).parent / "public"
    public_dir.mkdir(exist_ok=True)

    pub_interactions = public_dir / "interactions.feather"
    print(f"\n=== Writing {pub_interactions} ===")
    pd.DataFrame(rows).to_feather(pub_interactions)
    size_mb = pub_interactions.stat().st_size / 1e6
    print(f"  Wrote {len(rows)} interactions ({size_mb:.1f} MB)")

    pub_mapping = public_dir / "gene_uniprot_map.feather"
    print(f"\n=== Writing {pub_mapping} ===")
    mapping_rows = [
        {"gene_symbol": gene, "uniprot_id": uniprot}
        for gene, uniprot in sorted(gene_to_uniprot.items())
    ]
    pd.DataFrame(mapping_rows).to_feather(pub_mapping)
    size_mb = pub_mapping.stat().st_size / 1e6
    print(f"  Wrote {len(gene_to_uniprot)} mappings ({size_mb:.1f} MB)")


def main():
    print("=" * 60)
    print("Building Protein-Protein Interaction Database")
    print("=" * 60)

    # Step 1: Build ID mappings
    gene_map, uniprot_map = build_id_mapping()

    # Step 2: Load STRING
    string_data = load_string(gene_map, uniprot_map)

    # Step 3: Load BioGRID
    biogrid_data = load_biogrid(gene_map)

    # Step 4: Load HuRI
    huri_pairs = load_huri()
    ensg_to_gene = build_ensembl_gene_to_symbol()
    huri_data = resolve_huri_ids(huri_pairs, ensg_to_gene)

    # Step 5: Build gene->uniprot lookup
    gene_to_uniprot = build_gene_to_uniprot(uniprot_map, gene_map)

    # Step 6: Merge and write
    merge_and_write(string_data, biogrid_data, huri_data, gene_to_uniprot)

    print("\n=== Done! ===")


if __name__ == "__main__":
    main()
