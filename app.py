import marimo

__generated_with = "0.13.0"
app = marimo.App(width="medium", app_title="Protein Interaction Search")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
        # Protein Interaction Search

        Upload two lists of proteins to find predicted physical interactions between them.

        **Data sources:** [STRING v12.0](https://string-db.org/) (physical links) |
        [BioGRID](https://thebiogrid.org/) | [HuRI](https://www.interactome-atlas.org/)

        Interactions are shown when at least one database reports evidence for a physical
        interaction between a protein from List A and a protein from List B.
        """
    )
    return


@app.cell
def _(mo):
    import pandas as pd

    with mo.status.spinner("Loading interaction database..."):
        _base = mo.notebook_location()
        interactions_df = pd.read_csv(str(_base / "public" / "interactions.csv"))
        gene_map_df = pd.read_csv(str(_base / "public" / "gene_uniprot_map.csv"))
    return gene_map_df, interactions_df, pd


@app.cell
def _(gene_map_df, interactions_df, pd):
    def _build_indices(idf, mdf):
        """Build lookup dictionaries from the loaded DataFrames."""
        interactions = {}
        gene_to_uniprot = {}
        alias_to_gene = {}

        for _, row in idf.iterrows():
            ga, gb = row["gene_a"], row["gene_b"]
            interactions[(ga, gb)] = row.to_dict()

            for gene, up_key in [(ga, "uniprot_a"), (gb, "uniprot_b")]:
                up = row.get(up_key, "")
                alias_to_gene[gene.upper()] = gene
                alias_to_gene[gene.lower()] = gene
                if pd.notna(up) and up:
                    up = str(up)
                    gene_to_uniprot[gene] = up
                    alias_to_gene[up.upper()] = gene
                    alias_to_gene[up.lower()] = gene

        for _, row in mdf.iterrows():
            gene = str(row["gene_symbol"])
            up = str(row["uniprot_id"])
            gene_to_uniprot[gene] = up
            alias_to_gene[gene.upper()] = gene
            alias_to_gene[gene.lower()] = gene
            alias_to_gene[up.upper()] = gene
            alias_to_gene[up.lower()] = gene

        return interactions, gene_to_uniprot, alias_to_gene

    interactions, gene_to_uniprot, alias_to_gene = _build_indices(
        interactions_df, gene_map_df
    )
    db_size = len(interactions)
    return alias_to_gene, db_size, gene_to_uniprot, interactions


@app.cell
def _(db_size, mo):
    mo.md(
        f"**Database loaded:** {db_size:,} protein pairs from STRING, BioGRID, and HuRI."
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## Input

        Enter protein identifiers below — one per line. You can use **gene symbols**
        (e.g. `TP53`, `BRCA1`) or **UniProt accessions** (e.g. `P04637`).
        Alternatively, upload a text file (one identifier per line).
        """
    )
    return


@app.cell
def _(mo):
    EXAMPLE_A = "TP53\nBRCA1\nEGFR\nKRAS\nRB1\nBCL2"
    EXAMPLE_B = "MDM2\nMDM4\nATM\nCHEK2\nGRB2\nBRAF\nBAX\nE2F1"

    load_example = mo.ui.run_button(label="Load example (cancer signaling)")
    return EXAMPLE_A, EXAMPLE_B, load_example


@app.cell
def _(mo):
    file_a = mo.ui.file(filetypes=[".txt", ".csv", ".tsv"], label="Upload List A", kind="area")
    file_b = mo.ui.file(filetypes=[".txt", ".csv", ".tsv"], label="Upload List B", kind="area")
    return file_a, file_b


@app.cell
def _(EXAMPLE_A, EXAMPLE_B, load_example, mo):
    text_a = mo.ui.text_area(
        label="List A (one protein per line)",
        placeholder="TP53\nBRCA1\nEGFR\nKRAS",
        value=EXAMPLE_A if load_example.value else "",
        full_width=True,
        rows=8,
    )
    text_b = mo.ui.text_area(
        label="List B (one protein per line)",
        placeholder="MDM2\nATM\nGRB2\nBRAF",
        value=EXAMPLE_B if load_example.value else "",
        full_width=True,
        rows=8,
    )
    return text_a, text_b


@app.cell
def _(file_a, file_b, load_example, mo, text_a, text_b):
    mo.vstack([
        mo.hstack([load_example], justify="center"),
        mo.hstack(
            [
                mo.vstack([mo.md("### List A"), file_a, text_a]),
                mo.vstack([mo.md("### List B"), file_b, text_b]),
            ],
            widths="equal",
            gap=2,
        ),
    ])
    return


@app.cell
def _():
    def parse_proteins(text: str) -> list[str]:
        """Parse a text block into a list of protein identifiers."""
        proteins = []
        for line in text.strip().splitlines():
            line = line.strip()
            # Handle CSV/TSV: take first column
            for sep in ["\t", ",", " "]:
                if sep in line:
                    line = line.split(sep)[0].strip()
                    break
            if line and not line.startswith("#"):
                proteins.append(line)
        return proteins
    return (parse_proteins,)


@app.cell
def _(file_a, file_b, parse_proteins, text_a, text_b):
    def _get_list(file_widget, text_widget):
        if file_widget.value:
            content = file_widget.value[0].contents.decode("utf-8", errors="ignore")
            return parse_proteins(content)
        return parse_proteins(text_widget.value or "")

    list_a_raw = _get_list(file_a, text_a)
    list_b_raw = _get_list(file_b, text_b)
    return list_a_raw, list_b_raw


@app.cell
def _(alias_to_gene, list_a_raw, list_b_raw):
    def _resolve(raw_list):
        resolved = {}
        unresolved = []
        for identifier in raw_list:
            gene = alias_to_gene.get(identifier.upper()) or alias_to_gene.get(identifier)
            if gene:
                resolved[gene] = identifier
            else:
                unresolved.append(identifier)
        return resolved, unresolved

    resolved_a, unresolved_a = _resolve(list_a_raw)
    resolved_b, unresolved_b = _resolve(list_b_raw)
    return resolved_a, resolved_b, unresolved_a, unresolved_b


@app.cell
def _(mo, resolved_a, resolved_b, unresolved_a, unresolved_b):
    _parts = []

    if resolved_a or resolved_b:
        _parts.append(
            f"**Resolved:** {len(resolved_a)} proteins in List A, "
            f"{len(resolved_b)} proteins in List B"
        )

    if unresolved_a:
        _parts.append(
            mo.callout(
                f"**List A — unrecognized:** {', '.join(unresolved_a[:20])}"
                + (" ..." if len(unresolved_a) > 20 else ""),
                kind="warn",
            )
        )
    if unresolved_b:
        _parts.append(
            mo.callout(
                f"**List B — unrecognized:** {', '.join(unresolved_b[:20])}"
                + (" ..." if len(unresolved_b) > 20 else ""),
                kind="warn",
            )
        )

    mo.vstack(_parts) if _parts else None
    return


@app.cell
def _(gene_to_uniprot, interactions, resolved_a, resolved_b):
    def _search():
        if not resolved_a or not resolved_b:
            return []

        results = []
        genes_a = set(resolved_a.keys())
        genes_b = set(resolved_b.keys())

        for ga in genes_a:
            for gb in genes_b:
                if ga == gb:
                    continue
                key = (ga, gb) if ga < gb else (gb, ga)
                row = interactions.get(key)
                if row:
                    results.append({
                        "Gene A": ga,
                        "Gene B": gb,
                        "UniProt A": gene_to_uniprot.get(ga, ""),
                        "UniProt B": gene_to_uniprot.get(gb, ""),
                        "STRING Score": row.get("string_score", ""),
                        "STRING Experimental": row.get("string_experimental", ""),
                        "STRING Database": row.get("string_database", ""),
                        "BioGRID Systems": row.get("biogrid_systems", ""),
                        "HuRI": "Yes" if str(row.get("huri", "0")) == "1" else "",
                        "Sources": str(row.get("sources", "")).replace("|", ", "),
                        "# Sources": row.get("num_sources", ""),
                    })

        results.sort(
            key=lambda r: (
                -int(r["# Sources"]) if r["# Sources"] else 0,
                -int(r["STRING Score"]) if r["STRING Score"] else 0,
            )
        )
        return results

    results = _search()
    return (results,)


@app.cell
def _(mo, resolved_a, resolved_b, results):
    _has_input = bool(resolved_a and resolved_b)

    if _has_input:
        _total_possible = len(resolved_a) * len(resolved_b)
        mo.md(
            f"""
            ## Results

            Found **{len(results)}** interactions out of {_total_possible:,} possible
            pairs ({len(resolved_a)} x {len(resolved_b)} proteins).
            """
        )
    else:
        mo.md(
            """
            ## Results

            Enter proteins in both lists above to search for interactions.
            """
        )
    return


@app.cell
def _(mo, pd, resolved_a, resolved_b, results):
    if results:
        results_df = pd.DataFrame(results)
        results_table = mo.ui.table(
            results_df,
            selection=None,
            page_size=20,
            label="Interacting protein pairs",
        )
    else:
        results_df = pd.DataFrame()
        if resolved_a and resolved_b:
            results_table = mo.callout(
                mo.md(
                    "**No interactions found** between these two lists.\n\n"
                    "This means none of the databases (STRING, BioGRID, HuRI) "
                    "report a physical interaction for any pair across your two lists. "
                    "Try broadening your lists or checking for typos."
                ),
                kind="info",
            )
        else:
            results_table = mo.md(
                "*Enter proteins in both lists above to search for interactions.*"
            )

    results_table
    return (results_df,)


@app.cell
def _(mo, results, results_df):
    if results:
        _csv_data = results_df.to_csv(index=False)
        mo.download(
            data=_csv_data.encode("utf-8"),
            filename="protein_interactions.csv",
            mimetype="text/csv",
            label="Download results as CSV",
        )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ---

        ### About

        This tool searches for known and predicted physical protein-protein interactions
        across three complementary databases:

        | Source | Description | Evidence type |
        |--------|-------------|---------------|
        | **STRING** | Integrated interaction database (physical sub-network, score >= 400) | Combined computational + experimental |
        | **BioGRID** | Literature-curated interactions | Experimental methods (AP-MS, Y2H, etc.) |
        | **HuRI** | Human Reference Interactome | Systematic yeast two-hybrid screens |

        **STRING score interpretation:** 0-1000 scale.
        Low confidence >= 150, Medium >= 400, High >= 700, Highest >= 900.

        **How to cite:**
        - STRING: Szklarczyk et al., *Nucleic Acids Res.* 2023
        - BioGRID: Oughtred et al., *Nucleic Acids Res.* 2021
        - HuRI: Luck et al., *Nature* 2020
        """
    )
    return


if __name__ == "__main__":
    app.run()
