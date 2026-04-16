# Protein Interaction Search

A browser-based tool to find predicted physical interactions between two lists of proteins. Built with [marimo](https://marimo.io/) and deployed as a static WebAssembly app.

## How it works

1. **Upload** two lists of human proteins (gene symbols or UniProt IDs)
2. **View** which pairs from the two lists have evidence for physical interaction
3. **Download** results as CSV with evidence details

## Data sources

The interaction database is built from three complementary sources:

| Source | Type | License |
|--------|------|---------|
| [STRING v12.0](https://string-db.org/) physical links | Integrated (experimental + predicted), score >= 400 | CC BY 4.0 |
| [BioGRID](https://thebiogrid.org/) | Literature-curated experimental | MIT |
| [HuRI](https://www.interactome-atlas.org/) | Systematic yeast two-hybrid | CC BY 4.0 |

The database is rebuilt from source data on each deployment via GitHub Actions.

## Development

### Local development

```bash
pip install marimo pandas
marimo edit app.py
```

### Rebuild the interaction database

```bash
pip install requests
python build_db.py
```

This downloads data from STRING, BioGRID, and HuRI, then merges into `public/interactions.csv`.

### WASM export (for testing)

```bash
marimo export html-wasm app.py -o dist/ --mode run --no-show-code
cd dist && python -m http.server
```

## Deployment

The app is automatically deployed to GitHub Pages on push to `main` via the workflow in `.github/workflows/deploy.yml`.

## Citations

- STRING: Szklarczyk et al., *Nucleic Acids Res.* 2023
- BioGRID: Oughtred et al., *Nucleic Acids Res.* 2021
- HuRI: Luck et al., *Nature* 2020
