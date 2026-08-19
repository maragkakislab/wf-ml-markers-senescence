# Machine learning guided identification of senescence markers

This repository provides a reproducible workflow for identifying robust senescence markers from SenCat transcriptomic and proteomic data and for validating marker-based scoring in external IMR90 fibroblast datasets.

## Workflow

The workflow standardizes transcriptomic and proteomic measurements, applies consistent preprocessing, and uses a cross-cell-type machine learning strategy to identify markers that remain informative across biological contexts. A refined marker set is then used to derive stable marker weights and generate sample-level senescence scores in both reference and external validation datasets, with all primary outputs written to the analysis and plotting directories.

## Inputs

- **SenCat transcriptomic data**: primary RNA-level input used for marker discovery.
- **SenCat proteomic data**: primary protein-level input used for marker discovery.
- **External validation data**: IMR90 fibroblast datasets used to evaluate score transferability.
- **Workflow configuration**: centralized settings for inputs, analysis profiles, and output locations.

## ML markers

- **Transcriptomics markers**: `analysis/transcriptomics_loose5000f_tuned_common_features.csv`.
- **Proteomics markers**: `analysis/proteomics_loose5000f_tuned_common_features.csv`.

## Using ML markers for senescence scoring

You can use our senescence markers to score your data for senescence.

### Prepare data

The scoring is performed on `h5ad` file containing normalized transcriptomics or proteomics counts. Expected `h5ad` structure:

- `adata.X`: sample-by-feature expression matrix
- `adata.var_names`: feature identifiers matching the marker IDs in the marker CSV index

Both scripts below write a log file, so create the log directory first:

```bash
mkdir -p logs
```

If your data are not normalized, you can use `normalize_counts.py` script:

```python
python workflow/scripts/data/normalize_counts.py \
    --input-h5ad INPUT_H5AD \
    --design DESIGN_FACTORS \
    --output-h5ad NORMALIZED_H5AD \
    --log logs/my_data.normalize.log \
    --log-level INFO
```

- `INPUT_H5AD` specifies a path to your input `h5ad` file
- `DESIGN_FACTORS` specifies design factors for DESeq2, in the format `x + z` or `~x+z`.
- `NORMALIZED_H5AD` specifies a path where your normalized data will be saved

## Get senescence scores

```python
python workflow/scripts/cls/marker_classifier.py \
    --markers PATH_TO_ML_MARKERS \
    --input-h5ad NORMALIZED_H5AD \
    --output-results-csv OUTPUT_CSV \
    --log logs/my_data.score.log \
    --log-level INFO
```

- `PATH_TO_ML_MARKERS` specifies path to ML markers. Use `analysis/transcriptomics_loose5000f_tuned_common_features.csv` for transcriptomics and `analysis/proteomics_loose5000f_tuned_common_features.csv` for proteomics
- `NORMALIZED_H5AD` specifies path to your normalized `h5ad` data
- `OUTPUT_CSV` specifies path to output csv file with per-sample `score` values (higher values indicate stronger similarity to the senescence-associated signature).

Notes:

- `marker_classifier.py` applies `log1p` internally.
- Marker matching is based on `adata.var_names`; non-overlapping markers are skipped automatically.

## Manuscript

Anerillas, Carlos, et al. "SenCat: Cataloging human cell senescence through multiomic profiling of multiple senescent primary cell types." bioRxiv (2026): 2026-02. [https://doi.org/10.64898/2026.02.05.703986](https://doi.org/10.64898/2026.02.05.703986)