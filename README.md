# TRIUMPH-HD-Flow Repository

This repository contains the analysis code and metadata for high-dimensional immune profiling of mucosal flow cytometry data.

## Repository structure

- `scripts/`
  - `Mucosal_Master_Pipeline.R` — main pipeline for data import, quality control, clustering, UMAP, and figure generation.
  - `TRIUMPH_High_Dimensional_Flow_Data_Analysis.R` — alternate/high-dimensional analysis workflow.
- `metadata_corrected.csv` — sample metadata used by the analysis scripts.
- `.gitignore` — excludes raw FCS files, figures, results, and temporary files from version control.

## Notes

- Raw FCS data and generated outputs are intentionally excluded from the repository.
- The scripts are configured to run from the local working directory:
  `C:/Users/LENOVO/OneDrive - Malawi-Liverpool Wellcome Research Programme/Desktop/LAB WORK/Exported FCS Files for Mucosal Analysis/High Dimensionality reductionl-2026/DownSample`
- Running the scripts creates output folders such as `Results/` and `Figures/` locally, but these folders are ignored by `.gitignore`.

## How to use

1. Clone the repository locally or pull the latest changes.
2. Place raw FCS files in the local working directory configured by the scripts.
3. Open one of the scripts in RStudio or VS Code.
4. Run the pipeline sections in order (or adjust the `run_stages` / `RUN_STAGES` variables to rerun only the desired stages).

## GitHub

This repository is linked to `https://github.com/rknyirenda/TRIUMPH-HD-Flow`.
