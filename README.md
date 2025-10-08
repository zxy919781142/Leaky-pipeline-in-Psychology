Leaky Pipeline in Psychology






Overview

This repository contains data processing, analysis, and figure-generation code for the project “Failing in retention, not recruitment: Gender differences in academic progression in psychology.”

The project investigates gender differences in academic progression within the field of psychology, focusing on when gaps emerge and where retention shortfalls occur. Using large-scale bibliometric data, it examines academic age trajectories, publication performance, and cohort-based gender disparities.

Repository Structure

Leaky-pipeline-in-Psychology/
├─ R_code/ — R scripts for cleaning, modeling, and plotting
├─ data/ — Input and processed data (or placeholders)
├─ plots/ — Generated figures and outputs
├─ Leaky-pipeline-in-Psychology.Rproj
└─ README.md

Requirements

R version ≥ 4.2

Recommended packages:
tidyverse, data.table, janitor, lubridate, readr, stringr,
ggplot2, patchwork, broom, scales, and any modeling packages referenced in R_code/.

For reproducibility using renv:

renv::init()
renv::restore()

Workflow

Data Preparation
Scripts in R_code/ clean and organize raw bibliometric data stored in data/.

Modeling and Analysis
Time-varying effect models and cohort-based analyses are conducted to examine gendered publication and career patterns.

Visualization
Figures are created using ggplot2 and saved to plots/.

Each script includes detailed comments describing inputs, outputs, and parameters.

Data Notes

Raw data are not shared in this repository due to access restrictions.

Example or processed data can be added to the data/ folder for demonstration purposes.

All identifiers are anonymized to ensure confidentiality.

Example anonymization code:

data$id_anon <- digest::hmac("secret_salt", data$id)

Reproducibility

Random seeds are set with set.seed(12345) for consistent results.

Use renv or pak for package version control.

Optionally, automate the workflow using {targets} or {drake}.

Outputs

Cohort-level descriptive statistics

Gender-stratified career trajectories across academic age

Visualizations highlighting recruitment vs. retention patterns

All generated figures are stored in the plots/ folder.

Citation

If you use this repository, please cite:

Zhao, X. (2025). Failing in retention, not recruitment: Gender differences in academic progression in psychology. Max Planck Institute for Human Development.

(Replace with DOI or preprint link once available.)

Contributing

Contributions are welcome.
To contribute:

Open an issue describing your proposed change.

Reference relevant scripts or data.

Update documentation as needed.

License

This project is released under the MIT License.
See the LICENSE file for details.

Contact

Maintainer:
Xinyi Zhao
Max Planck Institute for Human Development

📧 Email: xinyi.zhao@mpib-berlin.mpg.de

🌐 GitHub: zxy919781142
