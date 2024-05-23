# Genetic-Survival-Analysis-in-UKB

This repository contains an scripts for performing detailed survival analysis using Cox proportional hazards and Kaplan-Meier plots. The script includes functions for loading and preparing phenotypic data, matching cases and controls, and generating survival plots. Key functionalities include:

- **Data Loading and Preparation:** Efficiently reads and processes phenotype and genotype data from UK Biobank
- **Survival Plotting:** Generates survival plots for time from diagnosis to cardiovascular outcomes using both Cox proportional hazards and Kaplan-Meier methods.
- **Matching Function:** Implements nearest neighbor matching to balance treatment and control groups based on various covariates.
- **Power Calculation:** Includes functions to perform power calculations for survival analysis.

## Table of Contents

- [Dependencies](##Dependencies)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Dependencies
Dependencies include `tidyverse`, `data.table`, `survminer`, `survival`, `MatchIt`, `expss`, and `cobalt`.

### Prerequisites

List any prerequisites needed to run the project.

```sh
# Example
pip install some_package
