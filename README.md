# qPCRAnalyzer

A Shiny App for easy calculation of qPCR data

## Installation

```R
# Install the development version from GitHub:
require("devtools")
devtools::install_github("j-y26/qPCRAnalyzer")
```

## Shiny App Usage

```R
require(qPCRAnalyzer)
qPCRAnalyzer::run_qPCRAnalyzer()
```

To run:
1. Load your data file (Excel, Tab-delimited, or CSV format)

2. Check for any possible invalid wells and select to exclude them from analysis

3. Run the analysis and you will get your results!