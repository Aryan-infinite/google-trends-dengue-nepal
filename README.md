
# Google Trends–Based Dengue Surveillance in Nepal

A reproducible Python-based research pipeline investigating the relationship between Google Trends search interest and laboratory-confirmed dengue cases in Nepal (2019–2024).

---

## Overview

This study explores whether Google Trends search activity can serve as an early indicator of dengue outbreaks in Nepal.

Using monthly Google Trends data and official EDCD dengue case reports, we performed correlation analyses to evaluate the association between online search behavior and confirmed dengue incidence.

---

## Research Questions

- Can Google Trends predict dengue outbreaks?
- Which search terms best correlate with dengue incidence?
- Can publicly available search data support digital disease surveillance?

---

## Dataset

**Period**

2019–2024

**Source**

- Google Trends
- Epidemiology and Disease Control Division (EDCD), Nepal

Variables:

- Dengue fever
- Fever
- Dengue Symptoms
- Confirmed Dengue Cases

---

## Statistical Analysis

Performed using Python.

Methods:

- Pearson Correlation
- Spearman Rank Correlation
- Scatter Plots with Linear Regression

---

## Repository Structure

```
google-trends-dengue-nepal/

├── data/
├── docs/
├── results/
│   ├── figures/
│   └── tables/
├── src/
├── README.md
├── requirements.txt
└── LICENSE
```

---

## Results

Outputs generated automatically:

- Correlation analysis table
- Publication-quality scatter plots

Figures are available inside:

```
results/figures/
```

Tables are available inside:

```
results/tables/
```

---

## Installation

```bash
git clone <repository-url>

cd google-trends-dengue-nepal

pip install -r requirements.txt
```

---

## Run Analysis

```bash
python src/analysis.py
```

---

## Authors

Aryan Bhusal

Independent Researcher

Bishal Neupane

Independent Researcher

Rita Neupane

Independent Researcher

---

## License

MIT License
=======
# google-trends-dengue-nepal
Correlation between Google Trends and Dengue Surveillance Data in Nepal (2019–2024)
>>>>>>> feee994d19baa61cba1f27aa1cf9a22f73963746
