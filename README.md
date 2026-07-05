# Dengue Fever in Nepal: A Correlation Between Google Trends and National Surveillance Reports

This repository contains the complete analytical workflow used in the study **"Dengue Fever in Nepal: A Correlation Between Google Trends and National Surveillance Reports."** The project investigates the relationship between Google Trends search activity and laboratory-confirmed dengue cases reported in Nepal between **January 2019 and December 2024**. Using official surveillance data from the Epidemiology and Disease Control Division (EDCD) and publicly available Google Trends data, the study evaluates whether internet search behavior can complement traditional disease surveillance and support earlier detection of dengue outbreaks.

---

## Research Objectives

This study aimed to:

- Quantify the relationship between Google Trends Relative Search Volume (RSV) and laboratory-confirmed dengue cases in Nepal.
- Compare the performance of three dengue-related Google search terms ("Dengue fever", "Dengue Symptoms", and "Fever").
- Evaluate whether Google search activity precedes officially reported dengue cases using time-lag correlation analysis.
- Assess the potential of Google Trends as a complementary digital epidemiology tool for strengthening dengue surveillance in Nepal.

---

## Dataset

### Study Period

January 2019 – December 2024 (72 monthly observations)

### Data Sources

- **Google Trends** (https://trends.google.com/)
- **Epidemiology and Disease Control Division (EDCD), Department of Health Services, Ministry of Health and Population, Government of Nepal**

### Variables

- Monthly Relative Search Volume (RSV) for:
  - Dengue fever
  - Dengue Symptoms
  - Fever
- Monthly laboratory-confirmed dengue cases reported through the national surveillance system.

---

## Statistical Analysis

The analytical workflow was designed to ensure full reproducibility.

Python (Version 3.13) was used for data processing, descriptive statistics, correlation analysis, and visualization using the following libraries:

- Pandas
- NumPy
- SciPy
- Matplotlib

The analyses included:

- Data preprocessing and cleaning
- Descriptive statistics
- Pearson correlation analysis
- Spearman rank correlation analysis
- Two-tailed significance testing (p-values)
- Scatter plots with fitted linear regression lines
- Pearson correlation heatmap

Time-lag correlation analysis was performed separately in **Microsoft Excel** by shifting Google Trends Relative Search Volume (RSV) against laboratory-confirmed dengue cases across predefined monthly lead and lag intervals.

---

## Repository Structure

```text
google-trends-dengue-nepal/

├── data/
│   ├── raw/
│   └── processed/
│
├── docs/
│
├── results/
│   ├── figures/
│   └── tables/
│
├── src/
│
├── README.md
├── requirements.txt
└── LICENSE
```

---

## Reproducing the Analysis

Clone the repository:

```bash
git clone <repository-url>

cd google-trends-dengue-nepal
```

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Run the complete analysis:

```bash
python src/analysis.py
```

Running the script will reproduce the statistical analyses, tables, and publication-quality figures included in the manuscript.

---

## Outputs

The repository generates:

- Descriptive statistics
- Pearson and Spearman correlation tables
- Scatter plots with fitted regression lines
- Pearson correlation heatmap
- Publication-ready figures
- Statistical summary tables

Generated outputs are saved in:

```text
results/
├── figures/
└── tables/
```

---

## Authors

**Aryan Bhusal**  
Independent Researcher

**Bishal Neupane**  
Independent Researcher

**Rita Neupane**  
Independent Researcher

---

## Citation

If you use this repository, please cite:

> Bhusal A, Neupane B, Neupane R. *Dengue Fever in Nepal: A Correlation Between Google Trends and National Surveillance Reports.* GitHub repository.

---

## License

This project is distributed under the **MIT License**.
