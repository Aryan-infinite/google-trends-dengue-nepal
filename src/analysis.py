from tables import save_correlation_table
from plotting import scatter_plot
from stats_analysis import correlation
import pandas as pd
from time_lag import lag_correlation
from lag_plot import lag_plot
from heatmap import correlation_heatmap

print("="*50)
print("Google Trends × Dengue Nepal")
print("="*50)

file_path = "data/raw/finaldata.xlsx"

data = pd.read_excel(file_path)

print("\nFirst 5 rows:")
print(data.head())

print("\nDataset Information")
print(data.info())

print("\nSummary Statistics")
print(data.describe())

print("\nMissing Values")
print(data.isnull().sum())

results = correlation(
    data["Dengue fever"],
    data["Confirmed Cases"]
)

print("\nCorrelation Results")
print(results)
print("\n")
print("=" * 50)
print("CORRELATION ANALYSIS")
print("=" * 50)

search_terms = [
    "Dengue fever",
    "Fever",
    "Dengue Symptoms"
]

results_list = []

for term in search_terms:

    print("\n-----------------------------------")
    print(term)
    print("-----------------------------------")

    results = correlation(
        data[term],
        data["Confirmed Cases"]
    )

    print(results)

    results["Search Term"] = term

    # THIS MUST BE INSIDE THE LOOP
    results_list.append(results)

# Save correlation table
save_correlation_table(results_list)

print("\nCreating publication-quality figures...")

for term in search_terms:
    scatter_plot(data, term)

print("Done.")
print("\nGenerating time-lag analysis...")

lag_results = []

for term in search_terms:
    df = lag_correlation(
        data,
        term,
        max_lag=3
    )
    lag_results.append(df)

lag_results = pd.concat(
    lag_results,
    ignore_index=True
)

lag_results.to_csv(
    "results/tables/Table_2_Time_Lag_Correlation.csv",
    index=False
)

for term in search_terms:
    lag_plot(
        lag_results,
        term
    )

print("✓ Time-lag analysis completed.")

print("\nCreating correlation heatmap...")

correlation_heatmap(data)

print("✓ Analysis complete.")