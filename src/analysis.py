from tables import save_correlation_table
from plotting import scatter_plot
from stats_analysis import correlation
import pandas as pd

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