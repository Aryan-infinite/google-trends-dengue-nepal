import pandas as pd

def save_correlation_table(results_list):

    df = pd.DataFrame(results_list)

    # Keep only required columns
    df = df[
        [
            "Search Term",
            "Pearson r",
            "Pearson p",
            "Spearman rho",
            "Spearman p"
        ]
    ]

    # Round correlation coefficients
    df["Pearson r"] = df["Pearson r"].round(3)
    df["Spearman rho"] = df["Spearman rho"].round(3)

    # Format p-values
    df["Pearson p"] = df["Pearson p"].apply(
        lambda x: "<0.001" if x < 0.001 else round(x, 3)
    )

    df["Spearman p"] = df["Spearman p"].apply(
        lambda x: "<0.001" if x < 0.001 else round(x, 3)
    )

    # Interpretation column
    def interpret(r):
        r = abs(r)
        if r >= 0.80:
            return "Very strong positive"
        elif r >= 0.60:
            return "Strong positive"
        elif r >= 0.40:
            return "Moderate positive"
        elif r >= 0.20:
            return "Weak positive"
        else:
            return "Very weak"

    df["Interpretation"] = df["Pearson r"].apply(interpret)

    # Save final publication-ready table
    df.to_csv(
        "results/tables/Table_1_Correlation_Results.csv",
        index=False
    )

    print("✓ Publication-ready Table 1 saved.")