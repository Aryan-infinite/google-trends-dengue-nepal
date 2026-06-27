import pandas as pd

def save_correlation_table(results_list):

    df = pd.DataFrame(results_list)

    df = df[
        [
            "Search Term",
            "Pearson r",
            "Pearson p",
            "Spearman rho",
            "Spearman p"
        ]
    ]

    df = df.round(4)

    df.to_csv(
        "results/tables/Table_1_Correlation.csv",
        index=False
    )

    print("✓ Table 1 saved.")