import matplotlib.pyplot as plt
import seaborn as sns


def correlation_heatmap(data):

    columns = [
        "Dengue fever",
        "Fever",
        "Dengue Symptoms",
        "Confirmed Cases"
    ]

    corr = data[columns].corr(method="pearson")

    plt.figure(figsize=(7,6))

    sns.heatmap(
        corr,
        annot=True,
        cmap="RdYlBu_r",
        vmin=-1,
        vmax=1,
        linewidths=0.5,
        square=True,
        fmt=".2f"
    )

    plt.title("Pearson Correlation Heatmap")

    plt.tight_layout()

    plt.savefig(
        "results/figures/correlation_heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("✓ Correlation heatmap saved.")