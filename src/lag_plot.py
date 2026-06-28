import matplotlib.pyplot as plt


def lag_plot(df, search_term):

    temp = df[df["Search Term"] == search_term]

    plt.figure(figsize=(7,4))

    plt.plot(
        temp["Lag"],
        temp["Pearson r"],
        marker="o",
        linewidth=2
    )

    plt.axhline(0,color="black",linewidth=0.8)

    plt.xlabel("Lag (Months)")
    plt.ylabel("Pearson Correlation")

    plt.title(f"Time-lag Correlation\n{search_term}")

    plt.grid(alpha=0.3)

    filename = (
        "results/figures/"
        + search_term.replace(" ","_").lower()
        + "_lag.png"
    )

    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()