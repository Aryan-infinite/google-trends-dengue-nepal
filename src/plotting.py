import matplotlib.pyplot as plt
import numpy as np

def scatter_plot(data, term):

    plt.figure(figsize=(6,5))

    plt.scatter(
        data[term],
        data["Confirmed Cases"]
    )

    m, b = np.polyfit(
        data[term],
        data["Confirmed Cases"],
        1
    )

    plt.plot(
        data[term],
        m*data[term]+b
    )

    plt.xlabel(term)
    plt.ylabel("Confirmed Cases")
    plt.title(f"{term} vs Confirmed Cases")

    filename = (
        "results/figures/"
        + term.replace(" ", "_").lower()
        + "_scatter.png"
    )

    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()