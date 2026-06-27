from scipy.stats import pearsonr, spearmanr


def correlation(x, y):
    """
    Calculates Pearson and Spearman correlation
    between two variables.
    """

    pearson_r, pearson_p = pearsonr(x, y)
    spearman_r, spearman_p = spearmanr(x, y)

    return {
        "Pearson r": float(round(pearson_r, 4)),
        "Pearson p": float(pearson_p),
        "Spearman rho": float(round(spearman_r, 4)),
        "Spearman p": float(spearman_p)
    }