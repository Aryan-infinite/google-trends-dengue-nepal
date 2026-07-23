from scipy.stats import pearsonr, spearmanr, shapiro



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


def shapiro_test(series):

    statistic, p_value = shapiro(series.dropna())

    if p_value < 0.001:
        p_display = "<0.001"
    else:
        p_display = round(p_value, 4)

    return {
        "W Statistic": round(statistic, 4),
        "p-value": p_display,
        "Normal": "Yes" if p_value > 0.05 else "No"
    }