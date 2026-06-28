import pandas as pd
from scipy.stats import pearsonr


def lag_correlation(data, search_term, max_lag=3):

    results = []

    for lag in range(-max_lag, max_lag + 1):

        shifted = data[search_term].shift(lag)

        temp = pd.DataFrame({
            "Search": shifted,
            "Cases": data["Confirmed Cases"]
        }).dropna()

        r, p = pearsonr(
            temp["Search"],
            temp["Cases"]
        )

        results.append({
            "Search Term": search_term,
            "Lag": lag,
            "Pearson r": round(r, 4),
            "p-value": round(p, 6)
        })

    return pd.DataFrame(results)