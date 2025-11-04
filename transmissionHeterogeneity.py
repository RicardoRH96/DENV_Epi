import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path("Scripts").resolve()))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import baltic as bt
from transmission_count import (
    tree_to_table,
    transmission_time_series,
    accumulated_indicator,
    get_window_size,
    count_spatial_transmission_linkages,
)

FOCAL_LOCATION = "Colombia"


def sss_from_table(table: pd.DataFrame, location: str) -> float:
    """SSS = (exports - imports) / (exports + imports),
    counting only true between-state events for the focal location."""
    if table.empty or not {"From","To"}.issubset(table.columns):
        return np.nan
    # Only transitions that change state
    between = table[table["From"] != table["To"]]
    # Focal exports/imports
    exp_mask = (between["From"] == location) & (between["To"] != location)
    imp_mask = (between["To"] == location) & (between["From"] != location)
    exports = exp_mask.sum()
    imports = imp_mask.sum()
    denom = exports + imports
    return (exports - imports) / denom if denom else np.nan


def _focal_edges(table: pd.DataFrame, location: str) -> pd.DataFrame:
    if table.empty or not {"From","To"}.issubset(table.columns):
        return table.iloc[0:0]
    between = table[table["From"] != table["To"]]
    exp_mask = (between["From"] == location) & (between["To"] != location)
    imp_mask = (between["To"] == location) & (between["From"] != location)
    return between[exp_mask | imp_mask].copy()

def bootstrap_sss_events(table: pd.DataFrame, location: str, *, iterations=5000, seed=123) -> np.ndarray:
    sub = _focal_edges(table, location)
    if sub.empty:
        return np.array([])
    rng = np.random.default_rng(seed)
    n = len(sub)
    boot = np.empty(iterations)
    for i in range(iterations):
        idx = rng.integers(0, n, size=n)
        boot[i] = sss_from_table(sub.iloc[idx], location)
    return boot

def bootstrap_sss_time_stratified(
    table: pd.DataFrame, location: str, *,
    time_col: str = "date", freq: str = "Q",
    iterations: int = 5000, seed: int = 456
) -> np.ndarray:
    sub = _focal_edges(table, location)
    if sub.empty or time_col not in sub.columns:
        return np.array([])
    sub = sub.copy()
    sub[time_col] = pd.to_datetime(sub[time_col], errors="coerce")
    sub["__time_bin__"] = sub[time_col].dt.to_period(freq)
    groups = [(b, g) for b, g in sub.groupby("__time_bin__", dropna=False) if not g.empty]
    if not groups:
        return np.array([])
    rng = np.random.default_rng(seed)
    boot = np.empty(iterations)
    for i in range(iterations):
        parts = []
        for _, g in groups:
            idx = rng.integers(0, len(g), size=len(g))
            parts.append(g.iloc[idx])
        boot_tbl = pd.concat(parts, ignore_index=True)
        boot[i] = sss_from_table(boot_tbl, location)
    return boot


def jackknife_sss_tips(
    tree_path: Path,
    trait: str,
    location: str,
    branch_cutoff: float,
    *,
    drop_fraction: float = 0.15,
    iterations: int = 400,
    seed: int = 789,
    absolute_time: bool = True,
    sort_branches: bool = True,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    results = np.empty(iterations)

    for i in range(iterations):
        tree = bt.loadNexus(
            str(tree_path),
            absoluteTime=absolute_time,
            sortBranches=sort_branches,
        )
        tips = tree.getExternal()
        if not tips:
            results[i] = np.nan
            continue

        keep_count = max(1, int(round(len(tips) * (1 - drop_fraction))))
        keep_indices = rng.choice(len(tips), size=keep_count, replace=False)
        keep_branches = [tips[j] for j in keep_indices]
        reduced = tree.reduceTree(keep_branches)
        nexus_payload = reduced.toString(nexus=True)

        with tempfile.NamedTemporaryFile("w", suffix=".nexus", delete=False) as tmp:
            tmp.write(nexus_payload)
            tmp_path = Path(tmp.name)

        try:
            resampled_table = tree_to_table(
                tmp_path,
                trait,
                branch_cutoff,
                absolute_time=absolute_time,
                sort_branches=sort_branches,
            )
            results[i] = sss_from_table(resampled_table, location)
        finally:
            tmp_path.unlink(missing_ok=True)

    return results


def summarize_bootstrap(name: str, point_estimate: float, samples: np.ndarray) -> None:
    if samples.size == 0 or np.all(np.isnan(samples)):
        print(f"{name}: insufficient data for bootstrap.")
        return
    clean = samples[~np.isnan(samples)]
    if clean.size == 0:
        print(f"{name}: all bootstrap iterations returned NaN.")
        return
    lo, hi = np.percentile(clean, [2.5, 97.5])
    print(f"{name}: point={point_estimate:.4f}, 95% bootstrap CI [{lo:.4f}, {hi:.4f}]")


def build_bootstrap_dataframe(
    datasets: dict[str, dict[str, np.ndarray]]
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset_label, methods in datasets.items():
        for method_label, samples in methods.items():
            if samples.size == 0:
                continue
            clean = samples[~np.isnan(samples)]
            if clean.size == 0:
                continue
            rows.extend(
                {
                    "dataset": dataset_label,
                    "method": method_label,
                    "SSS": value,
                }
                for value in clean
            )
    return pd.DataFrame(rows)




OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442"]

def plot_bootstrap_distributions(df: pd.DataFrame, output_path, *, width=13, height=5, dpi=300) -> None:
    if df.empty:
        print("No bootstrap samples available for visualization.")
        return

    sns.set_style("whitegrid")
    dataset_order = ["No subsampling", "Temporal + Geographical subsampling"]

    fig, axes = plt.subplots(1, len(dataset_order), figsize=(width, height), sharey=True)
    axes = np.atleast_1d(axes)

    for ax, dataset in zip(axes, dataset_order):
        sub = df[df["dataset"] == dataset]
        if sub.empty:
            ax.set_visible(False)
            continue

        # Horizontal violins: SSS on x, methods on y
        sns.violinplot(
            data=sub,
            x="SSS",
            y="method",
            hue="method",              # keep palette without deprecation
            legend=False,
            palette=OKABE_ITO[:sub["method"].nunique()],
            inner=None,                # we’ll draw our own summaries
            cut=0,
            density_norm="width",
            linewidth=0.8,
            ax=ax,
        )

        # Bold vertical line at zero (neutral SSS)
        ax.axvline(0, color="black", linestyle=(0, (4, 3)), linewidth=1.2, alpha=0.9)

        # Overlay mean and 95% CI per method
        for i, m in enumerate(sub["method"].unique()):
            vals = sub.loc[sub["method"] == m, "SSS"].to_numpy()
            if vals.size == 0:
                continue
            mean = np.mean(vals)
            lo, hi = np.quantile(vals, [0.025, 0.975])

            # y-position of the i-th category:
            y = i
            ax.plot([mean], [y], marker="o", markersize=5, color="k", zorder=5)
            ax.hlines(y, lo, hi, color="k", linewidth=1.5, zorder=5)
            ax.plot([lo, hi], [y, y], marker="|", color="k", linewidth=0, markersize=8, zorder=5)

        # Aesthetics
        ax.set_title(dataset, fontsize=14, pad=8)
        ax.set_ylabel("Bootstrap strategy", fontsize=12)
        ax.set_xlabel("Source–Sink Score (SSS)", fontsize=12)
        ax.tick_params(axis="y", labelsize=10)
        ax.tick_params(axis="x", labelsize=10)
        ax.set_xlim(min(-1.05, sub["SSS"].min()*1.1), max(1.05, sub["SSS"].max()*1.1))

        # Light grid; keep axes clean
        sns.despine(ax=ax)

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

tree_path_full = Path("Revision/dengue_denv2_genome_dates.nexus")
trait_full = "country"
#Define the cutoff using the whole epidemiological window (tree height * 365 days)
branch_cutoff_full = get_window_size(tree_path_full)
print(f"Using branch cutoff of {branch_cutoff_full} days")

migration_table_full = tree_to_table(tree_path_full, trait_full, branch_cutoff_full)
migration_table_full.to_csv("migration_table_full.csv", index=False)
print(migration_table_full)

series = transmission_time_series(migration_table_full, FOCAL_LOCATION)
print(series.groupby("Transmission")["num"].sum())          # totals per category
print(series[series["num"] > 0])                      # first non‑zero weeks

lis = accumulated_indicator(migration_table_full, [FOCAL_LOCATION], indicator="LIS")
sss = accumulated_indicator(migration_table_full, [FOCAL_LOCATION], indicator="SSS")
print("LIS:", lis[FOCAL_LOCATION])
print("SSS:", sss[FOCAL_LOCATION])

#add now the subsampled BEAST tree
tree_path_beast = Path("Revision/newMCCTree_F1112_Markov.mcc.tree")
trait_beast = "location.rate"
#Define the cutoff using the whole epidemiological window (tree height * 365 days)
branch_cutoff_beast = get_window_size(tree_path_beast)
print(f"Using branch cutoff of {branch_cutoff_beast} days")

migration_table_beast = tree_to_table(tree_path_beast, trait_beast, branch_cutoff_beast)
migration_table_beast.to_csv("migration_table.csv", index=False)
print(migration_table_beast)

series = transmission_time_series(migration_table_beast, FOCAL_LOCATION)
print(series.groupby("Transmission")["num"].sum())          # totals per category
print(series[series["num"] > 0])                      # first non‑zero weeks

lis = accumulated_indicator(migration_table_beast, [FOCAL_LOCATION], indicator="LIS")
sss = accumulated_indicator(migration_table_beast, [FOCAL_LOCATION], indicator="SSS")
print("LIS (BEAST):", lis[FOCAL_LOCATION])
print("SSS (BEAST):", sss[FOCAL_LOCATION])

print("\n=== Bootstrapping SSS (Full ML tree) ===")
sss_full_point = sss_from_table(migration_table_full, FOCAL_LOCATION)
event_boot_full = bootstrap_sss_events(
    migration_table_full,
    FOCAL_LOCATION,
    iterations=500,
    seed=123,
)
summarize_bootstrap("Event bootstrap (full tree)", sss_full_point, event_boot_full)

time_boot_full = bootstrap_sss_time_stratified(
    migration_table_full,
    FOCAL_LOCATION,
    freq="Q",
    iterations=500,
    seed=456,
)
summarize_bootstrap(
    "Time-stratified bootstrap (full tree)", sss_full_point, time_boot_full
)

jackknife_full = jackknife_sss_tips(
    tree_path_full,
    trait_full,
    FOCAL_LOCATION,
    branch_cutoff_full,
    drop_fraction=0.15,
    iterations=100,
    seed=789,
)
summarize_bootstrap("Tip jackknife (full tree)", sss_full_point, jackknife_full)

print("\n=== Bootstrapping SSS (BEAST MCC tree) ===")
sss_beast_point = sss_from_table(migration_table_beast, FOCAL_LOCATION)
event_boot_beast = bootstrap_sss_events(
    migration_table_beast,
    FOCAL_LOCATION,
    iterations=500,
    seed=321,
)
summarize_bootstrap("Event bootstrap (BEAST)", sss_beast_point, event_boot_beast)

time_boot_beast = bootstrap_sss_time_stratified(
    migration_table_beast,
    FOCAL_LOCATION,
    freq="Q",
    iterations=500,
    seed=654,
)
summarize_bootstrap(
    "Time-stratified bootstrap (BEAST)", sss_beast_point, time_boot_beast
)

jackknife_beast = jackknife_sss_tips(
    tree_path_beast,
    trait_beast,
    FOCAL_LOCATION,
    branch_cutoff_beast,
    drop_fraction=0.15,
    iterations=100,
    seed=987,
)
summarize_bootstrap("Tip jackknife (BEAST)", sss_beast_point, jackknife_beast)

bootstrap_df = build_bootstrap_dataframe(
    {
        "No subsampling": {
            "Event bootstrap": event_boot_full,
            "Time-stratified bootstrap": time_boot_full,
            "Tip jackknife": jackknife_full,
        },
        "Temporal + Geographical subsampling": {
            "Event bootstrap": event_boot_beast,
            "Time-stratified bootstrap": time_boot_beast,
            "Tip jackknife": jackknife_beast,
        },
    }
)

plot_bootstrap_distributions(bootstrap_df, Path("SSS_bootstrap_distributions.png"))
print("Saved SSS bootstrap distributions figure to SSS_bootstrap_distributions.png")
