"""
Utilities translated from `transmissionCount.R` for working with baltic trees.

The helpers expect a time-resolved tree in NEXUS format whose tips and
internal nodes carry the discrete trait used to track transitions
(`region`, `division`, etc.).  Branch lengths are assumed to be in years.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import baltic as bt

DEFAULT_COLOR = "#D3D3D3"
IMPORT_COLOR = "#22C55E"
LOCAL_COLOR = "#e09f3e"
EXPORT_COLOR = "#9c6ade"

NORTH_COUNTRIES: Set[str] = {
    "Mexico",
    "Guatemala",
    "Belize",
    "Honduras",
    "El Salvador",
    "Nicaragua",
    "Costa Rica",
    "Panama",
    "Jamaica",
    "Haiti",
    "Dominican Republic",
    "Cuba",
    "Puerto Rico",
    "Bahamas",
    "Trinidad and Tobago",
}

SOUTH_COUNTRIES: Set[str] = {
    "Colombia",
    "Venezuela",
    "Ecuador",
    "Peru",
    "Bolivia",
    "Brazil",
    "Paraguay",
    "Uruguay",
    "Argentina",
    "Chile",
    "Guyana",
    "Suriname",
    "French Guiana",
}


def _north_south_transition_category(
    origin: Optional[str],
    destination: Optional[str],
    target: str,
) -> Optional[str]:
    """
    Map an origin/destination pair into a north/south transition category.
    """
    if not origin or not destination:
        return None
    if destination == target:
        if origin in NORTH_COUNTRIES:
            return "north_import"
        if origin in SOUTH_COUNTRIES:
            return "south_import"
    if origin == target:
        if destination in NORTH_COUNTRIES:
            return "north_export"
        if destination in SOUTH_COUNTRIES:
            return "south_export"
    return None


def decimal_year_to_date(decimal_year: float) -> pd.Timestamp:
    """
    Convert a decimal year (e.g. 2024.25) into a pandas timestamp.
    """
    if pd.isna(decimal_year):
        return pd.NaT
    year = int(decimal_year)
    fraction = decimal_year - year
    start_of_year = pd.Timestamp(f"{year}-01-01")
    days_in_year = (pd.Timestamp(f"{year + 1}-01-01") - start_of_year).days
    return start_of_year + pd.to_timedelta(fraction * days_in_year, unit="D")

def get_window_size(tree_path: str | Path) -> int:
    """
    Calculate the epidemiological window size in days given the tree height in years.
    """
    tree = bt.loadNexus(str(tree_path))
    return int(tree.treeHeight * 365)


def _tree_to_records(tree: bt.tree, trait: str) -> List[Mapping[str, object]]:
    """
    Flatten a baltic tree into a list of records containing parent/child trait labels.
    """
    records = []
    for branch in tree.Objects:
        parent = branch.parent
        parent_name = getattr(parent, "name", None) if parent else None
        parent_trait = parent.traits.get(trait) if parent else None
        record = {
            "node": getattr(branch, "name", None) or getattr(branch, "index", None),
            "parent": parent_name,
            "num_date": getattr(branch, "absoluteTime", np.nan),
            trait: branch.traits.get(trait),
            "branch_length": branch.length if branch.length is not None else 0.0,
            "from_value": parent_trait,
        }
        records.append(record)
    return records


def tree_to_table(
    path: str | Path,
    trait: str,
    branch_cutoff: float | int = np.inf,
    *,
    absolute_time: bool = True,
    sort_branches: bool = True,
) -> pd.DataFrame:
    """
    Reproduce the R `treeToTable` helper using baltic.

    Parameters
    ----------
    path:
        Path to a NEXUS tree (`bt.loadNexus` compatible).
    trait:
        Discrete trait to extract (e.g. 'region', 'division').
    branch_cutoff:
        Maximum branch length (in days) to keep, analogous to the `bc`
        argument in the R implementation.
    absolute_time:
        Passed to `bt.loadNexus`.
    sort_branches:
        Passed to `bt.loadNexus`.

    Returns
    -------
    pandas.DataFrame
        Columns: `Epiweek`, `From`, `To`, `BranchLength`, `num_date`.
    """
    tree = bt.loadNexus(
        str(path),
        absoluteTime=absolute_time,
        sortBranches=sort_branches,
    )
    records = _tree_to_records(tree, trait)
    df = pd.DataFrame.from_records(records)
    df["num_date"] = pd.to_numeric(df["num_date"], errors="coerce")
    df["date"] = df["num_date"].apply(decimal_year_to_date)

    iso_weeks = df["date"].apply(
        lambda d: f"{d.isocalendar().year}{d.isocalendar().week:02d}"
        if pd.notna(d)
        else None
    )
    df["Epiweek"] = iso_weeks
    df.rename(columns={"from_value": "From", trait: "To"}, inplace=True)
    df["BranchLength"] = df["branch_length"] * 365.0
    df = df[["Epiweek", "From", "To", "BranchLength", "num_date", "date"]]
    df = df[df["BranchLength"] < branch_cutoff]
    df = df.dropna(subset=["Epiweek", "From", "To"])
    return df.reset_index(drop=True)


def count_spatial_transmission_linkages(
    migration_table: pd.DataFrame,
    category: str,
    location: Sequence[str] | str,
) -> pd.Series | pd.DataFrame:
    """
    Translate `countSaptialTransmissionLinkages` into Python.

    Parameters
    ----------
    migration_table:
        Output from ``tree_to_table``.
    category:
        One of ``"importation"``, ``"transmission linkage"``, ``"exportation"``,
        or ``"north2south"``.
    location:
        Regions to consider. For ``"north2south"`` this must be a single focal
        location.

    Returns
    -------
    pandas.Series or pandas.DataFrame
        Counts per epiweek for the requested category. When ``category`` is
        ``"north2south"`` a DataFrame is returned with columns describing whether
        the transition involved the focal location and a northern or southern
        partner (imports/exports).
    """
    epiweeks = sorted(migration_table["Epiweek"].unique())
    counts = pd.Series(0, index=epiweeks, dtype=int)

    if isinstance(location, str):
        locations = [location]
    else:
        locations = list(location)

    data = migration_table

    if category == "importation":
        if len(locations) != 1:
            raise ValueError("`location` must contain exactly one region for importation.")
        mask = (data["To"] == locations[0]) & (data["From"] != locations[0])
    elif category == "transmission linkage":
        if len(locations) != 2:
            raise ValueError(
                "`location` must contain two regions for a transmission linkage."
            )
        mask = (data["From"] == locations[0]) & (data["To"] == locations[1])
    elif category == "exportation":
        if len(locations) != 1:
            raise ValueError("`location` must contain exactly one region for exportation.")
        mask = (data["From"] == locations[0]) & (data["To"] != locations[0])
    elif category == "north2south":
        if len(locations) != 1:
            raise ValueError("`location` must contain exactly one region for north2south.")
        target = locations[0]
        epiweek_index = pd.Index(epiweeks, name="Epiweek")
        direction_labels = [
            "north_import",
            "south_import",
            "north_export",
            "south_export",
        ]
        result = pd.DataFrame(0, index=epiweek_index, columns=direction_labels, dtype=int)

        direction_masks = {
            "north_import": (data["To"] == target) & data["From"].isin(NORTH_COUNTRIES),
            "south_import": (data["To"] == target) & data["From"].isin(SOUTH_COUNTRIES),
            "north_export": (data["From"] == target) & data["To"].isin(NORTH_COUNTRIES),
            "south_export": (data["From"] == target) & data["To"].isin(SOUTH_COUNTRIES),
        }

        for label, label_mask in direction_masks.items():
            grouped = data.loc[label_mask].groupby("Epiweek").size()
            if not grouped.empty:
                result.loc[grouped.index, label] = grouped.astype(int)

        return result
    else:
        raise ValueError(
            "`category` must be importation, transmission linkage, exportation, or north2south."
        )

    filtered = data.loc[mask]
    if filtered.empty:
        return counts

    grouped = filtered.groupby("Epiweek").size()
    counts.update(grouped.astype(int))
    return counts


def accumulated_indicator(
    migration_table: pd.DataFrame,
    locations: Iterable[str],
    indicator: str,
) -> pd.Series:
    """
    Translate the R `accumulatedIndicator` helper.
    """
    values = {}
    for location in locations:
        imports = count_spatial_transmission_linkages(
            migration_table, "importation", location
        )
        if indicator == "LIS":
            local = count_spatial_transmission_linkages(
                migration_table, "transmission linkage", (location, location)
            )
            numerator = imports.sum()
            denominator = numerator + local.sum()
        elif indicator == "SSS":
            exports = count_spatial_transmission_linkages(
                migration_table, "exportation", location
            )
            numerator = exports.sum() - imports.sum()
            denominator = exports.sum() + imports.sum()
        else:
            raise ValueError("`indicator` must be either 'LIS' or 'SSS'.")

        values[location] = numerator / denominator if denominator else np.nan
    return pd.Series(values)


def epiweek_to_date(epiweek: str) -> pd.Timestamp:
    """
    Convert an `YYYYWW` epiweek string to the Monday of that ISO week.
    """
    if pd.isna(epiweek):
        return pd.NaT
    year = int(epiweek[:4])
    week = int(epiweek[4:])
    return pd.Timestamp.fromisocalendar(year, week, 1)


def transmission_time_series(
    migration_table: pd.DataFrame,
    location: str,
    time_unit: str = "week",
) -> pd.DataFrame:
    """
    Summarise transmission counts per category over time.

    Parameters
    ----------
    migration_table : pandas.DataFrame
        Output of ``tree_to_table``.
    location : str
        Focal location used to partition import/export/local counts.
    time_unit : str
        Either ``"week"`` (default) or ``"month"`` to aggregate to calendar months.
    """
    series = []
    for label in ("importation", "localtransmission", "exportation"):
        if label == "localtransmission":
            counts = count_spatial_transmission_linkages(
                migration_table, "transmission linkage", (location, location)
            )
        else:
            counts = count_spatial_transmission_linkages(
                migration_table,
                "importation" if label == "importation" else "exportation",
                location,
            )
        tmp = pd.DataFrame(
            {
                "Epiweek": counts.index,
                "num": counts.values,
                "Transmission": label,
            }
        )
        series.append(tmp)

    df = pd.concat(series, ignore_index=True)
    df["date"] = df["Epiweek"].apply(epiweek_to_date)

    if time_unit == "week":
        return df

    if time_unit == "month":
        df["month"] = df["date"].dt.to_period("M")
        monthly = (
            df.groupby(["Transmission", "month"], as_index=False)["num"].sum()
        )
        monthly["date"] = monthly["month"].dt.to_timestamp(how="start") + pd.offsets.Day(
            14
        )
        monthly["Epiweek"] = monthly["month"].astype(str)
        return monthly[["Epiweek", "num", "Transmission", "date"]]

    raise ValueError("`time_unit` must be either 'week' or 'month'.")


def plot_transmission_time_series(
    df: pd.DataFrame,
    *,
    ax: plt.Axes | None = None,
    title: str | None = None,
    colors: Mapping[str, str] | None = None,
) -> plt.Axes:
    """
    Generate the line/point plot from `epidemicInTexas.R`.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 3.5))

    color_map = {
        "importation": IMPORT_COLOR,
        "localtransmission": LOCAL_COLOR,
        "exportation": EXPORT_COLOR,
    }
    if colors:
        color_map.update(colors)

    pretty_labels = {
        "importation": "Importation",
        "localtransmission": "Local Transmission",
        "exportation": "Export",
    }

    for raw_label in ("importation", "localtransmission", "exportation"):
        color = color_map[raw_label]
        subset = df[df["Transmission"] == raw_label]
        ax.plot(
            subset["date"],
            subset["num"],
            marker="o",
            linewidth=1.5,
            markersize=4,
            label=pretty_labels[raw_label],
            color=color,
        )
    ax.set_ylabel("Transmission counts", size=16)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    if title:
        ax.set_title(title)
    ax.legend(title="Transmission Type", loc="upper left", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", alpha=0.2)
    return ax


def build_branch_color_lookup(
    tree: bt.tree,
    trait: str,
    *,
    location: str = "Texas",
) -> Mapping[int, str]:
    """
    Recreate the color-coding of branches used in `epidemicInTexas.R`.
    """
    colors = {}
    for branch in tree.Objects:
        parent = branch.parent
        origin = parent.traits.get(trait) if parent else None
        destination = branch.traits.get(trait)
        if origin == location and destination == location:
            colors[branch.index] = LOCAL_COLOR
        elif origin != location and destination == location:
            colors[branch.index] = IMPORT_COLOR
        elif origin == location and destination != location:
            colors[branch.index] = EXPORT_COLOR
        else:
            colors[branch.index] = DEFAULT_COLOR
    return colors


def north_south_descendant_sequences(
    tree_path: str | Path,
    trait: str,
    target: str,
    *,
    absolute_time: bool = True,
    sort_branches: bool = True,
) -> pd.Series:
    """
    Count unique descendant tips for north/south transitions involving ``target``.

    Each qualifying transition contributes all descendant tip names beneath the
    branch where it occurs. Tips are only counted once per category.
    """
    tree = bt.loadNexus(
        str(tree_path),
        absoluteTime=absolute_time,
        sortBranches=sort_branches,
    )

    tip_sets: Dict[str, Set[str]] = {
        "north_import": set(),
        "south_import": set(),
        "north_export": set(),
        "south_export": set(),
    }

    for branch in tree.Objects:
        parent = branch.parent
        if parent is None:
            continue
        origin = parent.traits.get(trait) if parent.traits else None
        destination = branch.traits.get(trait) if branch.traits else None
        category = _north_south_transition_category(origin, destination, target)
        if category is None:
            continue
        for tip_name in getattr(branch, "leaves", []):
            if tip_name:
                tip_sets[category].add(tip_name)

    counts = {key: len(value) for key, value in tip_sets.items()}
    return pd.Series(counts, name="descendant_tips")


__all__ = [
    "DEFAULT_COLOR",
    "tree_to_table",
    "count_spatial_transmission_linkages",
    "accumulated_indicator",
    "transmission_time_series",
    "plot_transmission_time_series",
    "build_branch_color_lookup",
    "decimal_year_to_date",
    "epiweek_to_date",
    "NORTH_COUNTRIES",
    "SOUTH_COUNTRIES",
    "north_south_descendant_sequences",
]
