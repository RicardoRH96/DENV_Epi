"""
Python translation of `epidemicInTexas.R`.

Run this module as a script to reproduce the migration summaries and plots
using the helpers defined in `transmission_count.py`.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Mapping

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch
from matplotlib.legend_handler import HandlerBase

import baltic as bt

from transmission_count import (
    DEFAULT_COLOR,
    build_branch_color_lookup,
    get_window_size,
    plot_transmission_time_series,
    transmission_time_series,
    tree_to_table,
    decimal_year_to_date,
)


# --- put these near the top of your module ---
DEFAULT_COLOR   = "#D3D3D3"  # non-location node (light gray)
LOC_COLOR       = "#e09f3e"  # location node (Texas blue)
BG_LINE         = "#BDBDBD"  # global/epidemic background (gray line)
IMPORT_LINE     = "#22C55E"  # import (green line)
LOCAL_LINE      = LOC_COLOR  # local transmission (blue line)
EXPORT_LINE     = "#9c6ade"  # export (soft purple line)
# ---------------------------------------------

def _edge_handle(line_color, left_color, right_color, label, lw=2, msize=6):
    """Return a (left_marker, line, right_marker) tuple for a legend entry."""
    left = Line2D([0], [0], linestyle="none", marker="o",
                  markerfacecolor=left_color, markeredgecolor="none",
                  markersize=msize)
    line = Line2D([0, 1], [0, 0], color=line_color, linewidth=lw)
    right = Line2D([0], [0], linestyle="none", marker="o",
                   markerfacecolor=right_color, markeredgecolor="none",
                   markersize=msize)
    left.set_zorder(3)
    right.set_zorder(3)
    line.set_zorder(2)
    left._legend_label = label  # attach label to tuple via first element
    return (left, line, right)


class EdgeLegendHandler(HandlerBase):
    """
    Custom legend handler that draws point-line-point entries with slight overlap.
    """

    def __init__(self, pad: float = 0.15):
        super().__init__()
        self.pad = pad

    def create_artists(
        self,
        legend,
        orig_handle,
        xdescent,
        ydescent,
        width,
        height,
        fontsize,
        trans,
    ):
        left_handle, line_handle, right_handle = orig_handle
        pad = self.pad * fontsize
        x_start = xdescent + pad
        x_end = xdescent + width - pad
        y = -ydescent + height / 2

        line_artist = Line2D(
            [x_start, x_end],
            [y, y],
            color=line_handle.get_color(),
            linewidth=line_handle.get_linewidth(),
            linestyle=line_handle.get_linestyle(),
        )
        line_artist.set_transform(trans)
        line_artist.set_zorder(2)

        left_artist = Line2D(
            [x_start],
            [y],
            marker=left_handle.get_marker(),
            linestyle="none",
            markerfacecolor=left_handle.get_markerfacecolor(),
            markeredgecolor=left_handle.get_markeredgecolor(),
            markersize=left_handle.get_markersize(),
        )
        left_artist.set_transform(trans)
        left_artist.set_zorder(3)

        right_artist = Line2D(
            [x_end],
            [y],
            marker=right_handle.get_marker(),
            linestyle="none",
            markerfacecolor=right_handle.get_markerfacecolor(),
            markeredgecolor=right_handle.get_markeredgecolor(),
            markersize=right_handle.get_markersize(),
        )
        right_artist.set_transform(trans)
        right_artist.set_zorder(3)

        return [line_artist, left_artist, right_artist]


class HeaderLegendHandler(HandlerBase):
    """
    Legend handler to suppress the marker for section headers
    so their text can sit near the edge of the legend box.
    """

    def create_artists(
        self,
        legend,
        orig_handle,
        xdescent,
        ydescent,
        width,
        height,
        fontsize,
        trans,
    ):
        pad = 0.05
        x = xdescent - width + pad
        y = -ydescent + height / 2
        marker = Line2D(
            [x],
            [y],
            marker=orig_handle.get_marker(),
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor="none",
            markersize=orig_handle.get_markersize(),
        )
        marker.set_transform(trans)
        return [marker]
def plot_transition_tree(
    tree: bt.tree,
    color_lookup: Mapping[int, str],
    *,
    trait: str,
    location: str,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 12), facecolor="w")

    def color_func(branch: bt.leaf | bt.node) -> str:
        return color_lookup.get(branch.index, DEFAULT_COLOR)

    x_attr = lambda branch: getattr(branch, "absoluteTime", getattr(branch, "height", 0.0))

    tree.plotTree(ax, colour=color_func, x_attr=x_attr, width=1.5)
    tree.plotPoints(ax, colour=color_func, x_attr=x_attr, size=20, zorder=5)
    ax.set_ylabel("")
    ax.set_yticks([])
    #Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)


    most_recent_year = getattr(tree, "mostRecent", None)
    if most_recent_year is not None:
        most_recent_date = decimal_year_to_date(most_recent_year)
        print("Most recent collection date:", most_recent_date.date())

    if most_recent_year is not None:
        earliest = min(getattr(branch, "absoluteTime", getattr(branch, "height", most_recent_year)) for branch in tree.Objects)
        years = range(int(earliest) - 1, int(most_recent_year) + 2)
        for year in years:
            if year % 2 == 0:
                ax.axvspan(year, year + 1, fc="k", ec="none", alpha=0.03, zorder=0)
    ax.tick_params(axis="x", labelsize=16)

    # ---------- Legend ----------
    header_nodes = Line2D([0], [0], marker="o", linestyle="none",
                          markerfacecolor="none", markeredgecolor="none", markersize=8,
                          label="Node Trait")
    header_links = Line2D([0], [0], marker="o", linestyle="none",
                          markerfacecolor="none", markeredgecolor="none", markersize=8,
                          label="Spatial Transmission Linkage")
    handles = [
        header_nodes,
        Line2D([0], [0], marker="o", linestyle="none",
               markerfacecolor=DEFAULT_COLOR, markeredgecolor="none", markersize=8,
               label=f"non-{location}"),
        Line2D([0], [0], marker="o", linestyle="none",
               markerfacecolor=LOC_COLOR, markeredgecolor="none", markersize=8,
               label=location),
        header_links,
        _edge_handle(BG_LINE, BG_LINE, BG_LINE, "Americas D2_II.F.1.1.2 Background"),
        _edge_handle(IMPORT_LINE, DEFAULT_COLOR, LOC_COLOR, "Import"),
        _edge_handle(LOCAL_LINE, LOC_COLOR, LOC_COLOR, "Local Transmission"),
        _edge_handle(EXPORT_LINE, LOC_COLOR, DEFAULT_COLOR, "Export"),
    ]

    labels = [
        "Node Trait",
        f"non-{location}",
        location,
        "Spatial Transmission Linkage",
        "II.F.1.1.2 Background",
        "Import",
        "Local Transmission",
        "Export",
    ]

    leg = ax.legend(
        handles,
        labels,
        loc="upper left",
        frameon=False,
        handlelength=2.4,
        handletextpad=0.8,
        labelspacing=0.8,
        handler_map={
            tuple: EdgeLegendHandler(),
            header_nodes: HeaderLegendHandler(),
            header_links: HeaderLegendHandler(),
        },
    )

    for text, handle in zip(leg.get_texts(), leg.legend_handles):
        label = text.get_text()
        if label in {"Node Trait", "Spatial Transmission Linkage"}:
            text.set_fontweight("bold")
            text.set_x(text.get_position()[0] - 35)

    return ax


def run_pipeline(
    tree_path: Path,
    trait: str,
    location: str,
    branch_cutoff: float | None,
    *,
    output_prefix: Path | None = None,
) -> None:
    """
    Execute the migration table aggregation and produce two plots.
    """
    if branch_cutoff is None:
        try:
            branch_cutoff = get_window_size(tree_path)
        except FileNotFoundError:
            branch_cutoff = 15.0
    migration_table = tree_to_table(tree_path, trait, branch_cutoff)
    series = transmission_time_series(migration_table, location, time_unit="month")

    tree = bt.loadNexus(str(tree_path), absoluteTime=True, sortBranches=False)
    tree.sortBranches(descending=False)
    color_lookup = build_branch_color_lookup(tree, trait, location=location)

    fig, (ax_tree, ax_ts) = plt.subplots(
        2,
        1,
        sharex=False,
        figsize=(14, 12),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.1},
        facecolor="w",
    )

    plot_transition_tree(tree, color_lookup, trait=trait, location=location, ax=ax_tree)
    plot_transmission_time_series(series, ax=ax_ts, title="")
    ax_ts.set_xlabel("Time (decimal year)")

    if output_prefix is not None:
        pdf_path = output_prefix.with_suffix(".combined.pdf")
        jpeg_path = output_prefix.with_suffix(".combined.jpeg")
        fig.savefig(pdf_path, bbox_inches="tight")
        fig.savefig(jpeg_path, bbox_inches="tight", dpi=600)
    else:
        plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Translate epidemicInTexas.R pipeline into Python.",
    )
    parser.add_argument(
        "tree",
        type=Path,
        help="Path to the time-scaled NEXUS tree.",
    )
    parser.add_argument(
        "--trait",
        default="country",
        help="Discrete trait to track (default: country).",
    )
    parser.add_argument(
        "--location",
        default="Colombia",
        help="Focal location (default: Colombia).",
    )
    parser.add_argument(
        "--branch-cutoff",
        type=float,
        default=None,
        help="Maximum branch length in days to keep; defaults to tree window size.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=None,
        help="When provided, save the figures using <prefix>.timeseries.pdf and <prefix>.tree.pdf.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_pipeline(
        args.tree,
        args.trait,
        args.location,
        args.branch_cutoff,
        output_prefix=args.output_prefix,
    )
