import sys
import inspect
from pathlib import Path

sys.path.insert(0, str(Path("Scripts").resolve()))
from transmission_count import (
    tree_to_table,
    transmission_time_series,
    accumulated_indicator,
    get_window_size,
    count_spatial_transmission_linkages,
    north_south_descendant_sequences,
)


def summarise_tree(
    tree_path: Path,
    trait: str,
    location: str,
    *,
    output_prefix: str,
    include_descendants: bool = True,
) -> dict[str, Path]:
    """
    Run the transmission summary workflow for a single trait-annotated tree.
    """
    window_size = get_window_size(tree_path)
    branch_cutoff = window_size  # days, mirrors the R pipeline
    print(f"[{output_prefix}] Using branch cutoff of {branch_cutoff} days")

    migration_table = tree_to_table(tree_path, trait, branch_cutoff)
    csv_path = Path(f"{output_prefix}_migration_table.csv")
    migration_table.to_csv(csv_path, index=False)
    print(migration_table)

    series = transmission_time_series(migration_table, location)
    print(series.groupby("Transmission")["num"].sum())  # totals per category
    print(series[series["num"] > 0])  # first non-zero weeks

    lis = accumulated_indicator(migration_table, [location], indicator="LIS")
    sss = accumulated_indicator(migration_table, [location], indicator="SSS")
    print(f"[{output_prefix}] LIS: {lis[location]}")
    print(f"[{output_prefix}] SSS: {sss[location]}")

    north_south_counts = count_spatial_transmission_linkages(
        migration_table, "north2south", location
    )
    counts_path = Path(f"{output_prefix}_north_south_counts.csv")
    north_south_counts.to_csv(counts_path)
    print(f"[{output_prefix}] North/South epiweek totals:")
    print(north_south_counts.sum().to_dict())

    outputs: dict[str, Path] = {
        "migration_table": csv_path,
        "north_south_counts": counts_path,
    }

    if include_descendants:
        descendant_counts = north_south_descendant_sequences(tree_path, trait, location)
        descendants_path = Path(f"{output_prefix}_north_south_descendants.csv")
        descendant_counts.to_frame(name="descendant_tips").to_csv(descendants_path)
        print(f"[{output_prefix}] North/South descendant tips:")
        print(descendant_counts.to_dict())
        outputs["north_south_descendants"] = descendants_path

    return outputs


def main() -> None:
    location = "Colombia"

    # Full dengue genome tree (Nextstrain build)
    full_tree_path = Path("Revision/dengue_denv2_genome_dates.nexus")
    full_trait = "country"
    summarise_tree(
        full_tree_path,
        full_trait,
        location,
        output_prefix="full_tree",
        include_descendants=True,
    )

    # Subsampled BEAST MCC tree
    beast_tree_path = Path("Revision/newMCCTree_F1112_Markov.mcc.tree")
    beast_trait = "location.rate"
    results_beast = summarise_tree(
        beast_tree_path,
        beast_trait,
        location,
        output_prefix="beast_tree",
        include_descendants=False,
    )

    # Ensure dyad-level migration table is available for the BEAST tree if supported.
    try:
        sig = inspect.signature(tree_to_table)
        if any(
            p.name in ("mode", "return_dyads", "per_transition")
            for p in sig.parameters.values()
        ):
            try:
                migration_table_beast = tree_to_table(
                    beast_tree_path, beast_trait, get_window_size(beast_tree_path), mode="dyads"
                )
            except TypeError:
                try:
                    migration_table_beast = tree_to_table(
                        beast_tree_path,
                        beast_trait,
                        get_window_size(beast_tree_path),
                        return_dyads=True,
                    )
                except TypeError:
                    migration_table_beast = tree_to_table(
                        beast_tree_path,
                        beast_trait,
                        get_window_size(beast_tree_path),
                        per_transition=True,
                    )
            dyad_path = Path("beast_tree_migration_table_dyads.csv")
            migration_table_beast.to_csv(dyad_path, index=False)
            results_beast["migration_table_dyads"] = dyad_path
            print("[beast_tree] Rebuilt migration_table in DYAD mode.")
            print(migration_table_beast.head())
    except Exception as exc:
        print(f"[beast_tree] Warning: could not rebuild migration table in dyad mode: {exc}")


if __name__ == "__main__":
    main()
