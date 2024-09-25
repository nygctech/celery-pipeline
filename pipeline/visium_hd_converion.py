import argparse
import random
import os
from collections import deque
import numpy as np
from scipy.sparse import coo_matrix
import pandas as pd
import shapely
import shapely.affinity
import shapely.geometry
import geojson


def rebin_positions(positions, bin_size, in_tissue_pct=0.5):
    if isinstance(positions, str):
        old_tissue_positions = pd.read_parquet(positions)
    elif isinstance(positions, pd.DataFrame):
        old_tissue_positions = positions.copy()
    else:
        raise ValueError("Argument `positions` must be a string or a DataFrame")

    # barcode indices are of the form s_[binsize]um_00000_00000-1
    old_bin_size = int(old_tissue_positions.loc[0, "barcode"].split("_")[1].strip("um"))

    assert (
        bin_size >= old_bin_size
    ), f"New bin size {bin_size} must be greater than the original bin size: {old_bin_size}"
    assert (
        bin_size % old_bin_size == 0
    ), f"New bin size {bin_size} must be a multiple of original bin size: {old_bin_size}"

    if old_bin_size == bin_size:
        return positions

    print(f"Rebinning spots based on bin size {bin_size}")

    max_old_row = old_tissue_positions["array_row"].max()
    max_old_col = old_tissue_positions["array_col"].max()

    old_tissue_positions = old_tissue_positions.set_index(["array_row", "array_col"])

    # could do a for loop

    bin_ratio = bin_size // old_bin_size
    max_new_row = max_old_row // bin_ratio
    max_new_col = max_old_col // bin_ratio

    data_list = []
    for r in range(max_new_row):
        start_r = r * bin_ratio
        end_r = start_r + bin_ratio - 1
        for c in range(max_new_col):
            start_c = c * bin_ratio
            end_c = start_c + bin_ratio - 1

            barcode = (
                f"s_{str(bin_size).zfill(3)}um_{str(r).zfill(5)}_{str(c).zfill(5)}-1"
            )
            # loc with slice is inclusive of the stop range...
            inner_spots = old_tissue_positions.loc[
                (slice(start_r, end_r), slice(start_c, end_c)), :
            ]

            in_tissue = int(inner_spots["in_tissue"].mean() >= in_tissue_pct)
            pxl_row = inner_spots["pxl_row_in_fullres"].mean()
            pxl_col = inner_spots["pxl_col_in_fullres"].mean()

            data_list.append((barcode, in_tissue, r, c, pxl_row, pxl_col))

    new_tissue_positions = pd.DataFrame(
        data_list,
        columns=[
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ],
    )

    return new_tissue_positions


def flood_fill_annotations(df):
    # city block style flood fill/expand from annotations on the outline of regions

    df = df.copy()

    df["labeled_in_tissue"] = df["in_tissue"] & (df["AARs"] != "")
    # convert to numpy array
    data = df["labeled_in_tissue"].to_numpy()
    row = df["array_row"].to_numpy()
    col = df["array_col"].to_numpy()
    arr = coo_matrix((data, (row, col)), (row.max() + 1, col.max() + 1)).toarray()

    # find outlines of annotated regions
    ones = arr == 1

    # Check for adjacent 0s by shifting the array and performing logical OR
    adjacent_to_zero = (
        (np.roll(arr, 1, axis=0) == 0)  # Shift up
        | (np.roll(arr, -1, axis=0) == 0)  # Shift down
        | (np.roll(arr, 1, axis=1) == 0)  # Shift left
        | (np.roll(arr, -1, axis=1) == 0)  # Shift right
    )
    border = ones & adjacent_to_zero
    border_row_col = np.where(border)

    # list of already annotated, but on the border, spots
    to_visit = deque([(r, c) for (r, c) in zip(border_row_col[0], border_row_col[1])])
    visited = set()

    # index by rows and columns now
    df = df.set_index(["array_row", "array_col"])

    # bfs/floodfill type algorithm
    while len(to_visit) > 0:
        r, c = to_visit.popleft()

        # in case we've already visited
        if (r, c) in visited:
            continue
        else:
            visited.add((r, c))
        current_aar = df.loc[(r, c), "AARs"]
        neighbors = [(r + 1, c), (r - 1, c), (r, c + 1), (r, c - 1)]

        # print(len(to_visit))
        for neighbor_r, neighbor_c in neighbors:
            if (neighbor_r, neighbor_c) not in df.index:
                continue

            row = df.loc[(neighbor_r, neighbor_c)]
            # set the aar annotation and add it
            if row["in_tissue"] and row["AARs"] == "":
                df.at[(neighbor_r, neighbor_c), "AARs"] = current_aar
                to_visit.append((neighbor_r, neighbor_c))

    return df


def get_polygons_aars(geoj, x_translation=0, y_translation=0, scaling=1):
    """
    Return list of pairs of shapely polygons and respsective aars
    """
    l = []
    for feat in geoj["features"]:
        coords = feat["geometry"]["coordinates"]
        polygon = shapely.geometry.Polygon(coords)
        polygon = shapely.affinity.translate(
            polygon, xoff=x_translation, yoff=y_translation
        )
        polygon = shapely.affinity.scale(
            polygon, xfact=scaling, yfact=scaling, origin=(0, 0)
        )
        aar = feat["properties"]["annotation"]
        l.append((polygon, aar))
    return l


def get_pt_aar(pt, polygons_aars, seed=0):
    """
    Get AAR of region pt is within. Return random if within multiple.
    Return empty string if in none
    """
    random.seed(hash((pt.x, pt.y, seed)))
    aars = []
    for polygon, aar in polygons_aars:
        if polygon.contains(pt):
            aars.append(aar)

    return random.choice(aars) if len(aars) > 0 else ""


def get_spatial_aar_df(
    positions,
    annotation_path,
    x_translation=0,
    y_translation=0,
    scaling=1,
    bin_size=None,
    in_tissue_pct=0.5,
):
    if isinstance(positions, str):
        df = pd.read_parquet(positions)
    elif isinstance(positions, pd.DataFrame):
        df = positions.copy()
    else:
        raise ValueError("Argument `positions` must be a string or a DataFrame")

    if bin_size is not None:
        df = rebin_positions(df, bin_size=bin_size, in_tissue_pct=in_tissue_pct)

    with open(annotation_path, "rb") as f:
        geoj = geojson.load(f)

    polygons_aars = get_polygons_aars(geoj, x_translation, y_translation, scaling)

    df["pt"] = [
        shapely.geometry.Point(x, y)
        for x, y in df[["pxl_col_in_fullres", "pxl_row_in_fullres"]].itertuples(
            index=False, name=None
        )
    ]

    print("Getting barcode annotations from json regions")
    get_aar = lambda pt: get_pt_aar(pt, polygons_aars=polygons_aars, seed=0)
    df["AARs"] = df["pt"].apply(get_aar)

    print("Assigning unannotated spots to nearby contigous regions")
    df = flood_fill_annotations(df)
    return df


def get_bin_options(spaceranger_path):
    binned_dir = os.path.join(spaceranger_path, "outs/binned_outputs")

    paths = [p for p in os.listdir(binned_dir) if p.startswith("square_")]
    bins = [int(p.split("_")[-1].strip("um")) for p in paths]

    return bins


def process_visium_data(
    spaceranger_path,
    annotation_path,
    target_bin_size,
    out_path,
    save_parquet,
    xy_translation,
    scalefactor,
):

    bin_options = get_bin_options(spaceranger_path)

    assert (
        len(bin_options) > 0
    ), f"Could not find any binned outputs in spaceranger path {spaceranger_path}"

    largest_divisible_bin = 0
    for bin in bin_options:
        if target_bin_size % bin == 0 and bin > largest_divisible_bin:
            largest_divisible_bin = bin

    assert (
        largest_divisible_bin > 0
    ), f"None of the binned outputs in spaceranger, {bin_options}, are a multiple of target_bin_size {target_bin_size}"

    positions_path = os.path.join(
        spaceranger_path,
        f"outs/binned_outputs/square_{str(largest_divisible_bin).zfill(3)}um/spatial/tissue_positions.parquet",
    )

    x_translation, y_translation = xy_translation
    processed_df = get_spatial_aar_df(
        positions_path,
        annotation_path,
        x_translation,
        y_translation,
        scalefactor,
        target_bin_size,
    )
    processed_df = processed_df.reset_index()

    aar_df = processed_df[processed_df["in_tissue"] == 1][["barcode", "AARs"]].rename(
        {"barcode": "Barcode"}, axis=1
    )
    aar_df.to_csv(out_path, index=False)

    if save_parquet:
        path = f"custom_bin_{str(target_bin_size).zfill(3)}.parquet"
        print(f"Saving spatial file to {path}")
        parquet_df = processed_df[
            [
                "array_row",
                "array_col",
                "barcode",
                "in_tissue",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres",
            ]
        ]
        parquet_df.to_parquet(path, index=False)


def main():
    parser = argparse.ArgumentParser(description="CLI for processing Visium HD data")

    # Required arguments
    parser.add_argument(
        "spaceranger_path", type=str, help="Path to the Spaceranger output files"
    )
    parser.add_argument(
        "annotation_path",
        type=str,
        help="Path to the json annotation file exported by Celery",
    )
    parser.add_argument(
        "target_bin_size",
        type=int,
        help="Target bin size, will automatically rebin if size does not exist",
    )

    # Optional arguments
    parser.add_argument(
        "--out_path",
        "-o",
        type=str,
        default="AARs.csv",
        help='Path to csv file with barcodes and AARs (default "AARs.csv")',
    )
    parser.add_argument(
        "--save_parquet",
        "-s",
        action="store_true",
        help="Save the rebinned data as a parquet file",
    )
    parser.add_argument(
        "--xy_translation",
        "-t",
        type=float,
        nargs=2,
        metavar=("x", "y"),
        default=(0, 0),
        help="Translation applied to annotations in x and y directions (default 0 0)",
    )
    parser.add_argument(
        "--scalefactor",
        "-f",
        type=float,
        default=1,
        help="Scale factor applied to annotations after translation (default 1)",
    )

    args = parser.parse_args()

    # Call the function to process the data
    process_visium_data(
        args.spaceranger_path,
        args.annotation_path,
        args.target_bin_size,
        args.out_path,
        args.save_parquet,
        args.xy_translation,
        args.scalefactor,
    )


if __name__ == "__main__":
    main()
