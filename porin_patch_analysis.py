import os
from pathlib import Path
import numpy as np
import pandas as pd
import json
import yaml

import matplotlib.pyplot as plt
import seaborn as sns

from shapely.geometry import Point, MultiPoint, Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely.plotting import plot_polygon, plot_line
from shapely import simplify, delaunay_triangles

def analyse_trimer_free_regions(coords, roi, max_edge_separation):
    """Identify and characterize OMP-trimer-free regions
    using Delaunay triangles method (may overestimate patch area)
    Parameters:
    - coords: np.ndarray
        Shape (N, 2) array of [Y,X] coordinates of trimers (nm)
    - roi:
    - buffer_radius: float
        Radius around trimer detections in nm used as part of LPS patch
        definition
    """
    # Switch to [X,Y] coords
    coords_xy = coords.copy()
    coords_xy[:,[1,0]] = coords[:,[0,1]]

    # Create triangles from trimer vertices
    points = MultiPoint([Point(coords_xy[row,:]) for row in range(coords_xy.shape[0])])
    triangles = delaunay_triangles(points)
    
    # Remove triangles with edges larger than max_edge_separation
    valid_triangles = []
    for triangle in triangles.geoms:
        tri_coords = list(triangle.exterior.coords)[:-1] # Remove duplicate last point
        edge_lengths = [
            Point(tri_coords[0]).distance(Point(tri_coords[1])),
            Point(tri_coords[1]).distance(Point(tri_coords[2])),
            Point(tri_coords[2]).distance(Point(tri_coords[0]))
        ]

        if max(edge_lengths) <= max_edge_separation:
            valid_triangles.append(triangle)

    # Make area from trimer triangles
    occupied_union = unary_union(valid_triangles)
    
    # LPS patches are the complement
    lps_patches = roi.difference(occupied_union)
    
    # Extract individual patches
    patches = []
    if lps_patches.geom_type == 'Polygon':
        patches = [lps_patches]
    elif lps_patches.geom_type == 'MultiPolygon':
        patches = list(lps_patches.geoms)
    
    # Analyze patch properties
    patch_data = []
    roi_boundary = roi.boundary
    for patch in patches:
        # Ignore patch if it is on edge of image (roi)
        if patch.touches(roi_boundary):
            pass
        else:
            # Simplify area
            # patch = simplify(patch, tolerance=1)

            area = patch.area
            perimeter = patch.length
            centroid = patch.centroid
            
            # Shape metrics
            circularity = 4 * np.pi * area / (perimeter ** 2) if perimeter > 0 else 0
            
            # Approximate equivalent diameter
            equiv_diameter = 2 * np.sqrt(area / np.pi)
            
            patch_data.append({
                'geometry': patch,
                'area': area,
                'perimeter': perimeter,
                'centroid': (centroid.x, centroid.y),
                'circularity': circularity,
                'equiv_diameter': equiv_diameter
            })
    
    # Sort by area
    patch_data.sort(key=lambda x: x['area'], reverse=True)

    lps_patches = MultiPolygon([patch["geometry"] for patch in patch_data])
    
    return patch_data, lps_patches, occupied_union

if __name__ == "__main__":
    cwd = Path(__file__).parent.resolve() # Script directory
    results_dir = os.path.join(cwd, "results") # Directory to dump results
    localisation_dir = os.path.join(cwd, "trimer-localisations") # Directory to load localisation x-y coords from
    os.makedirs(results_dir, exist_ok=True)

    # Grab hyperparams
    param_path = os.path.join(cwd, "image-analysis-hyperparams.yaml")
    with open(param_path) as f:
        params = yaml.safe_load(f)
    params.update({"hyperparameter_file_path": param_path})

    # Paths to trimer localisations in different CSV files
    all_data = {
        "WT": {
            "base_dir": os.path.join(localisation_dir, "WT"),
            "localisation_csvs": [
                "MG_C1_save-2025.07.21-15.21.04.930_trimers_yx.csv",
                "MG_C2_save-2025.07.21-15.47.02.798_trimers_yx.csv",
                "MG_C3_save-2025.07.21-15.52.28.580_trimers_yx.csv",
                "MG1_C1_save-2025.10.02-12.14.40.128_trimers_yx.csv",
                "MG1_C2_save-2025.10.02-12.22.29.107_trimers_yx.csv",
                "MG1_C3_save-2025.10.02-12.49.44.440_trimers_yx.csv",
                "MG3_C1_save-2025.10.03-14.24.42.873_trimers_yx.csv",
                "MG3_C2_save-2025.10.03-15.06.23.276_trimers_yx.csv",
                "MG3_C3_save-2025.10.03-15.19.52.655_trimers_yx.csv",
            ],
        },
        "WbbL": {
            "base_dir": os.path.join(localisation_dir, "WbbL"),
            "localisation_csvs": [
                "G_WbbL_1_save-2025.07.21-16.31.11.637_0_0_trimers_yx.csv",
                "G_WbbL_2_save-2025.07.21-16.37.02.142_0_0_trimers_yx.csv",
                "G_WbbL_C3_save-2025.07.21-16.52.10.880_trimers_yx.csv",
                "WbbL1_C1_save-2025.10.03-12.37.56.169_trimers_yx.csv",
                "WbbL1_C2save-2025.10.03-12.57.56.038_trimers_yx.csv",
                "WbbL1_C3_save-2025.10.03-13.06.50.330_trimers_yx.csv",
                "WbbL2_C1_save-2025.10.03-13.23.23.848_trimers_yx.csv",
                "WbbL2_C2_save-2025.10.03-13.38.45.750_trimers_yx.csv",
                "WbbL2_C3_save-2025.10.03-13.58.04.246_trimers_yx.csv",
            ],
        },
    }

    patch_params = params["patch_params"]
    patch_infos = []
    for strain, strain_data in all_data.items():
        for i, csv_fname in enumerate(strain_data["localisation_csvs"]):

            # Locate and load trimer localisations from individual CSV files
            localisation_path = os.path.join(strain_data["base_dir"], csv_fname)
            with open(localisation_path, "r") as f: # First, grab image extent from CSV comment
                for i, line in enumerate(f):
                    if i == 0:
                        image_name = line[line.find(": ")+2:-1]
                    if i == 2:
                        line = line[1:]
                        line = line.replace("(", "[")
                        line = line.replace(")", "]")
                        (xmin, xmax, ymin, ymax) = json.loads(line)
            trimers_yx_nm = pd.read_csv(localisation_path, sep=",", usecols=["x_nm", "y_nm"], comment="#").to_numpy()

            roi = Polygon([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])
            patch_data, lps_patches, occupied_union = analyse_trimer_free_regions(trimers_yx_nm, roi, patch_params["max_edge_separation"])

            # Store data for writing to CSV later
            for idx, patch in enumerate(patch_data):
                x, y = patch["centroid"]

                patch_info = {
                    "strain": strain,
                    "image_name": image_name,
                    "patch_index": idx,
                    "patch_centroid_x_nm": x,
                    "patch_centroid_y_nm": y,
                    "patch_area_nm2": patch["area"],
                    "patch_circularity": patch["circularity"],
                }
                patch_infos.append(patch_info)

    # Convert patch data to pandas dataframe
    patch_info = pd.DataFrame(patch_infos)

    # Plot data from each cell type
    fig, axes = plt.subplots(ncols=2)
    sns.violinplot(patch_info, x="strain", y="patch_area_nm2", color="magenta", ax=axes[0])
    sns.violinplot(patch_info, x="strain", y="patch_circularity", color="cyan", ax=axes[1])

    for ax in axes:
        ax.set_xlabel("Strain")
    axes[0].set_ylabel("Patch area (nm$^3$)")
    axes[1].set_ylabel("Patch circularity")
    plt.suptitle("LPS Patch Characteristics")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "summary_histogram.png"), dpi=450)
    plt.show()

    # Save dataframe as a csv
    patch_info.to_csv(os.path.join(results_dir, "all_patch_info.csv"))