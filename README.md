# bacteria-imaging

Analysis code accompanying the manuscript:

> **Nanometer-resolution imaging of living bacteria across species and cellular orientations**
>
> Anna Scheeder, Joanna Szczepaniak, Yulianna Koziy, Naomi Mara Claro De Oliveira, William Trewby, Renata Kaminska, Jeremy Brown, Colin Kleanthous, Bart W. Hoogenboom
>
> _Nature Communications_ (in revision).

This repository contains the code and the derived data needed to reproduce the
outer-membrane-protein (OMP) trimer-free ("LPS patch") analysis reported in Fig. 3 of the
manuscript, together with the trimer localisations extracted from the AFM images
and the resulting patch statistics.

---

## Repository contents

| Path | Description |
| --- | --- |
| `porin_patch_analysis.py` | Main analysis script. Identifies and characterises OMP-trimer-free regions from trimer localisations. |
| `image-analysis-hyperparams.yaml` | All tunable analysis parameters (currently the Delaunay edge-length threshold). |
| `trimer-localisations/WT/*.csv` | Trimer centre coordinates for the wild-type strain, one CSV per AFM image (9 images). |
| `trimer-localisations/WbbL/*.csv` | Trimer centre coordinates for the WbbL strain, one CSV per AFM image (9 images). |
| `results/all_patch_info.csv` | Per-patch measurements produced by the script (526 patches). |
| `results/summary_histogram.png` | Violin plots of patch area and circularity by strain. |
| `requirements.txt`, `environment.yml` | Pinned software environments (pip and conda respectively). |
| `LICENSE` | MIT licence. |
| `CITATION.cff` | Citation metadata for this repository. |

The AFM image stacks from which the trimer localisations were extracted are
deposited separately; see **Data availability** below.

---

## System requirements

### Software dependencies and versions

The analysis is a pure-Python script. It has been tested with the following
package versions:

| Package | Tested version | Minimum required |
| --- | --- | --- |
| Python | 3.12.5 | 3.9 |
| numpy | 2.2.3 | 1.20 |
| pandas | 2.2.3 | 1.3 |
| shapely | 2.1.1 | 2.0 (`delaunay_triangles` requires Shapely ≥ 2.0) |
| matplotlib | 3.10.0 | 3.5 |
| seaborn | 0.13.2 | 0.12 |
| PyYAML | 6.0.2 | 5.4 |

No other software is required. The code uses no compiled extensions of its own
and no GPU libraries.

### Operating systems tested

- Windows 11 Enterprise (build 26200) — primary development and analysis platform.

The code contains no platform-specific calls (paths are constructed with
`os.path`/`pathlib`), and is expected to run unchanged on Linux and macOS with
the dependency versions listed above.

### Hardware requirements

No non-standard hardware is required. A standard desktop or laptop computer is
sufficient; the full analysis uses well under 1 GB of RAM and a single CPU core.
Timings quoted below were measured on an Intel Core i7-13800H (14 cores) with
32 GB RAM.

---

## Installation guide

### Option A — conda (recommended)

```bash
git clone https://github.com/hoogenboom-lab/bacteria-imaging.git
cd bacteria-imaging
conda env create -f environment.yml
conda activate porin-patch-analysis
```

### Option B — pip / virtualenv

```bash
git clone https://github.com/hoogenboom-lab/bacteria-imaging.git
cd bacteria-imaging
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

**Typical install time:** 2–5 minutes on a normal desktop computer with a
broadband connection (dominated by downloading numpy, pandas, matplotlib and
shapely wheels). Cloning the repository takes a few seconds.

---

## Demo

The trimer localisations included in `trimer-localisations/` are the data used
in the manuscript and double as the demo dataset — no separate download or
configuration is needed.

**Run:**

```bash
python porin_patch_analysis.py
```

**Expected output** (written to `results/`, overwriting the copies shipped in
the repository; the regenerated CSV reproduces the shipped one exactly):

- `all_patch_info.csv` — one row per detected trimer-free patch, with columns
  `strain`, `image_name`, `patch_index`, `patch_centroid_x_nm`,
  `patch_centroid_y_nm`, `patch_area_nm2`, `patch_circularity`.
- `summary_histogram.png` — a two-panel figure of patch area and patch
  circularity by strain.

**Expected numerical result** for the included data and the default
`max_edge_separation: 27.5` nm:

| Strain | Images | Patches detected | Median patch area (nm²) | Median circularity |
| --- | --- | --- | --- | --- |
| WT | 9 | 282 | 624.2 | 0.765 |
| WbbL | 9 | 244 | 629.7 | 0.770 |
| **Total** | 18 | **526** | | |

---

## Instructions for use

### Running on your own data

1. **Prepare one CSV per image** of trimer localisations, in the format written
   by the upstream detection step (see *Input format* below).
2. **Place the CSVs** under `trimer-localisations/<STRAIN>/`.
3. **Register the files** in the `all_data` dictionary near the top of the
   `if __name__ == "__main__":` block of
   [`porin_patch_analysis.py`](porin_patch_analysis.py), adding one entry per
   strain with its `base_dir` and the list of `localisation_csvs`.
4. **Adjust parameters** in `image-analysis-hyperparams.yaml` if required.
5. Run `python porin_patch_analysis.py`.

### Input format

Each localisation CSV has three comment lines followed by a standard CSV table:

```
# Processed AFM image name: MG_C1_save-2025.07.21-15.21.04.930
# Image extent (xmin, xmax, ymin, ymax) [nm]:
#(0.0,500.00000000000006,0.0,500.00000000000006)
trimer_index,y_nm,x_nm
0,50.78125000000001,217.77343750000003
...
```

- Line 1 supplies the image name recorded in the output table.
- Line 3 supplies the image extent in nm, which defines the region of interest
  (ROI) used for the patch analysis.
- The table gives the centre of each detected OMP trimer in nanometres.

All three comment lines are required and are parsed positionally, so their order
must be preserved.

### Parameters

`image-analysis-hyperparams.yaml`:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `patch_params.max_edge_separation` | 27.5 nm | Maximum Delaunay edge length. Triangles with any edge longer than this are not treated as trimer-occupied membrane, so the space they cover becomes part of a trimer-free patch. |

The path to the parameter file is recorded alongside the loaded parameters at
run time so that every analysis can be traced back to the settings used.

---

## Data availability

The trimer localisation tables required to reproduce the analysis are included
in this repository. The raw AFM image stacks and the trimer-detection step that
produced these localisations are described in the manuscript and deposited at
[REPOSITORY / ACCESSION / DOI].

---

## License

This code is released under the MIT License — see [LICENSE](LICENSE).

## Citation

If you use this code, please cite the manuscript above. Citation metadata for
the software itself is in [CITATION.cff](CITATION.cff); the archived version of
this repository corresponding to the published article is available at
[ZENODO DOI].

## Contact
- Anna Scheeder, a.scheeder@ucl.ac.uk
- Bart Hoogenboom, b.hoogenboom@ucl.ac.uk


Issues and questions may also be raised at
<https://github.com/hoogenboom-lab/bacteria-imaging/issues>.
