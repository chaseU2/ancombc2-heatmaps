# ANCOM-BC2 Heatmap Plotter

Reusable Python package for generating **ANCOM-BC2 log fold change heatmaps** across timepoints.

This tool was developed for microbiome time-series analyses (e.g. 16S data) and is designed to be **flexible, reusable, and easy to adapt to new datasets**.

---

## ✨ Features

* Plot ANCOM-BC2 log fold change heatmaps across multiple timepoints
* Supports **any comparison variable** (e.g. sex, treatment, genotype)
* Flexible **subset selection** (e.g. WT vs Apc, sham vs irradiated, combinations)
* Optional **cell annotations**:

  * mean relative abundance
  * log fold change
  * none
* **Dynamic row heights** based on abundance
* Automatic handling of:

  * taxonomy labels
  * significance filtering
  * effect direction (positive vs negative group)

---

## 📦 Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/chaseU2/ancombc2-heatmaps.git
cd ancombc2-heatmaps
pip install -e .
```

---

## 🚀 Quick Start

```python
from ancombc2_heatmaps import (
    ANCOMBC2HeatmapPlotter,
    HeatmapConfig,
    MetadataConfig,
    ComparisonConfig,
    PathConfig,
    SubsetSpec
)

config = HeatmapConfig(
    metadata=MetadataConfig(
        sample_col="sample_name",
        timepoint_col="time_point",
        comparison_col="sex",
        timepoints=["baseline1", "baseline2", "baseline3", "day1", "day3", "day7", "day14"]
    ),
    comparison=ComparisonConfig(
        variable_name="sex",
        positive_class="male",
        negative_class="female"
    ),
    paths=PathConfig(
        base_table_dir="/path/to/qza_tables",
        base_ancom_dir="/path/to/ancom_exports",
        metadata_path="/path/to/metadata.tsv",
        output_dir="/path/to/output"
    )
)

subset = SubsetSpec(
    label="all_genus_ANCOM",
    title="all samples",
    filters={}
)

plotter = ANCOMBC2HeatmapPlotter(config)
meta = plotter.load_metadata()

plotter.plot_subset(meta, subset)
```

---

## 📁 Expected Data Structure

### QIIME2 tables

```
{timepoint}/table_{timepoint}_{subset_label}.qza
```

### ANCOM-BC2 exports

```
{timepoint}/table_{timepoint}_{subset_label}_<variable>_ANCOMB_exported/
    ├── lfc.jsonl
    ├── q.jsonl
    └── diff.jsonl
```

Example:

```
day7/table_day7_WT_sham_genus_ANCOM.qza
day7/table_day7_WT_sham_genus_ANCOM_sex_ANCOMB_exported/
```

---

## 🔧 Configuration Overview

### MetadataConfig

Defines how metadata is interpreted:

* `sample_col` → sample IDs
* `timepoint_col` → timepoint column
* `comparison_col` → grouping variable
* `timepoints` → order of timepoints
* `timepoint_map` → optional mapping of names

---

### ComparisonConfig

Defines the ANCOM comparison:

* `variable_name` → name used in ANCOM (e.g. "sex")
* `positive_class` → shown in **red**
* `negative_class` → shown in **blue**

---

### PathConfig

Defines file structure:

* `base_table_dir` → QZA tables
* `base_ancom_dir` → ANCOM outputs
* `metadata_path` → metadata file
* `output_dir` → plot output

---

### SubsetSpec

Defines subsets of your data:

```python
SubsetSpec(
    label="WT_sham_genus_ANCOM",
    title="WT | sham",
    filters={
        "mice_model": "WT",
        "description_of_treatment": "sham"
    }
)
```

---

## 📊 Multiple Plots

```python
subsets = [
    SubsetSpec(label="WT_sham_genus_ANCOM", title="WT | sham", filters={...}),
    SubsetSpec(label="WT_irradiated_genus_ANCOM", title="WT | irradiated", filters={...}),
]

plotter.plot_all_subsets(subsets)
```

---

## 🎨 Customization

### Cell annotations

```python
config.cell_text_mode = "relative_abundance"  # default
config.cell_text_mode = "lfc"
config.cell_text_mode = "none"
```

---

### Highlight taxa

```python
config.style.highlight_taxa = ["Akkermansia", "Alistipes"]
```

---

### Significance filtering

```python
config.min_sig_cells_per_taxon = 2
```

---

## ⚠️ Important Notes

* `subset.label` **must match your file names**
* Metadata sample IDs must match QZA table sample IDs
* ANCOM output must contain:

  * `lfc.jsonl`
  * `q.jsonl`
  * `diff.jsonl`

---

## 🧪 Typical Workflow

1. Run ANCOM-BC2 in QIIME2
2. Export results
3. Define config + subsets
4. Generate heatmaps

---

## 📌 Use Cases

* Microbiome time series analysis
* Genotype vs treatment comparisons
* Sex-specific microbiome effects
* Longitudinal intervention studies

---

## 📄 License

Add your preferred license here (e.g. MIT).

---

## 👤 Author

Karl Balzer
