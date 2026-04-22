# ANCOM-BC2 Heatmap & Trajectory Plotter

A reusable Python package to visualize **ANCOM-BC2 results** and **taxon trajectories** from microbiome time-series data.

Designed for flexible workflows in **QIIME2-based longitudinal microbiome analyses**.

---

## ✨ Features

### 📊 Heatmaps

- ANCOM-BC2 log fold change heatmaps
- Multiple timepoints support
- Flexible subgroup filtering (genotype, treatment, sex, etc.)
- Supports multiple taxonomic levels:
  - genus
  - family
  - phylum
- Optional cell annotations:
  - relative abundance
  - log fold change
  - none
- Dynamic row heights based on abundance
- Automatic significance filtering

---

<img width="1735" height="870" alt="grafik" src="https://github.com/user-attachments/assets/474d49aa-2b39-4143-aa7d-5f68f42b687e" />


---

### 📈 Trajectory Plots

- Time-series plots for individual taxa
- Two plotting modes:
  - `comparison` → direct comparison (e.g. male vs female)
  - `combo` → stratified subgroups (e.g. WT vs Apc within treatment)
- Mean or median trajectories
- Error bars:
  - IQR
  - bootstrap CI
- Optional individual sample trajectories
- ANCOM significance overlay


---

<img width="1166" height="766" alt="grafik" src="https://github.com/user-attachments/assets/acd08698-38da-488f-ad63-b444411231d2" />


---

### 🔁 Workflow Integration

- Unified interface combining:
  - heatmaps → overview
  - trajectory plots → detailed temporal dynamics
- Minimal notebook code required
- Designed for high-throughput subset analysis

---

## 📦 Installation

Install directly from GitHub:

```bash
pip install "git+https://github.com/chaseU2/ancombc2-heatmaps.git"
```

---

## 🚀 Quick Start

### 1. Import

```python
from ancombc2_heatmaps import (
    PlotWorkflow,
    HeatmapConfig,
    MetadataConfig,
    ComparisonConfig,
    PathConfig,
    SubsetSpec,
    TrajectoryConfig,
    TrajectoryMetadataConfig,
    TrajectoryPathConfig,
    TrajectoryPlotConfig,
)
```

---

### 2. Create Configs

#### Heatmap

```python
heatmap_config = HeatmapConfig(
    metadata=MetadataConfig(
        sample_col="sample_name",
        timepoint_col="time_point",
        comparison_col="sex",
        timepoints=[...],
        timepoint_map={...},
    ),
    comparison=ComparisonConfig(
        variable_name="sex",
        positive_class="male",
        negative_class="female",
    ),
    paths=PathConfig(
        base_table_dir="...",
        base_ancom_dir="...",
        metadata_path="...",
        output_dir="...",
    ),
)
```

---

#### Trajectory

```python
traj_config = TrajectoryConfig(
    metadata=TrajectoryMetadataConfig(
        sample_col="sample_name",
        timepoint_col="time_point",
        mouse_col="host_subject_id",
        comparison_col="sex",

        # Required for combo plots
        genotype_col="mice_model",
        treatment_col="description_of_treatment",

        timepoint_order=[...],
        timepoint_numeric_map={...},
        timepoint_label_map={...},
    ),
    paths=TrajectoryPathConfig(
        metadata_path="...",
        table_base="...",
        ancom_base="...",
    ),
)
```

---

### 3. Create Workflow

```python
workflow = PlotWorkflow(
    heatmap_config=heatmap_config,
    trajectory_config=traj_config,
)
```

---

## 📊 Example Usage

### Heatmap

```python
subset = SubsetSpec(
    label="WT_sham_phylum_ANCOM",
    title="WT | sham",
    filters={
        "mice_model": "WT",
        "description_of_treatment": "sham",
    }
)

workflow.plot_heatmap(subset, show=True)
```

---

### Trajectory Plot (recommended start)

```python
workflow.plot_trajectory(
    taxon_query="d__Bacteria;p__Verrucomicrobiota",
    plot_mode="comparison",
    comparison_levels=["female", "male"],
)
```

---

### Combo Plot (advanced)

```python
workflow.plot_trajectory(
    taxon_query="d__Bacteria;p__Verrucomicrobiota",
    plot_mode="combo",
    comparison_levels=["female", "male"],
    combo_groups=[("WT", "sham")],
)
```

---

### Combined Workflow

```python
workflow.plot_heatmap_and_trajectory(
    subset=subset,
    taxon_query="d__Bacteria;p__Verrucomicrobiota",
    trajectory_plot_mode="comparison",
    comparison_levels=["female", "male"],
)
```

---

## 📁 Expected Data Structure

### QIIME2 Tables

```
collapsed_tables/
└── phylum/
    └── by_sex/
        └── baseline1/
            table_baseline1_WT_sham_phylum_ANCOM.qza
```

---

### ANCOM-BC2 Exports

```
exported/
└── phylum/
    └── sex/
        └── baseline1/
            table_baseline1_WT_sham_phylum_ANCOM_sex_ANCOMBC2_exported/
                ├── lfc.jsonl
                ├── q.jsonl
                └── diff.jsonl
```

---

## ⚙️ Customization

### Heatmap cell text

```python
heatmap_config.cell_text_mode = "relative_abundance"
heatmap_config.cell_text_mode = "lfc"
heatmap_config.cell_text_mode = "none"
```

---

### Trajectory options

```python
traj_config.plot.estimator = "median"
traj_config.plot.error_style = "iqr"
traj_config.plot.show_individual_lines = False
```

---

## 🧠 Typical Workflow

1. Run ANCOM-BC2 in QIIME2  
2. Export results (`qiime tools export`)  
3. Generate heatmaps  
4. Identify taxa of interest  
5. Generate trajectory plots  

---

## ⚠️ Notes

- `subset.label` must match your file naming exactly
- metadata sample IDs must match QZA tables
- ANCOM export must contain:
  - `lfc.jsonl`
  - `q.jsonl`
  - `diff.jsonl`

---

## 🧩 Known Pitfall (BIOM TSV Import)

Exported BIOM tables contain a header like:

```
# Constructed from biom file
#OTU ID    sample1    sample2 ...
```

When loading with pandas, you **must skip the first line**:

```python
pd.read_csv(tsv_fp, sep="\t", skiprows=1)
```

Using `comment="#"` will break sample ID parsing.

---

## 📌 Use Cases

- Microbiome time-series analysis
- Treatment vs control comparisons
- Genotype-dependent effects
- Sex-specific microbiome dynamics
- Multi-factor experimental designs

---

## 👤 Author

Karl Balzer
