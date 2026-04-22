# ANCOM-BC2 Heatmap & Trajectory Plotter

A reusable Python package to visualize **ANCOM-BC2 results** and **taxon trajectories** from microbiome time-series data.

Designed for flexible workflows in QIIME2-based analyses.

---

## ✨ Features

### 📊 Heatmaps

* ANCOM-BC2 log fold change heatmaps
* Multiple timepoints
* Flexible grouping (genotype, treatment, sex, etc.)
* Optional cell annotations:

  * relative abundance
  * log fold change
  * none
* Dynamic row heights based on abundance
* Automatic significance filtering

### 📈 Trajectory Plots

* Time-series plots for individual taxa
* Supports any grouping variable (e.g. sex, treatment)
* Mean or median trajectories
* Error bars:

  * IQR
  * bootstrap CI
* Optional individual sample trajectories
* ANCOM significance overlay

### 🔁 Workflow Integration

* Unified interface to combine:

  * heatmaps → overview
  * trajectory plots → detailed view
* Minimal notebook code required

---

## 📦 Installation

Install directly from GitHub:

```bash
pip install "git+https://github.com/chaseU2/ancombc2-heatmaps.git"
```



---

<img width="1565" height="1826" alt="grafik" src="https://github.com/user-attachments/assets/91e99085-4bd3-4c23-8c1b-fe6f2a6cbc63" />

---

<img width="1166" height="766" alt="grafik" src="https://github.com/user-attachments/assets/76a6784d-dac5-4802-821e-9a926781b104" />


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

#### Trajectory

```python
traj_config = TrajectoryConfig(
    metadata=TrajectoryMetadataConfig(
        sample_col="sample_name",
        timepoint_col="time_point",
        mouse_col="host_subject_id",
        comparison_col="sex",
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
    label="WT_sham_genus_ANCOM",
    title="WT | sham",
    filters={
        "mice_model": "WT",
        "description_of_treatment": "sham",
    }
)

workflow.plot_heatmap(subset, show=True)
```

---

### Trajectory Plot

```python
workflow.plot_trajectory(
    taxon_query="g_Akkermansia",
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
    taxon_query="g_Akkermansia",
    trajectory_plot_mode="combo",
    comparison_levels=["female", "male"],
    combo_groups=[("WT", "sham")],
)
```

---

## 📁 Expected Data Structure

### QIIME2 Tables

```
{timepoint}/table_{timepoint}_{subset_label}.qza
```

### ANCOM-BC2 Exports

```
{timepoint}/table_{timepoint}_{subset_label}_<variable>_ANCOMB_exported/
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

### Trajectory options

```python
traj_config.plot.estimator = "median"
traj_config.plot.error_style = "ci"
traj_config.plot.show_individual_lines = False
```

---

## 🧠 Typical Workflow

1. Run ANCOM-BC2 in QIIME2
2. Export results
3. Generate heatmaps
4. Select taxa of interest
5. Generate trajectory plots

---

## ⚠️ Notes

* `subset.label` must match your file naming
* metadata sample IDs must match QZA tables
* ANCOM export must contain:

  * `lfc.jsonl`
  * `q.jsonl`
  * `diff.jsonl`

---

## 📌 Use Cases

* Microbiome time-series analysis
* Treatment vs control comparisons
* Genotype-dependent effects
* Sex-specific microbiome dynamics

---

## 👤 Author

Karl Balzer
