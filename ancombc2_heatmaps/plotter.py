from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import qiime2
from biom import Table

warnings.filterwarnings("ignore")


# =========================================================
# CONFIG CLASSES
# =========================================================

@dataclass
class HeatmapStyleConfig:
    title_fontsize: int = 13
    axis_label_fontsize: int = 11
    xtick_fontsize: int = 8
    top_tick_fontsize: int = 8
    ytick_fontsize: int = 7
    celltext_fontsize: int = 8
    cbar_tick_fontsize: int = 8
    cbar_label_fontsize: int = 9

    base_cell_height: float = 0.25
    base_cell_width: float = 0.78
    height_scale: float = 6.0
    rowheight_ra_min: float = 0.00
    rowheight_ra_max: float = 0.05

    figure_min_height: float = 4.0
    figure_max_height: float = 100.0
    figure_min_width: float = 9.0
    figure_max_width: float = 14.0

    left_margin: float = 0.52
    right_margin: float = 0.95
    top_margin: float = 0.82
    bottom_margin: float = 0.10

    heatmap_cmap: str = "RdBu_r"
    missing_color: str = "white"
    edge_color: str = "lightgray"
    edge_linewidth: float = 0.35
    split_line_color: str = "black"
    split_linewidth: float = 1.2

    highlight_taxa: List[str] = field(default_factory=list)
    sns_style: str = "white"
    sns_context: str = "talk"


@dataclass
class MetadataConfig:
    sample_col: str
    timepoint_col: str
    comparison_col: Optional[str] = None

    timepoint_map: Dict[str, str] = field(default_factory=dict)
    timepoints: List[str] = field(default_factory=list)
    allowed_values: Dict[str, List[str]] = field(default_factory=dict)


@dataclass
class ComparisonConfig:
    variable_name: str
    positive_class: str
    negative_class: Optional[str] = None
    effect_column: Optional[str] = None
    invert_sign: Optional[bool] = None


@dataclass
class PathConfig:
    base_table_dir: str
    base_ancom_dir: str
    metadata_path: str
    output_dir: str

    table_template: str = "{timepoint}/table_{timepoint}_{subset_label}.qza"
    ancom_template: str = "{timepoint}/table_{timepoint}_{subset_label}_ANCOMB_exported"

    lfc_filename: str = "lfc.jsonl"
    q_filename: str = "q.jsonl"
    diff_filename: str = "diff.jsonl"


@dataclass
class PlotTextConfig:
    title_template: str = (
        "ANCOM-BC2 log fold change ({positive_label} vs {negative_label})\n"
        "{subset_title} | number in cell = {cell_text_description} | "
        "row height = globally scaled mean relative abundance across significant cells | "
        "only taxa with ≥{min_sig_cells} significant cells"
    )
    x_label: str = "Timepoint"
    y_label: str = "Taxon"
    colorbar_template: str = (
        "ANCOM-BC2 log fold change\n"
        "red = higher in {positive_label}\n"
        "blue = higher in {negative_label}"
    )
    top_axis_count_template: str = "{positive_label}={positive_n} | {negative_label}={negative_n}"


@dataclass
class TaxonomyConfig:
    normalizer: Optional[Callable[[str], str]] = None
    formatter: Optional[Callable[[str], str]] = None


@dataclass
class SubsetSpec:
    label: str
    title: str
    filters: Dict[str, Union[str, List[str]]] = field(default_factory=dict)


@dataclass
class HeatmapConfig:
    metadata: MetadataConfig
    comparison: ComparisonConfig
    paths: PathConfig
    text: PlotTextConfig = field(default_factory=PlotTextConfig)
    style: HeatmapStyleConfig = field(default_factory=HeatmapStyleConfig)
    taxonomy: TaxonomyConfig = field(default_factory=TaxonomyConfig)

    q_cutoff: float = 0.05
    min_sig_cells_per_taxon: int = 1
    split_after_timepoint: Optional[str] = None
    remove_empty_rows_after_masking: bool = True

    # "relative_abundance", "lfc", "none"
    cell_text_mode: str = "relative_abundance"
    lfc_decimals: int = 2


# =========================================================
# HELPERS
# =========================================================

def read_table_auto(filepath: str) -> pd.DataFrame:
    if filepath.endswith(".tsv") or filepath.endswith(".txt"):
        return pd.read_csv(filepath, sep="\t")
    if filepath.endswith(".csv"):
        return pd.read_csv(filepath)
    return pd.read_csv(filepath, sep="\t")


def read_export_tsv(tsv_fp: str) -> pd.DataFrame:
    """
    Robust reader for exported TSV files.
    Handles both:
    - old style: header directly in first row
    - new style: one extra leading line before real header

    Also normalizes '#OTU ID' -> 'feature'
    """
    df = pd.read_csv(tsv_fp, sep="\t")

    if "#OTU ID" not in df.columns and "feature" not in df.columns:
        df = pd.read_csv(tsv_fp, sep="\t", skiprows=1)

    df = df.rename(columns={"#OTU ID": "feature"})
    return df


def load_qza_table_as_df(qza_fp: str) -> pd.DataFrame:
    art = qiime2.Artifact.load(qza_fp)
    table = art.view(Table)
    df = table.to_dataframe(dense=True)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def clean_tax_piece(x) -> Optional[str]:
    if x is None:
        return None
    x = str(x).strip()
    return x if x else None


def is_empty_tax(x) -> bool:
    if x is None:
        return True
    stripped = re.sub(r"^[a-z]__+", "", str(x)).strip().lower()
    return stripped in {"", "uncultured", "unclassified", "unknown", "ambiguous_taxa"}


def default_normalize_taxon_label(raw_tax) -> str:
    if raw_tax is None:
        return "_; _; _"

    if isinstance(raw_tax, (list, tuple)):
        parts = [clean_tax_piece(x) for x in raw_tax]
    else:
        parts = [clean_tax_piece(x) for x in str(raw_tax).split(";")]

    parts = [x for x in parts if x is not None]

    tax_map = {}
    for p in parts:
        p = p.strip()
        if "__" in p:
            prefix = p.split("__", 1)[0].lower()
            if not is_empty_tax(p):
                tax_map[prefix] = p

    top_rank = None
    for rank in ["o", "c", "p", "k", "d"]:
        if rank in tax_map:
            top_rank = tax_map[rank]
            break

    family = tax_map.get("f")
    genus = tax_map.get("g")

    if top_rank is None:
        top_rank = "_"
    if family is None:
        family = "_"
    if genus is None:
        genus = "_"

    return f"{top_rank}; {family}; {genus}"


def default_taxon_formatter(label: str) -> str:
    return str(label)


def make_unique_labels(labels: Sequence[str]) -> List[str]:
    counts = {}
    out = []
    for lab in labels:
        if lab not in counts:
            counts[lab] = 1
            out.append(lab)
        else:
            counts[lab] += 1
            out.append(f"{lab} ({counts[lab]})")
    return out


def format_percent_value(x: float) -> str:
    if pd.isna(x):
        return ""
    pct = x * 100
    if pct == 0:
        return "0%"
    if pct >= 10:
        return f"{pct:.0f}%"
    if pct >= 1:
        return f"{pct:.2f}%"
    if pct >= 0.01:
        return f"{pct:.2f}%"
    return "<0.01%"


def format_lfc_value(x: float, decimals: int = 2) -> str:
    if pd.isna(x):
        return ""
    return f"{x:.{decimals}f}"


def is_highlight_taxon(label: str, highlight_list: Sequence[str]) -> bool:
    label_lower = label.lower()
    return any(h.lower() in label_lower for h in highlight_list)


def compute_row_heights_from_relative_abundance(
    row_sig_ra_mean: pd.Series,
    base_cell_h: float = 0.25,
    height_scale: float = 6.0,
    ra_min: float = 0.00,
    ra_max: float = 0.05,
) -> np.ndarray:
    row_sig_ra_mean = row_sig_ra_mean.fillna(0).clip(lower=0)

    if ra_max <= ra_min:
        raise ValueError("rowheight_ra_max must be greater than rowheight_ra_min.")

    scaled = ((row_sig_ra_mean - ra_min) / (ra_max - ra_min)).clip(0, 1)
    row_heights = base_cell_h * (1 + height_scale * scaled.to_numpy())
    return row_heights


# =========================================================
# MAIN PLOTTER CLASS
# =========================================================

class ANCOMBC2HeatmapPlotter:
    def __init__(self, config: HeatmapConfig):
        self.config = config

        sns.set(style=config.style.sns_style, context=config.style.sns_context)

        if self.config.taxonomy.normalizer is None:
            self.config.taxonomy.normalizer = default_normalize_taxon_label
        if self.config.taxonomy.formatter is None:
            self.config.taxonomy.formatter = default_taxon_formatter

        valid_cell_text_modes = {"relative_abundance", "lfc", "none"}
        if self.config.cell_text_mode not in valid_cell_text_modes:
            raise ValueError(
                f"cell_text_mode must be one of {valid_cell_text_modes}, "
                f"got {self.config.cell_text_mode}"
            )

        ensure_dir(self.config.paths.output_dir)

    # -----------------------------------------------------
    # Metadata
    # -----------------------------------------------------
    def load_metadata(self) -> pd.DataFrame:
        fp = self.config.paths.metadata_path
        if not os.path.exists(fp):
            raise FileNotFoundError(f"Metadata file not found:\n{fp}")

        meta = read_table_auto(fp)
        self._validate_metadata_columns(meta)
        meta = self._prepare_metadata(meta)
        return meta

    def _validate_metadata_columns(self, meta: pd.DataFrame) -> None:
        required = [
            self.config.metadata.sample_col,
            self.config.metadata.timepoint_col,
        ]
        if self.config.metadata.comparison_col is not None:
            required.append(self.config.metadata.comparison_col)

        missing = [c for c in required if c not in meta.columns]
        if missing:
            raise ValueError(f"Missing metadata columns: {missing}")

    def _prepare_metadata(self, meta: pd.DataFrame) -> pd.DataFrame:
        meta = meta.copy()

        sample_col = self.config.metadata.sample_col
        time_col = self.config.metadata.timepoint_col

        meta[sample_col] = meta[sample_col].astype(str)
        meta[time_col] = meta[time_col].astype(str)

        if self.config.metadata.timepoint_map:
            meta[time_col] = meta[time_col].map(
                lambda x: self.config.metadata.timepoint_map.get(x, x)
            )

        for col in meta.columns:
            if meta[col].dtype == object:
                meta[col] = meta[col].astype(str).str.strip()

        for col, allowed in self.config.metadata.allowed_values.items():
            if col not in meta.columns:
                raise ValueError(f"allowed_values references unknown metadata column: {col}")
            meta = meta[meta[col].isin(allowed)].copy()

        if self.config.metadata.timepoints:
            meta = meta[meta[time_col].isin(self.config.metadata.timepoints)].copy()

        return meta

    def filter_metadata(self, meta: pd.DataFrame, subset: SubsetSpec) -> pd.DataFrame:
        out = meta.copy()

        for col, val in subset.filters.items():
            if col not in out.columns:
                raise ValueError(f"Subset filter references unknown metadata column: {col}")

            if isinstance(val, (list, tuple, set)):
                out = out[out[col].isin(list(val))].copy()
            else:
                out = out[out[col] == val].copy()

        return out

    # -----------------------------------------------------
    # Paths
    # -----------------------------------------------------
    def get_table_qza_path(self, timepoint: str, subset: SubsetSpec) -> str:
        rel = self.config.paths.table_template.format(
            timepoint=timepoint,
            subset_label=subset.label,
        )
        return os.path.join(self.config.paths.base_table_dir, rel)

    def get_ancom_export_dir(self, timepoint: str, subset: SubsetSpec) -> str:
        rel = self.config.paths.ancom_template.format(
            timepoint=timepoint,
            subset_label=subset.label,
        )
        return os.path.join(self.config.paths.base_ancom_dir, rel)

    # -----------------------------------------------------
    # Effect column handling
    # -----------------------------------------------------
    def detect_effect_column(self, df: pd.DataFrame) -> str:
        cfg = self.config.comparison

        if cfg.effect_column is not None:
            if cfg.effect_column not in df.columns:
                raise ValueError(
                    f"Configured effect_column '{cfg.effect_column}' not found. "
                    f"Available columns: {list(df.columns)}"
                )
            return cfg.effect_column

        prefix = f"{cfg.variable_name}::"
        candidates = [c for c in df.columns if c.startswith(prefix)]

        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            raise ValueError(
                f"Multiple candidate effect columns found for prefix '{prefix}': {candidates}. "
                f"Please specify comparison.effect_column explicitly."
            )

        raise ValueError(
            f"No effect column found for prefix '{prefix}'. "
            f"Available columns: {list(df.columns)}"
        )

    def infer_negative_class_from_effect_column(self, effect_col: str) -> str:
        cfg = self.config.comparison

        if cfg.negative_class is not None:
            return cfg.negative_class

        if "::" in effect_col:
            level = effect_col.split("::", 1)[1]
            if level != cfg.positive_class:
                return level

        return "reference"

    def determine_sign_inversion(self, effect_col: str) -> bool:
        cfg = self.config.comparison

        if cfg.invert_sign is not None:
            return cfg.invert_sign

        if "::" in effect_col:
            level = effect_col.split("::", 1)[1]
            return level != cfg.positive_class

        return False

    def convert_lfc_to_positive_class(self, raw_lfc: float, effect_col: str) -> float:
        invert = self.determine_sign_inversion(effect_col)
        return -raw_lfc if invert else raw_lfc

    # -----------------------------------------------------
    # Data builders
    # -----------------------------------------------------
    def build_mean_relative_abundance(
        self,
        meta_df: pd.DataFrame,
        subset: SubsetSpec,
    ) -> Tuple[pd.DataFrame, Dict[str, Dict[str, int]]]:
        cfg_meta = self.config.metadata
        comparison_cfg = self.config.comparison

        sample_col = cfg_meta.sample_col
        time_col = cfg_meta.timepoint_col
        comparison_col = cfg_meta.comparison_col

        meta = self.filter_metadata(meta_df, subset)

        mean_dict: Dict[str, pd.Series] = {}
        group_counts: Dict[str, Dict[str, int]] = {}

        for tp in cfg_meta.timepoints:
            qza_fp = self.get_table_qza_path(tp, subset)
            if not os.path.exists(qza_fp):
                continue

            tab = load_qza_table_as_df(qza_fp)

            tab.index = [self.config.taxonomy.normalizer(x) for x in tab.index]
            tab = tab.groupby(tab.index).sum()

            col_sums = tab.sum(axis=0)
            col_sums[col_sums == 0] = np.nan
            rel_tab = tab.div(col_sums, axis=1)

            tp_meta = meta[meta[time_col] == tp].copy()
            tp_samples = [s for s in tp_meta[sample_col] if s in rel_tab.columns]

            if len(tp_samples) == 0:
                continue

            mean_vals = rel_tab[tp_samples].mean(axis=1, skipna=True)
            mean_dict[tp] = mean_vals

            tp_meta_used = tp_meta[tp_meta[sample_col].isin(tp_samples)].copy()

            if comparison_col is not None and comparison_col in tp_meta_used.columns:
                pos_n = int((tp_meta_used[comparison_col] == comparison_cfg.positive_class).sum())

                neg_label = self.config.comparison.negative_class or "other"
                if neg_label == "reference":
                    neg_n = int((tp_meta_used[comparison_col] != comparison_cfg.positive_class).sum())
                else:
                    neg_n = int((tp_meta_used[comparison_col] == neg_label).sum())

                group_counts[tp] = {
                    comparison_cfg.positive_class: pos_n,
                    neg_label: neg_n,
                }
            else:
                group_counts[tp] = {}

        return pd.DataFrame(mean_dict), group_counts

    def read_ancom_for_subset(
        self,
        subset: SubsetSpec,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[str]]:
        lfc_matrix = {}
        sig_matrix = {}
        detected_effect_col = None

        for tp in self.config.metadata.timepoints:
            export_dir = self.get_ancom_export_dir(tp, subset)
            if not os.path.isdir(export_dir):
                continue

            lfc_fp = os.path.join(export_dir, self.config.paths.lfc_filename)
            q_fp = os.path.join(export_dir, self.config.paths.q_filename)
            diff_fp = os.path.join(export_dir, self.config.paths.diff_filename)

            if not (os.path.exists(lfc_fp) and os.path.exists(q_fp) and os.path.exists(diff_fp)):
                continue

            lfc = pd.read_json(lfc_fp, lines=True)
            q = pd.read_json(q_fp, lines=True)
            diff = pd.read_json(diff_fp, lines=True)

            tax_col_lfc = "taxon" if "taxon" in lfc.columns else lfc.columns[0]
            tax_col_q = "taxon" if "taxon" in q.columns else q.columns[0]
            tax_col_diff = "taxon" if "taxon" in diff.columns else diff.columns[0]

            effect_col = self.detect_effect_column(lfc)
            detected_effect_col = effect_col

            if effect_col not in q.columns or effect_col not in diff.columns:
                raise ValueError(
                    f"Effect column '{effect_col}' not found in q or diff for {export_dir}."
                )

            lfc_sub = lfc[[tax_col_lfc, effect_col]].rename(
                columns={tax_col_lfc: "taxon_raw", effect_col: "lfc_raw"}
            )
            q_sub = q[[tax_col_q, effect_col]].rename(
                columns={tax_col_q: "taxon_raw", effect_col: "q"}
            )
            diff_sub = diff[[tax_col_diff, effect_col]].rename(
                columns={tax_col_diff: "taxon_raw", effect_col: "diff"}
            )

            tmp = (
                lfc_sub
                .merge(q_sub, on="taxon_raw", how="inner")
                .merge(diff_sub, on="taxon_raw", how="inner")
                .copy()
            )

            tmp["taxon"] = tmp["taxon_raw"].apply(self.config.taxonomy.normalizer)
            tmp["lfc"] = tmp["lfc_raw"].apply(
                lambda x: self.convert_lfc_to_positive_class(x, effect_col)
            )
            tmp["significant"] = (tmp["q"] < self.config.q_cutoff) & (tmp["diff"] == True)

            tmp_lfc = tmp.groupby("taxon", as_index=True)["lfc"].mean()
            tmp_sig = tmp.groupby("taxon", as_index=True)["significant"].any()

            lfc_matrix[tp] = tmp_lfc
            sig_matrix[tp] = tmp_sig

        lfc_df = pd.DataFrame(lfc_matrix)
        sig_df = pd.DataFrame(sig_matrix)

        return lfc_df, sig_df, detected_effect_col

    # -----------------------------------------------------
    # Plot helpers
    # -----------------------------------------------------
    def _get_positive_negative_labels(self, effect_col: Optional[str]) -> Tuple[str, str]:
        pos = self.config.comparison.positive_class
        neg = self.config.comparison.negative_class

        if neg is None:
            if effect_col is not None:
                neg = self.infer_negative_class_from_effect_column(effect_col)
            else:
                neg = "reference"

        return pos, neg

    def _get_cell_text_description(self) -> str:
        mode = self.config.cell_text_mode
        if mode == "relative_abundance":
            return "mean relative abundance"
        if mode == "lfc":
            return "log fold change"
        return "no annotation"

    def _get_cell_text_value(
        self,
        row_name: str,
        col_name: str,
        mean_rel_df: pd.DataFrame,
        plot_df: pd.DataFrame,
    ) -> str:
        mode = self.config.cell_text_mode

        if mode == "relative_abundance":
            return format_percent_value(mean_rel_df.loc[row_name, col_name])

        if mode == "lfc":
            return format_lfc_value(plot_df.loc[row_name, col_name], self.config.lfc_decimals)

        if mode == "none":
            return ""

        raise ValueError(f"Unknown cell_text_mode: {mode}")

    # -----------------------------------------------------
    # Plotting
    # -----------------------------------------------------
    def plot_subset(
        self,
        meta_df: pd.DataFrame,
        subset: SubsetSpec,
        save_png: bool = True,
        save_pdf: bool = True,
        show: bool = True,
    ) -> None:
        mean_rel_df, group_counts = self.build_mean_relative_abundance(meta_df, subset)
        if mean_rel_df.empty:
            print(f"[{subset.label}] No mean relative abundance data available.")
            return

        lfc_df, sig_df, effect_col = self.read_ancom_for_subset(subset)
        if lfc_df.empty:
            print(f"[{subset.label}] No LFC data found.")
            return

        positive_label, negative_label = self._get_positive_negative_labels(effect_col)

        present_tps = [tp for tp in self.config.metadata.timepoints if tp in lfc_df.columns]
        lfc_df = lfc_df[present_tps]
        sig_df = sig_df[present_tps]

        sig_df = sig_df.astype("boolean").fillna(False).astype(bool)

        sig_counts = sig_df.sum(axis=1)
        keep_taxa = sig_counts >= self.config.min_sig_cells_per_taxon

        lfc_df = lfc_df.loc[keep_taxa].copy()
        sig_df = sig_df.loc[keep_taxa].copy()

        if lfc_df.empty:
            print(f"[{subset.label}] No taxa remaining after significance filter.")
            return

        common_taxa = lfc_df.index.intersection(sig_df.index).intersection(mean_rel_df.index)

        lfc_df = lfc_df.loc[common_taxa].copy()
        sig_df = sig_df.loc[common_taxa].copy()
        mean_rel_df = mean_rel_df.loc[common_taxa].copy()

        plot_df = lfc_df.mask(~sig_df)

        if self.config.remove_empty_rows_after_masking:
            plot_df = plot_df.loc[plot_df.notna().any(axis=1)]

        mean_rel_df = mean_rel_df.reindex(plot_df.index)
        sig_df = sig_df.reindex(plot_df.index).fillna(False)
        lfc_df = lfc_df.reindex(plot_df.index)

        if plot_df.empty:
            print(f"[{subset.label}] No plottable values remain.")
            return

        row_order = pd.DataFrame({
            "n_sig": sig_df.sum(axis=1),
            "mean_abs_lfc": lfc_df.abs().mean(axis=1),
        }).sort_values(["n_sig", "mean_abs_lfc"], ascending=[False, False]).index

        plot_df = plot_df.loc[row_order]
        mean_rel_df = mean_rel_df.loc[row_order]
        sig_df = sig_df.loc[row_order]
        lfc_df = lfc_df.loc[row_order]

        display_labels = [
            self.config.taxonomy.formatter(x) for x in plot_df.index
        ]
        display_labels = make_unique_labels(display_labels)

        sig_ra_df = mean_rel_df.where(sig_df)
        row_sig_ra_mean = sig_ra_df.mean(axis=1, skipna=True).fillna(0)

        row_heights = compute_row_heights_from_relative_abundance(
            row_sig_ra_mean=row_sig_ra_mean,
            base_cell_h=self.config.style.base_cell_height,
            height_scale=self.config.style.height_scale,
            ra_min=self.config.style.rowheight_ra_min,
            ra_max=self.config.style.rowheight_ra_max,
        )

        cmap = mpl.colormaps[self.config.style.heatmap_cmap].copy()
        cmap.set_bad(self.config.style.missing_color)

        vals = plot_df.to_numpy(dtype=float)
        if vals.size == 0 or np.isnan(vals).all():
            print(f"[{subset.label}] No valid numeric values to plot.")
            return

        max_abs = np.nanmax(np.abs(vals))
        if not np.isfinite(max_abs) or max_abs == 0:
            max_abs = 1

        masked_vals = np.ma.masked_invalid(vals)

        n_rows, n_cols = plot_df.shape

        y_edges = np.concatenate([[0], np.cumsum(row_heights)])
        y_centers = (y_edges[:-1] + y_edges[1:]) / 2

        x_edges = np.arange(n_cols + 1)
        x_centers = np.arange(n_cols) + 0.5

        total_height = row_heights.sum()
        fig_h = total_height + 2.2
        fig_w = n_cols * self.config.style.base_cell_width + 7.5

        fig_h = min(max(fig_h, self.config.style.figure_min_height), self.config.style.figure_max_height)
        fig_w = min(max(fig_w, self.config.style.figure_min_width), self.config.style.figure_max_width)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        mesh = ax.pcolormesh(
            x_edges,
            y_edges,
            masked_vals,
            cmap=cmap,
            vmin=-max_abs,
            vmax=max_abs,
            shading="flat",
            edgecolors=self.config.style.edge_color,
            linewidth=self.config.style.edge_linewidth,
        )

        ax.set_ylim(y_edges[-1], 0)

        cbar_ax = fig.add_axes([0.96, 0.6, 0.015, 0.18])
        cbar_label = self.config.text.colorbar_template.format(
            positive_label=positive_label,
            negative_label=negative_label,
        )

        cbar = fig.colorbar(mesh, cax=cbar_ax, label=cbar_label)
        cbar.ax.tick_params(labelsize=self.config.style.cbar_tick_fontsize)
        cbar.ax.yaxis.label.set_size(self.config.style.cbar_label_fontsize)

        if self.config.cell_text_mode != "none":
            white_text_threshold = 0.72 * max_abs

            for i, r in enumerate(plot_df.index):
                for j, c in enumerate(plot_df.columns):
                    val = plot_df.loc[r, c]
                    if pd.notna(val):
                        txt = self._get_cell_text_value(
                            row_name=r,
                            col_name=c,
                            mean_rel_df=mean_rel_df,
                            plot_df=plot_df,
                        )
                        text_color = "white" if abs(val) >= white_text_threshold else "black"
                        ax.text(
                            x_centers[j],
                            y_centers[i],
                            txt,
                            ha="center",
                            va="center",
                            fontsize=self.config.style.celltext_fontsize,
                            color=text_color,
                        )

        if self.config.split_after_timepoint is not None:
            if self.config.split_after_timepoint in plot_df.columns:
                split_pos = plot_df.columns.get_loc(self.config.split_after_timepoint) + 1
                ax.axvline(
                    split_pos,
                    color=self.config.style.split_line_color,
                    linewidth=self.config.style.split_linewidth,
                )

        title = self.config.text.title_template.format(
            positive_label=positive_label,
            negative_label=negative_label,
            subset_title=subset.title,
            min_sig_cells=self.config.min_sig_cells_per_taxon,
            cell_text_description=self._get_cell_text_description(),
        )

        ax.set_title(title, fontsize=self.config.style.title_fontsize, pad=12)
        ax.set_xlabel(self.config.text.x_label, fontsize=self.config.style.axis_label_fontsize)
        ax.set_ylabel(self.config.text.y_label, fontsize=self.config.style.axis_label_fontsize)

        ax.set_xticks(x_centers)
        ax.set_xticklabels(
            plot_df.columns,
            rotation=45,
            ha="right",
            fontsize=self.config.style.xtick_fontsize,
        )

        top_ax = ax.secondary_xaxis("top")
        top_ax.set_xticks(x_centers)

        top_labels = []
        for tp in plot_df.columns:
            counts = group_counts.get(tp, {})
            pos_n = counts.get(positive_label, "?")
            neg_n = counts.get(negative_label, "?")
            top_labels.append(
                self.config.text.top_axis_count_template.format(
                    positive_label=positive_label,
                    negative_label=negative_label,
                    positive_n=pos_n,
                    negative_n=neg_n,
                )
            )

        top_ax.set_xticklabels(
            top_labels,
            rotation=45,
            ha="left",
            fontsize=self.config.style.top_tick_fontsize,
        )
        top_ax.tick_params(axis="x", pad=6)

        ax.set_yticks(y_centers)
        ax.set_yticklabels(display_labels, rotation=0, fontsize=self.config.style.ytick_fontsize)

        for tick_label in ax.get_yticklabels():
            if is_highlight_taxon(tick_label.get_text(), self.config.style.highlight_taxa):
                tick_label.set_fontweight("bold")

        plt.subplots_adjust(
            left=self.config.style.left_margin,
            right=self.config.style.right_margin,
            top=self.config.style.top_margin,
            bottom=self.config.style.bottom_margin,
        )

        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

        out_base = (
            f"ancombc2_lfc_heatmap_"
            f"{self.config.comparison.variable_name}_"
            f"{subset.label}_"
            f"min{self.config.min_sig_cells_per_taxon}sig"
        )

        if save_png:
            out_png = os.path.join(self.config.paths.output_dir, out_base + ".png")
            plt.savefig(out_png, dpi=300, bbox_inches="tight")
            print(f"Saved: {out_png}")

        if save_pdf:
            out_pdf = os.path.join(self.config.paths.output_dir, out_base + ".pdf")
            plt.savefig(out_pdf, bbox_inches="tight")
            print(f"Saved: {out_pdf}")

        if show:
            plt.show()
        else:
            plt.close(fig)

    def plot_all_subsets(
        self,
        subsets: Sequence[SubsetSpec],
        save_png: bool = True,
        save_pdf: bool = True,
        show: bool = True,
    ) -> None:
        meta = self.load_metadata()
        for subset in subsets:
            self.plot_subset(
                meta_df=meta,
                subset=subset,
                save_png=save_png,
                save_pdf=save_pdf,
                show=show,
            )