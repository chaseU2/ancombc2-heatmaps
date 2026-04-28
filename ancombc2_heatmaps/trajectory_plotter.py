from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import qiime2
from biom import Table
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

from .plotter import SubsetSpec


@dataclass
class TrajectoryMetadataConfig:
    sample_col: str
    timepoint_col: str
    mouse_col: str
    comparison_col: str
    genotype_col: Optional[str] = None
    treatment_col: Optional[str] = None
    timepoint_order: List[str] = field(default_factory=list)
    timepoint_numeric_map: Dict[str, int] = field(default_factory=dict)
    timepoint_label_map: Dict[int, str] = field(default_factory=dict)
    allowed_values: Dict[str, List[str]] = field(default_factory=dict)


@dataclass
class TrajectoryPathConfig:
    metadata_path: str
    table_base: str
    ancom_base: str

    table_template: str = "{timepoint}/table_{timepoint}_{subset_label}.qza"
    ancom_template: str = "{timepoint}/table_{timepoint}_{subset_label}_{variable_name}_ANCOMB_exported"


@dataclass
class TrajectoryPlotConfig:
    estimator: str = "mean"
    error_style: str = "iqr"
    show_individual_lines: bool = False
    merge_baselines: bool = False
    y_lim: Union[str, Tuple[float, float]] = "auto_fix"
    show_significance: bool = False
    q_cutoff: float = 0.05
    figsize: Tuple[float, float] = (12, 8)
    line_styles: Dict[str, Union[str, Tuple[int, int]]] = field(default_factory=dict)
    sns_style: str = "whitegrid"
    sns_context: str = "talk"
    y_label: str = "relative abundance"


@dataclass
class TrajectoryConfig:
    metadata: TrajectoryMetadataConfig
    paths: TrajectoryPathConfig
    plot: TrajectoryPlotConfig = field(default_factory=TrajectoryPlotConfig)


RANK_TO_PREFIX = {
    "d": "d__", "k": "k__", "p": "p__", "c": "c__",
    "o": "o__", "f": "f__", "g": "g__",
}


def clean_float_formatter(x, pos):
    return f"{x:.4f}".rstrip("0").rstrip(".")


def load_qza_table_as_df(qza_fp: str) -> pd.DataFrame:
    art = qiime2.Artifact.load(qza_fp)
    table = art.view(Table)
    df = table.to_dataframe(dense=True)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def clean_tax_piece(x):
    if x is None:
        return None
    x = str(x).strip()
    return x if x else None


def is_empty_tax(x):
    if x is None:
        return True
    stripped = re.sub(r"^[a-z]__+", "", str(x)).strip().lower()
    return stripped in {"", "uncultured", "unclassified", "unknown", "ambiguous_taxa"}


def normalize_taxon_label(raw_tax):
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

    family = tax_map.get("f", "_")
    genus = tax_map.get("g", "_")
    top_rank = top_rank or "_"

    return f"{top_rank}; {family}; {genus}"


def parse_label(label):
    p = [x.strip() for x in str(label).split(";")]
    while len(p) < 3:
        p.append("_")
    return {"top": p[0], "family": p[1], "genus": p[2]}


def normalize_query(q):
    q = q.strip()

    if ";" in q:
        return {"mode": "exact", "value": q}

    m = re.match(r"([a-zA-Z])_(.+)", q)
    if m:
        r, name = m.group(1).lower(), m.group(2)
        if r not in RANK_TO_PREFIX:
            raise ValueError(f"Unsupported rank prefix: {r}")
        return {"mode": "rank", "rank": r, "value": RANK_TO_PREFIX[r] + name}

    raise ValueError(
        "Invalid taxon_query. Use e.g. 'g_Akkermansia', "
        "'f_Akkermansiaceae', 'p_Verrucomicrobiota', or a full exact label."
    )


def match_taxa(index, query):
    spec = normalize_query(query)
    matches = []

    for lab in index:
        parsed = parse_label(lab)

        if spec["mode"] == "exact":
            if lab == spec["value"]:
                matches.append(lab)
        else:
            r = spec["rank"]
            if r == "f" and parsed["family"] == spec["value"]:
                matches.append(lab)
            elif r == "g" and parsed["genus"] == spec["value"]:
                matches.append(lab)
            elif r in ["o", "c", "p", "k", "d"] and parsed["top"] == spec["value"]:
                matches.append(lab)

    return matches


def detect_effect_col(df, variable_name):
    candidates = [c for c in df.columns if c.startswith(f"{variable_name}::")]
    if len(candidates) >= 1:
        return candidates[0]
    return None


def sig_to_label(is_sig):
    return "*" if is_sig else "ns"


class TaxonTrajectoryPlotter:
    def __init__(self, config: TrajectoryConfig):
        self.config = config
        sns.set_theme(style=config.plot.sns_style, context=config.plot.sns_context)

    def load_metadata(self):
        meta = pd.read_csv(self.config.paths.metadata_path, sep="\t", dtype=str)
        meta = meta.copy()

        tp_col = self.config.metadata.timepoint_col
        meta["tp"] = meta[tp_col].map(self.config.metadata.timepoint_numeric_map)

        for col in meta.columns:
            if meta[col].dtype == object:
                meta[col] = meta[col].astype(str).str.strip()

        for col, allowed in self.config.metadata.allowed_values.items():
            if col not in meta.columns:
                raise ValueError(f"allowed_values references unknown metadata column: {col}")
            meta = meta[meta[col].isin(allowed)].copy()

        return meta

    def filter_metadata(self, meta: pd.DataFrame, subset: Optional[SubsetSpec]) -> pd.DataFrame:
        out = meta.copy()

        if subset is None:
            return out

        for col, val in subset.filters.items():
            if col not in out.columns:
                raise ValueError(f"Subset filter references unknown metadata column: {col}")

            if isinstance(val, (list, tuple, set)):
                out = out[out[col].isin(list(val))].copy()
            else:
                out = out[out[col] == val].copy()

        return out

    def get_table_qza(self, tp: str, subset: SubsetSpec) -> str:
        rel = self.config.paths.table_template.format(
            timepoint=tp,
            subset_label=subset.label,
            variable_name=self.config.metadata.comparison_col,
        )
        return os.path.join(self.config.paths.table_base, rel)

    def get_ancom_export_dir(self, tp: str, subset: SubsetSpec, variable_name: str) -> Optional[str]:
        rel = self.config.paths.ancom_template.format(
            timepoint=tp,
            subset_label=subset.label,
            variable_name=variable_name,
        )
        p = os.path.join(self.config.paths.ancom_base, rel)
        return p if os.path.isdir(p) else None

    def build_df(
        self,
        meta: pd.DataFrame,
        taxon_query: str,
        subset: SubsetSpec,
        comparison_levels: Sequence[str],
    ) -> pd.DataFrame:
        rows = []
        meta_group = self.filter_metadata(meta, subset)

        sample_col = self.config.metadata.sample_col
        mouse_col = self.config.metadata.mouse_col
        comp_col = self.config.metadata.comparison_col

        for tp in self.config.metadata.timepoint_order:
            fp = self.get_table_qza(tp, subset)

            if not os.path.exists(fp):
                continue

            df = load_qza_table_as_df(fp)
            df.index = [normalize_taxon_label(x) for x in df.index]
            df = df.groupby(df.index).sum()

            col_sums = df.sum(axis=0)
            col_sums[col_sums == 0] = np.nan
            rel = df.div(col_sums, axis=1)

            matches = match_taxa(rel.index, taxon_query)
            if not matches:
                continue

            vals = rel.loc[matches].sum()

            tp_num = self.config.metadata.timepoint_numeric_map.get(tp)
            m = meta_group.copy()
            m = m[m["tp"] == tp_num]
            m = m[m[sample_col].isin(vals.index)]
            m = m[m[comp_col].isin(comparison_levels)]

            for _, r in m.iterrows():
                rows.append([
                    r[mouse_col],
                    r[comp_col],
                    r["tp"],
                    vals.get(r[sample_col], np.nan),
                ])

        return pd.DataFrame(rows, columns=["mouse", "group", "tp", "abundance"]).dropna()

    def build_ancom_significance_map(
        self,
        taxon_query: str,
        subset: SubsetSpec,
        variable_name: str,
    ) -> Dict[int, bool]:
        sig_map = {}

        for tp in self.config.metadata.timepoint_order:
            export_dir = self.get_ancom_export_dir(tp, subset=subset, variable_name=variable_name)

            if export_dir is None:
                continue

            q_fp = os.path.join(export_dir, "q.jsonl")
            diff_fp = os.path.join(export_dir, "diff.jsonl")

            if not (os.path.exists(q_fp) and os.path.exists(diff_fp)):
                continue

            q = pd.read_json(q_fp, lines=True)
            diff = pd.read_json(diff_fp, lines=True)

            effect_col = detect_effect_col(q, variable_name)

            if effect_col is None or effect_col not in diff.columns:
                continue

            tax_col_q = "taxon" if "taxon" in q.columns else q.columns[0]
            tax_col_diff = "taxon" if "taxon" in diff.columns else diff.columns[0]

            q_sub = q[[tax_col_q, effect_col]].rename(
                columns={tax_col_q: "taxon_raw", effect_col: "q"}
            )
            diff_sub = diff[[tax_col_diff, effect_col]].rename(
                columns={tax_col_diff: "taxon_raw", effect_col: "diff"}
            )

            tmp = q_sub.merge(diff_sub, on="taxon_raw", how="inner").copy()
            tmp["taxon"] = tmp["taxon_raw"].apply(normalize_taxon_label)
            tmp["significant"] = (tmp["q"] < self.config.plot.q_cutoff) & (tmp["diff"] == True)

            sig_by_taxon = tmp.groupby("taxon", as_index=True)["significant"].any()

            matches = match_taxa(sig_by_taxon.index, taxon_query)
            tp_num = self.config.metadata.timepoint_numeric_map.get(tp)

            if tp_num is None:
                continue

            sig_map[tp_num] = bool(sig_by_taxon.loc[matches].any()) if matches else False

        return sig_map

    def apply_baseline(self, df):
        df = df.copy()

        if self.config.plot.merge_baselines:
            df["tp_plot"] = df["tp"].replace({
                -7: 0,
                -4: 0,
                -1: 0,
                1: 1,
                3: 3,
                7: 7,
                14: 14,
            })
        else:
            df["tp_plot"] = df["tp"]

        return df

    def compute_global_ylim(self, dfs):
        upper_bounds = []

        for df in dfs:
            if df.empty:
                continue

            df_plot = self.apply_baseline(df)

            for _, g in df_plot.groupby(["tp_plot", "group"]):
                vals = pd.to_numeric(g["abundance"], errors="coerce").dropna().to_numpy()

                if len(vals) == 0:
                    continue

                center = np.median(vals) if self.config.plot.estimator == "median" else np.mean(vals)

                if self.config.plot.error_style == "iqr":
                    upper = np.percentile(vals, 75)
                else:
                    rng = np.random.default_rng(42)
                    boots = []
                    func = np.median if self.config.plot.estimator == "median" else np.mean

                    for _ in range(1000):
                        sample = rng.choice(vals, size=len(vals), replace=True)
                        boots.append(func(sample))

                    upper = np.percentile(boots, 97.5)

                upper_bounds.append(max(center, upper))

        if not upper_bounds:
            return (0, 1)

        ymax = max(upper_bounds)

        if not np.isfinite(ymax) or ymax <= 0:
            return (0, 1)

        return (0, ymax * 1.03)

    def _get_visible_upper(self, g):
        uppers = []

        for _, gs in g.groupby("group"):
            vals = pd.to_numeric(gs["abundance"], errors="coerce").dropna().to_numpy()

            if len(vals) == 0:
                continue

            if self.config.plot.error_style == "iqr":
                upper = np.percentile(vals, 75)
            else:
                rng = np.random.default_rng(42)
                func = np.median if self.config.plot.estimator == "median" else np.mean
                boots = []

                for _ in range(1000):
                    sample = rng.choice(vals, size=len(vals), replace=True)
                    boots.append(func(sample))

                upper = np.percentile(boots, 97.5)

            center = np.median(vals) if self.config.plot.estimator == "median" else np.mean(vals)
            uppers.append(max(center, upper))

        return max(uppers) if uppers else np.nan

    def plot_single(self, df, title, comparison_levels, sig_map=None, ylim=None):
        if df.empty:
            print(f"No data: {title}")
            return

        df = self.apply_baseline(df)

        fig, ax = plt.subplots(figsize=self.config.plot.figsize)

        dashes = self.config.plot.line_styles if self.config.plot.line_styles else True

        if self.config.plot.show_individual_lines:
            sns.lineplot(
                data=df,
                x="tp_plot",
                y="abundance",
                hue="group",
                style="group",
                dashes=dashes,
                units="mouse",
                estimator=None,
                alpha=0.35,
                linewidth=1,
                legend=False,
                ax=ax,
            )

        error_setting = ("ci", 95) if self.config.plot.error_style == "ci" else ("pi", 50)

        sns.lineplot(
            data=df,
            x="tp_plot",
            y="abundance",
            hue="group",
            style="group",
            dashes=dashes,
            estimator=self.config.plot.estimator,
            errorbar=error_setting,
            marker="o",
            linewidth=2.5,
            hue_order=comparison_levels,
            style_order=comparison_levels,
            ax=ax,
        )

        estimator_label = "median" if self.config.plot.estimator == "median" else "mean"
        err_label = "± 95% CI" if self.config.plot.error_style == "ci" else "± IQR"

        ax.set_title(f"{title}\n({estimator_label} {err_label})", pad=15)
        ax.set_ylabel(self.config.plot.y_label)
        ax.set_xlabel("")
        ax.yaxis.set_major_formatter(FuncFormatter(clean_float_formatter))

        tp_vals = sorted(df["tp_plot"].dropna().unique())
        tp_labels = [
            self.config.metadata.timepoint_label_map.get(tp, str(tp))
            for tp in tp_vals
        ]

        ax.set_xticks(tp_vals)
        ax.set_xticklabels(tp_labels, rotation=35, ha="right")

        if self.config.plot.y_lim == "auto_fix" and ylim is not None:
            ax.set_ylim(ylim)
        elif isinstance(self.config.plot.y_lim, tuple):
            ax.set_ylim(self.config.plot.y_lim)

        if self.config.plot.show_significance and sig_map is not None:
            y0, y1 = ax.get_ylim()
            y_range = y1 - y0
            offset = y_range * 0.012

            for tp_plot in sorted(df["tp_plot"].dropna().unique()):
                g = df[df["tp_plot"] == tp_plot]

                visible_upper = self._get_visible_upper(g)

                if pd.isna(visible_upper):
                    continue

                original_tps = sorted(pd.to_numeric(g["tp"], errors="coerce").dropna().unique())

                if tp_plot == 0 and self.config.plot.merge_baselines:
                    is_sig = any(sig_map.get(tp_num, False) for tp_num in [-7, -4, -1])
                else:
                    is_sig = any(sig_map.get(int(tp_num), False) for tp_num in original_tps)

                y = min(visible_upper + offset, y1 - y_range * 0.03)

                ax.text(
                    tp_plot,
                    y,
                    sig_to_label(is_sig),
                    ha="center",
                    va="bottom",
                    fontsize=11,
                )

        handles, labels = ax.get_legend_handles_labels()
        seen = set()
        clean_handles = []
        clean_labels = []

        for h, l in zip(handles, labels):
            if l not in seen:
                clean_handles.append(h)
                clean_labels.append(l)
                seen.add(l)

        if self.config.plot.show_significance:
            star_legend = Line2D(
                [0],
                [0],
                color="black",
                lw=0,
                label="* ANCOM significant\nns : not significant",
            )
            clean_handles.append(star_legend)
            clean_labels.append(star_legend.get_label())

        ax.legend(
            handles=clean_handles,
            labels=clean_labels,
            loc="upper left",
            fontsize=10,
            frameon=True,
            framealpha=0.9,
        )

        plt.tight_layout()
        plt.show()

    def plot_taxon(
        self,
        taxon_query: str,
        subset: SubsetSpec,
        comparison_levels: Optional[Sequence[str]] = None,
    ):
        meta = self.load_metadata()

        if comparison_levels is None:
            comparison_levels = sorted(
                meta[self.config.metadata.comparison_col].dropna().unique()
            )

        df = self.build_df(
            meta=meta,
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
        )

        sig_map = self.build_ancom_significance_map(
            taxon_query=taxon_query,
            subset=subset,
            variable_name=self.config.metadata.comparison_col,
        )

        ylim = self.compute_global_ylim([df]) if self.config.plot.y_lim == "auto_fix" else None

        self.plot_single(
            df=df,
            title=f"{taxon_query} — {subset.title}",
            comparison_levels=comparison_levels,
            sig_map=sig_map,
            ylim=ylim,
        )

    def list_available_queries(self, qza_fp: str) -> Dict[str, List[str]]:
        df = load_qza_table_as_df(qza_fp)
        taxa = [normalize_taxon_label(x) for x in df.index]

        families = set()
        phyla = set()
        genera = set()

        for t in taxa:
            parts = [p.strip() for p in str(t).split(";")]

            top = parts[0] if len(parts) > 0 else "_"
            family = parts[1] if len(parts) > 1 else "_"
            genus = parts[2] if len(parts) > 2 else "_"

            if family.startswith("f__"):
                name = re.sub(r"^f__", "", family)
                if name and name.lower() not in {
                    "uncultured",
                    "unclassified",
                    "unknown",
                    "ambiguous_taxa",
                    "_",
                }:
                    families.add(f"f_{name}")

            if genus.startswith("g__"):
                name = re.sub(r"^g__", "", genus)
                if name and name.lower() not in {
                    "uncultured",
                    "unclassified",
                    "unknown",
                    "ambiguous_taxa",
                    "_",
                }:
                    genera.add(f"g_{name}")

            if top.startswith("p__"):
                name = re.sub(r"^p__", "", top)
                if name and name.lower() not in {
                    "uncultured",
                    "unclassified",
                    "unknown",
                    "ambiguous_taxa",
                    "_",
                }:
                    phyla.add(f"p_{name}")

        return {
            "family_queries": sorted(families),
            "phylum_queries": sorted(phyla),
            "genus_queries": sorted(genera),
        }