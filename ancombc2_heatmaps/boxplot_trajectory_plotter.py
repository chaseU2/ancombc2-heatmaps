from __future__ import annotations

from typing import Optional, List

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from .trajectory_plotter import (
    TaxonTrajectoryPlotter,
    TrajectoryConfig,
    sig_to_label,
    clean_float_formatter,
)


class TaxonBoxplotTrajectoryPlotter(TaxonTrajectoryPlotter):
    """
    Boxplot-based trajectory plotter.

    Uses the same config, metadata loading, taxon matching,
    table loading and ANCOM significance logic as TaxonTrajectoryPlotter.

    Instead of connected mean/median lines, this plot shows:
    - per timepoint and group: boxplots
    - all individual sample points
    - optional smooth trend line per group
    - optional ANCOM significance labels
    """

    def plot_single_boxplot(
        self,
        df: pd.DataFrame,
        title: str,
        comparison_levels: List[str],
        sig_map: Optional[dict] = None,
        ylim=None,
        show_trend: bool = True,
        trend_order: int = 2,
        jitter: float = 0.18,
    ):
        if df.empty:
            print(f"No data: {title}")
            return

        df = self.apply_baseline(df).copy()
        df["tp_plot"] = pd.to_numeric(df["tp_plot"], errors="coerce")
        df["abundance"] = pd.to_numeric(df["abundance"], errors="coerce")
        df = df.dropna(subset=["tp_plot", "abundance", "group"])

        if df.empty:
            print(f"No usable numeric data: {title}")
            return

        fig, ax = plt.subplots(figsize=self.config.plot.figsize)

        palette = dict(
            zip(
                comparison_levels,
                sns.color_palette(n_colors=len(comparison_levels))
            )
        )

        # -------------------------------------------------
        # Boxplots
        # -------------------------------------------------
        sns.boxplot(
            data=df,
            x="tp_plot",
            y="abundance",
            hue="group",
            hue_order=comparison_levels,
            palette=palette,
            dodge=True,
            showfliers=False,
            width=0.65,
            linewidth=1.2,
            ax=ax,
        )

        # -------------------------------------------------
        # Sample points
        # -------------------------------------------------
        sns.stripplot(
            data=df,
            x="tp_plot",
            y="abundance",
            hue="group",
            hue_order=comparison_levels,
            palette=palette,
            dodge=True,
            jitter=jitter,
            alpha=0.65,
            size=3.5,
            linewidth=0.3,
            edgecolor="black",
            ax=ax,
        )

        # -------------------------------------------------
        # Smooth approximate trend line per group
        # -------------------------------------------------
        if show_trend:
            grouped = (
                df.groupby(["tp_plot", "group"], as_index=False)["abundance"]
                .median()
            )

            tp_sorted = sorted(df["tp_plot"].unique())
            x_positions = {tp: i for i, tp in enumerate(tp_sorted)}

            for group in comparison_levels:
                g = grouped[grouped["group"] == group].copy()
                if g.empty:
                    continue

                g["x_pos"] = g["tp_plot"].map(x_positions)
                g = g.dropna(subset=["x_pos", "abundance"])

                if len(g) < 2:
                    continue

                x = g["x_pos"].to_numpy(dtype=float)
                y = g["abundance"].to_numpy(dtype=float)

                order = min(trend_order, len(g) - 1)

                try:
                    coeffs = np.polyfit(x, y, deg=order)
                    poly = np.poly1d(coeffs)

                    x_smooth = np.linspace(x.min(), x.max(), 200)
                    y_smooth = poly(x_smooth)

                    ax.plot(
                        x_smooth,
                        y_smooth,
                        color=palette[group],
                        linewidth=2.2,
                        alpha=0.9,
                        label=f"{group} trend",
                    )
                except Exception:
                    ax.plot(
                        x,
                        y,
                        color=palette[group],
                        linewidth=2.2,
                        alpha=0.9,
                        label=f"{group} trend",
                    )

        # -------------------------------------------------
        # Axis labels and title
        # -------------------------------------------------
        ax.set_title(f"{title}\n(boxplots + sample points + approximate trend)", pad=15)
        ax.set_ylabel(self.config.plot.y_label)
        ax.set_xlabel("")
        ax.yaxis.set_major_formatter(FuncFormatter(clean_float_formatter))

        tp_vals = sorted(df["tp_plot"].dropna().unique())
        tp_labels = [
            self.config.metadata.timepoint_label_map.get(tp, str(tp))
            for tp in tp_vals
        ]

        ax.set_xticks(range(len(tp_vals)))
        ax.set_xticklabels(tp_labels, rotation=35, ha="right")

        if self.config.plot.y_lim == "auto_fix" and ylim is not None:
            ax.set_ylim(ylim)
        elif isinstance(self.config.plot.y_lim, tuple):
            ax.set_ylim(self.config.plot.y_lim)

        # -------------------------------------------------
        # Significance labels
        # -------------------------------------------------
        if self.config.plot.show_significance and sig_map is not None:
            y0, y1 = ax.get_ylim()
            y_range = y1 - y0
            offset = y_range * 0.015

            for i, tp_plot in enumerate(tp_vals):
                g = df[df["tp_plot"] == tp_plot]
                if g.empty:
                    continue

                visible_upper = g["abundance"].quantile(0.75)
                original_tps = sorted(
                    pd.to_numeric(g["tp"], errors="coerce").dropna().unique()
                )

                if tp_plot == 0 and self.config.plot.merge_baselines:
                    is_sig = any(sig_map.get(tp_num, False) for tp_num in [-7, -4, -1])
                else:
                    is_sig = any(sig_map.get(int(tp_num), False) for tp_num in original_tps)

                label = sig_to_label(is_sig)
                y = min(visible_upper + offset, y1 - y_range * 0.03)

                ax.text(
                    i,
                    y,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    color="black",
                )

        # -------------------------------------------------
        # Clean legend duplicates
        # -------------------------------------------------
        handles, labels = ax.get_legend_handles_labels()

        clean_handles = []
        clean_labels = []
        seen = set()

        for h, l in zip(handles, labels):
            if l not in seen and not l.endswith("trend"):
                clean_handles.append(h)
                clean_labels.append(l)
                seen.add(l)

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

    def plot_taxon_boxplot(
        self,
        taxon_query,
        plot_mode="full",
        comparison_levels=None,
        partial_groups=None,
        combo_groups=None,
        show_trend: bool = True,
        trend_order: int = 2,
    ):
        meta = self.load_metadata()

        if comparison_levels is None:
            comparison_levels = sorted(
                meta[self.config.metadata.comparison_col].dropna().unique()
            )

        if partial_groups is None:
            partial_groups = []

        if combo_groups is None:
            combo_groups = []

        jobs = []

        if plot_mode == "full":
            df = self.build_df(
                meta,
                taxon_query,
                plot_mode="full",
                comparison_levels=comparison_levels,
            )

            sig_map = self.build_ancom_significance_map(
                taxon_query,
                plot_mode="full",
                variable_name=self.config.metadata.comparison_col,
            )

            jobs.append((df, f"{taxon_query} — all samples", sig_map))

        elif plot_mode == "partial":
            for group in partial_groups:
                df = self.build_df(
                    meta,
                    taxon_query,
                    plot_mode="partial",
                    group=group,
                    comparison_levels=comparison_levels,
                )

                sig_map = self.build_ancom_significance_map(
                    taxon_query,
                    plot_mode="partial",
                    variable_name=self.config.metadata.comparison_col,
                    group=group,
                )

                jobs.append((df, f"{taxon_query} — {group}", sig_map))

        elif plot_mode == "combo":
            for group in combo_groups:
                df = self.build_df(
                    meta,
                    taxon_query,
                    plot_mode="combo",
                    group=group,
                    comparison_levels=comparison_levels,
                )

                sig_map = self.build_ancom_significance_map(
                    taxon_query,
                    plot_mode="combo",
                    variable_name=self.config.metadata.comparison_col,
                    group=group,
                )

                jobs.append((df, f"{taxon_query} — {group[0]} | {group[1]}", sig_map))

        else:
            raise ValueError("plot_mode must be 'full', 'partial', or 'combo'")

        ylim = (
            self.compute_global_ylim([j[0] for j in jobs])
            if self.config.plot.y_lim == "auto_fix"
            else None
        )

        for df, title, sig_map in jobs:
            self.plot_single_boxplot(
                df=df,
                title=title,
                comparison_levels=comparison_levels,
                sig_map=sig_map,
                ylim=ylim,
                show_trend=show_trend,
                trend_order=trend_order,
            )