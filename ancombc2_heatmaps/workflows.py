from dataclasses import dataclass

from .plotter import ANCOMBC2HeatmapPlotter
from .trajectory_plotter import TaxonTrajectoryPlotter
from .boxplot_trajectory_plotter import TaxonBoxplotTrajectoryPlotter


@dataclass
class PlotWorkflow:
    heatmap_config: object
    trajectory_config: object

    def __post_init__(self):
        self.heatmap_plotter = ANCOMBC2HeatmapPlotter(self.heatmap_config)
        self.trajectory_plotter = TaxonTrajectoryPlotter(self.trajectory_config)
        self.boxplot_trajectory_plotter = TaxonBoxplotTrajectoryPlotter(self.trajectory_config)

        self.heatmap_meta = self.heatmap_plotter.load_metadata()
        self.trajectory_meta = self.trajectory_plotter.load_metadata()

    def plot_heatmap(self, subset, **kwargs):
        return self.heatmap_plotter.plot_subset(
            meta_df=self.heatmap_meta,
            subset=subset,
            **kwargs,
        )

    def plot_heatmaps(self, subsets, **kwargs):
        return self.heatmap_plotter.plot_all_subsets(
            subsets=subsets,
            **kwargs,
        )

    def plot_trajectory(
        self,
        taxon_query,
        subset,
        comparison_levels=None,
    ):
        return self.trajectory_plotter.plot_taxon(
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
        )

    def plot_boxplot_trajectory(
        self,
        taxon_query,
        subset,
        comparison_levels=None,
        show_trend=True,
        trend_order=2,
    ):
        return self.boxplot_trajectory_plotter.plot_taxon_boxplot(
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
            show_trend=show_trend,
            trend_order=trend_order,
        )

    def plot_heatmap_and_trajectory(
        self,
        subset,
        taxon_query,
        comparison_levels=None,
        heatmap_kwargs=None,
    ):
        if heatmap_kwargs is None:
            heatmap_kwargs = {}

        self.plot_heatmap(
            subset=subset,
            **heatmap_kwargs,
        )

        self.plot_trajectory(
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
        )

    def plot_heatmap_and_boxplot_trajectory(
        self,
        subset,
        taxon_query,
        comparison_levels=None,
        heatmap_kwargs=None,
        show_trend=True,
        trend_order=2,
    ):
        if heatmap_kwargs is None:
            heatmap_kwargs = {}

        self.plot_heatmap(
            subset=subset,
            **heatmap_kwargs,
        )

        self.plot_boxplot_trajectory(
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
            show_trend=show_trend,
            trend_order=trend_order,
        )