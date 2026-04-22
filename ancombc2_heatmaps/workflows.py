from dataclasses import dataclass
from .plotter import ANCOMBC2HeatmapPlotter
from .trajectory_plotter import TaxonTrajectoryPlotter


@dataclass
class PlotWorkflow:
    heatmap_config: object
    trajectory_config: object

    def __post_init__(self):
        self.heatmap_plotter = ANCOMBC2HeatmapPlotter(self.heatmap_config)
        self.trajectory_plotter = TaxonTrajectoryPlotter(self.trajectory_config)

        self.heatmap_meta = self.heatmap_plotter.load_metadata()
        self.trajectory_meta = self.trajectory_plotter.load_metadata()

    def plot_heatmap(self, subset, **kwargs):
        return self.heatmap_plotter.plot_subset(
            meta_df=self.heatmap_meta,
            subset=subset,
            **kwargs
        )

    def plot_heatmaps(self, subsets, **kwargs):
        return self.heatmap_plotter.plot_all_subsets(
            subsets=subsets,
            **kwargs
        )

    def plot_trajectory(self, taxon_query, plot_mode="full", comparison_levels=None,
                        partial_groups=None, combo_groups=None):
        return self.trajectory_plotter.plot_taxon(
            taxon_query=taxon_query,
            plot_mode=plot_mode,
            comparison_levels=comparison_levels,
            partial_groups=partial_groups,
            combo_groups=combo_groups,
        )


    def plot_heatmap_and_trajectory(
        self,
        subset,
        taxon_query,
        trajectory_plot_mode="full",
        comparison_levels=None,
        partial_groups=None,
        combo_groups=None,
        heatmap_kwargs=None,
    ):
        if heatmap_kwargs is None:
            heatmap_kwargs = {}
    
        self.plot_heatmap(subset=subset, **heatmap_kwargs)
    
        self.plot_trajectory(
            taxon_query=taxon_query,
            plot_mode=trajectory_plot_mode,
            comparison_levels=comparison_levels,
            partial_groups=partial_groups,
            combo_groups=combo_groups,
        )