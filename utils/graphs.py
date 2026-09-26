import colorsys

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from .metrics import normalized_auc

cmap_hex = dict()
cmap_hex["gold"] = "#C99D0B"
cmap_hex["red"] = "#C62828"
cmap_hex["blue"] = "#20498D"
cmap_hex["dark_blue"] = "#112445"

blend_cmap = LinearSegmentedColormap.from_list(
    "blend", ["#C99D0B", "#C62828","#20498D", "#112445"])
def standardize_graph(
    ax,
    title=None,
    xlabel=None,
    ylabel=None,
    title_size=24,
    label_size=18,
    tick_size=16,
    pad=20,
    include_secondary=True,
):
    """
    Apply standard formatting to a Matplotlib Axes object,
    including any secondary (twinx/twiny) axes.
    Gridlines are only applied to the primary axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object to style.
    title, xlabel, ylabel : str, optional
        Titles and labels for the axes.
    title_size : int, default=16
        Font size for the title.
    label_size : int, default=12
        Font size for axis labels.
    tick_size : int, default=11
        Font size for tick labels.
    pad : int, default=15
        Padding for title from the plot.
    include_secondary : bool, default=True
        If True, also applies settings to any secondary twin axes.
    """

    def _apply(ax_, draw_grid=False):

        ax_.set_title(ax_.get_title(), fontsize=title_size, pad=pad)

        # --- Axis labels ---
        xlbl = xlabel or ax_.get_xlabel()
        ylbl = ylabel or ax_.get_ylabel()
        ax_.set_xlabel(xlbl, fontsize=label_size)
        ax_.set_ylabel(ylbl, fontsize=label_size)

        # --- Tick labels ---
        ax_.tick_params(axis="both", which="major", labelsize=tick_size)
        ax_.tick_params(axis="both", which="minor", labelsize=tick_size - 1)

        # --- Grid only on the first axis ---
        if draw_grid:
            ax_.grid(alpha=0.3)
        else:
            ax_.grid(False)


    # Apply to main axis
    _apply(ax, draw_grid=True)

    # Apply to secondary twin axes (if any)
    if include_secondary:
        for other_ax in ax.figure.axes:
            if other_ax is not ax and other_ax.bbox.bounds == ax.bbox.bounds:
                _apply(other_ax, draw_grid=False)

    plt.draw()
    return ax

def saturation_curve(sims, y, title):
    fig, ax = plt.subplots(figsize=(16/2,9/2))
    ax.plot(sims,y,c=cmap_hex['blue'])
    ax.grid(alpha=0.2)

    ax.set_xlabel("Tanimoto Similarity")
    ax.fill_between(sims, np.min(y), y, color=cmap_hex['blue'], alpha=0.1)  # shade under the curve
    ax.fill_between(sims, y, np.max(y), color=cmap_hex['red'], alpha=0.1)

    auc = round(normalized_auc([i for i in range(len(sims))], y),2)

    ax.set_ylabel("Number of Clusters")
    ax.set_title(f"B{title}\n AUC : {auc}",fontdict={"size":11})
    standardize_graph(ax)
    sns.despine()



def heatmap(sims, sorted=True,xticks=None, yticks = None,xlabel=None,ylabel= None,save=False,save_path=None,save_name = None, save_format="pdf",prov_ax = None,title="", title_font=12,cmap=None,tick_font=6):
    if xticks is None:
        xticks = []
    if yticks is None:
        yticks = []
    if cmap is not None:
        cmap = cmap.reversed()
    if sorted:
        # sort rows and columns by descending max similarity
        row_order = np.argsort(-sims.max(axis=1))  # indices of rows
        col_order = np.argsort(-sims.max(axis=0))  # indices of cols

        sims_sorted = sims[row_order][:, col_order]
        if len(xticks) > 0:
            xticks = xticks[col_order]
        if len(yticks) > 0:
            yticks = yticks[col_order]
        sims = sims_sorted
    # plt.figure(figsize=(6, 5))
    if prov_ax is  None:
        print("no ax")
        ax = sns.heatmap(
            sims, vmin=0, vmax=1, cmap=cmap,  xticklabels=xticks, yticklabels=yticks
        )
    else:
        sns.heatmap(
            sims, vmin=0,vmax=1, cmap=cmap,  xticklabels=xticks, yticklabels=yticks,ax=prov_ax
        )
        ax = prov_ax
    standardize_graph(ax)
    ax.set_title(title,fontsize=title_font)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="both", which="major", labelsize=tick_font)
    ax.grid(False)

    plt.tight_layout()
    plt.draw()
    if save:
        plt.savefig(
            save_path.joinpath(f"{save_name}.{save_format}"),
            dpi=300,
            bbox_inches="tight"
        )

def get_stratified_colors(n, n_tones=3, sat=0.85, light_range=(0.35, 0.75)):
    """
    Return n maximally distinct colors, stratified by tone (lightness).

    - n_tones: how many lightness levels to use (e.g. 3 → dark/medium/light)
    - sat: saturation in [0,1]
    - light_range: (min_L, max_L) for lightness in HLS space
    """
    n_hues = int(np.ceil(n / n_tones))
    hues = np.linspace(0, 1, n_hues, endpoint=False)
    lights = np.linspace(light_range[0], light_range[1], n_tones)

    colors = []
    for L in lights:          # stratify by tone
        for H in hues:        # spread hues around the wheel
            r, g, b = colorsys.hls_to_rgb(H, L, sat)
            colors.append(to_hex((r, g, b)))
    return colors[:n]


def get_indices_from_bit_array(fps_subset,length=576):
    # subset is a 2D array of shape (N, 576)
    # np.where returns the indices where the value is 1
    rows, cols = np.where(fps_subset == 1)

    # We want to group the column indices (the bit positions) by row
    results = []
    for i in range(len(fps_subset)):
        # Get all column indices where the row index matches i
        row_bits = cols[rows == i].tolist()
        results.append(tuple(row_bits))

    return results
