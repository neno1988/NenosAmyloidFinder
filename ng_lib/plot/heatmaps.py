import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.axes

np.random.seed(19680801)

def gradient_image(ax, direction=0.3, cmap_range=(0, 1), **kwargs):
    """
    Draw a gradient image based on a colormap.

    Parameters
    ----------
    ax : Axes
        The Axes to draw on.
    direction : float
        The direction of the gradient. This is a number in
        range 0 (=vertical) to 1 (=horizontal).
    cmap_range : float, float
        The fraction (cmin, cmax) of the colormap that should be
        used for the gradient, where the complete colormap is (0, 1).
    **kwargs
        Other parameters are passed on to `.Axes.imshow()`.
        In particular, *cmap*, *extent*, and *transform* may be useful.
    """
    phi = direction * np.pi / 2
    v = np.array([np.cos(phi), np.sin(phi)])
    X = np.array([[v @ [1, 0], v @ [1, 1]],
                  [v @ [0, 0], v @ [0, 1]]])
    a, b = cmap_range
    X = a + (b - a) / X.max() * X
    im = ax.imshow(X, interpolation='bicubic', clim=(0, 1),
                   aspect='auto', **kwargs)
    return im

def nice_heatmap_plot(vector, threshold = 0, name="This is the default name", cmap_range=(0.4, 0.8)):
    fig, ax = plt.subplots(figsize=(6, 6))

    # Background image (heatmap)
    vector = vector.reshape(1, -1)
    N = vector.shape[1]

    # Normalize the data
    if threshold==0:
        norm = Normalize(vmin=np.nanmin(vector), vmax=np.nanmax(vector))
        vector = norm(vector)
        # Scale the normalized data to the desired cmap_range
        heatmap = ax.imshow(vector, vmin = 0, vmax=1,cmap='jet', interpolation='nearest', extent=(0, N, 0, 10))

    else:
        # cut at treshold
        vector = np.nan_to_num(vector, nan=0)
        vector[vector>threshold] = 0
        vector[vector<threshold] = 1

        # Scale the normalized data to the desired cmap_range
        #heatmap = ax.imshow(vector, vmin=0, vmax=1, cmap='custom_cmap', interpolation='nearest', extent=(0, len(vector), 0, 10))
        heatmap = ax.imshow(vector, vmin = -1, vmax=1,cmap='bwr', interpolation='nearest', extent=(0, N, 0, 10))

    #heatmap = ax.imshow(vector, cmap='YlOrRd', interpolation=None, extent=(0, N, 0, 10))

    # Gradient image
    gradient_image(ax, direction=0, extent=(0, N, 5, 10), transform=ax.transData,
                cmap=plt.cm.binary, cmap_range=(0.1, 0.9), alpha=0.4)
    gradient_image(ax, direction=0, extent=(0, N, 0, 5), transform=ax.transData,
                cmap=plt.cm.binary, cmap_range=(0.9, 0.1), alpha=0.4)
    add_protein_data_string(ax, name, threshold, N)

    plt.axis([0, N, 0, 10])
    #plt.show()
    return fig, ax

def heatmaps_binary_non_binary(binary_data, non_binary_data, threshold = 0, name="This is the default name", force_cmap=None, bar_height=10, seg_bar_height=5, xticks=0):
    fig, ax = nice_heatmap_plot(non_binary_data, threshold=threshold, name=name)
    
    N = len(non_binary_data)

    binary_data = binary_data.reshape(1, -1)
    # Define a custom colormap with light gray and lighter blue
    if force_cmap is None:
        colors = ["lightgray", "black"]
    else:
        colors = force_cmap

    cmap = LinearSegmentedColormap.from_list("custom_binary", colors, N=2)
    ax.imshow(binary_data, vmin=0, vmax=1, cmap=cmap, interpolation='nearest', extent=(0, N, bar_height, bar_height +seg_bar_height), alpha=1)
    plt.axis([0, N, 0, bar_height+seg_bar_height])
    
     # Remove Y-ticks
    ax.set_yticks([])
    if xticks>0:
        xticks_vec = np.arange(0, ax.get_xlim()[1], xticks)
        ax.set_xticks(xticks_vec)

    return fig
def add_protein_data_string(ax: matplotlib.axes.Axes, name, threshold, protein_length):

    # Add text above the figure
    plt.suptitle(f"Name: {name}, Threshold: {threshold}, Total Length: {protein_length}", fontsize=12, y=0.35)

def set_default_figure_size(ax):
    # Set the axis limits
    ax.set_ylim(-30, 30 )

if __name__ == "__main__":
    # Example usage
    vector = np.random.random(100)
    vectorB = np.round(np.random.random(100))
    
    #nice_heatmap_plot(vector, cmap_range=(0.4, 0.8))
    heatmaps_binary_non_binary(vectorB, vector)
    
