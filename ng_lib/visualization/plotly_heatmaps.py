import plotly.graph_objects as go
import numpy as np

class HeatmapElement:
    data: np.ndarray
    number_of_rows: int
    min_value: float
    max_value: float
    color_segments: int
    legend: str
    colors: list[str]

    def __init__(self, data, number_of_rows, min_value, max_value, colors, legend, color_segments=2):
        self.data = data
        self.number_of_rows = number_of_rows
        self.min_value = min_value
        self.max_value = max_value
        self.colors = colors
        self.legend = legend
        self.color_segments = color_segments
    
    def get_gui_parameters(self):
        return {
            "number_of_rows": "integer",
            "colors": "List[str]",
            "color_segments": "integer"
        }

def get_xticks_interval(data_length):
        # Natural xtick interval
        if data_length > 0:
            if data_length <= 20:
                tick_interval = 1
            elif data_length <= 50:
                tick_interval = 5
            elif data_length <= 100:
                tick_interval = 10
            elif data_length <= 250:
                tick_interval = 25
            elif data_length <= 500:
                tick_interval = 50
            else:
                tick_interval = 100

            return tick_interval

def make_axes(fig, data_length, xticks):
    fig.update_xaxes(title_text="Amino Acid Position")

    # Dynamic xticks or fixed
    if xticks > 0:
        x_ticks_interval = xticks
    else:
        x_ticks_interval = get_xticks_interval(data_length)

    fig.update_layout(
            xaxis = dict(
                tickmode = 'linear',
                tick0 = 0,
                dtick = x_ticks_interval,
            ),
            font=dict(size=18, color="black"))


def make_annotations(fig, name, name_position, seg_treshold, zipperdb_threshold):
        # Keep heatmap domain fixed; leave blank space on the right for annotations.
    fig.update_xaxes(domain=[0, 0.7])

    # Name and annotation placement. TODO: make these available in a configuration file
    name_x = -0.05 if name_position == "left" else 0.5
    name_y = 0.5
    fig.add_annotation(x=name_x, y=name_y, xref='paper', yref='paper', text=name, showarrow=False, font=dict(size=14), textangle=-90 if name_position == "left" else 0)

    parameter_annotations_x = 0.85
    parameter_annotations_y = 0.2
    annotations_text = f"SEG Threshold: {seg_treshold}<br>ZipperDB Threshold: {zipperdb_threshold}"
    fig.add_annotation(x=parameter_annotations_x, y=parameter_annotations_y, 
                       xref='paper', yref='paper', text=annotations_text, 
                       showarrow=False, 
                       font=dict(size=14))
    return fig

def normalize_data(data, min=None, max=None ):
    if min is None:
        min = np.min(data)
    if max is None:
        max = np.max(data)
    normalized = (data - min) / (max - min)
    return normalized

def prepare_colormap_matrix(colormap_data: list[HeatmapElement], data_length, heatmap_height):

    # initializations
    colormap_matrix = np.zeros((heatmap_height, data_length))
    hover = np.zeros((heatmap_height, data_length))
    current_row = 0    
    colormap_offset = 0

    for cm in colormap_data:
        from_row = current_row
        to_row = current_row + cm.number_of_rows
        norm_data = normalize_data(cm.data, max=cm.max_value, min=cm.min_value)
        colormap_matrix[from_row:to_row] = colormap_offset + norm_data # shift colors to leave the first spots to previous
        hover[from_row:to_row] = cm.data
        current_row += cm.number_of_rows
        colormap_offset += cm.color_segments
    return colormap_matrix, hover

def create_interactive_heatmaps(heatmaps: list[HeatmapElement],
                                name, 
                                xticks=0, 
                                legend_position="top", 
                                name_position="left",
                                image_height=200, 
                                image_width=1800):
    """
    Create interactive heatmaps using Plotly.
    maps are normalized from 0 to 1 and then mapped to color segments.
    color segments ate specified for each heatmap elements. 
    in this funcion, all elements are combined into a single ta matrix, hover matrix, and color scale.
    """
  
    colormap_height = sum([heatmap.number_of_rows for heatmap in heatmaps])
    colormap_matrix, hover = prepare_colormap_matrix(
        heatmaps, 
        data_length=len(heatmaps[0].data),
        heatmap_height=colormap_height)
   
    all_colors = [c for heatmap in heatmaps for c in heatmap.colors]
    number_of_bars = sum([heatmap.color_segments for heatmap in heatmaps])

    # colorscale: generates the color vecctot and correct mapping for the colorbars
    colorscale = [[i / (len(all_colors) - 1), color] for i, color in enumerate(all_colors)]

    fig = go.Figure(data=go.Heatmap(
        z=colormap_matrix,
        zmin=0,
        zmax=number_of_bars - 1,
        colorscale=colorscale,
        showscale=False,
        hoverongaps=False,
        name='Combined',
        text=hover,
        hovertemplate = "%{text}<extra></extra>"
    ))

    # Legend using plotly legend, each line: dict containing {'label': 'XXX', 'color': 'YYY'}
    legend_items = [{'label': heatmap.legend, 'color': heatmap.colors[-1]} for heatmap in heatmaps if heatmap.legend]

    for item in legend_items:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color=item['color']), name=item['label']))

    fig.update_layout(
        showlegend=True,
        legend_orientation="h" if legend_position == "top" or legend_position == "bottom" else "v",
        legend=dict(x=0.3 if legend_position == "top" or legend_position == "bottom" else 1.02, y=1.1 if legend_position == "top" else 0.5, xanchor="center", yanchor="bottom" if legend_position == "top" else "middle"),
        title_text="Protein Analysis",
        title_x=0.3,
        yaxis=dict(showticklabels=False, fixedrange=True),
        # grow right margin for annotations without stretching plot content
        margin=dict(l=150 if name_position == "left" else 50, r=200, t=100, b=50),
        height=image_height, # fixed height
        width=image_width, # fixed width
        
    )
    

    return fig

def test_create_interactive_heatmaps():
    # Test data
    zipperdb_threshold = -23
    lcr_data = np.random.randint(0, 2, size=(1, 100))
    zipperdb_data = np.random.uniform(-30, 10, size=(1, 100))
    amylpred_data = np.random.randint(0, 7, size=(1, 100))
    zipperdb_data[0] = np.where(zipperdb_data[0] <= zipperdb_threshold, 1, 0)
    heatmaps: list[HeatmapElement] = []
    heatmaps.append(HeatmapElement(lcr_data[0], number_of_rows=2, min_value=0, max_value=1, colors = ["White", "orange"], legend = "LCR"))
    heatmaps.append(HeatmapElement(amylpred_data[0], number_of_rows=7, min_value=0, max_value=6, colors = ["White", "green"], legend = "Amylpred"))
    heatmaps.append(HeatmapElement(zipperdb_data[0], number_of_rows=2, min_value=0, max_value=1, colors = ["White", "darkgray"], legend = "ZipperDB"))
    name = "Test Protein"
    name_position = "left"
    seg_treshold=25
    fig = create_interactive_heatmaps(
        heatmaps=heatmaps,
        name="Test Protein",
        xticks=10
    )
    make_annotations(fig, name, name_position, seg_treshold, zipperdb_threshold)
    fig.show()

if __name__ == "__main__":
    test_create_interactive_heatmaps()  
