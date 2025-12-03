import plotly.graph_objects as go
import numpy as np


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

def create_interactive_heatmaps(lcr_data, seg_treshold, zipperdb_data, amylpred_data, zipperdb_threshold, name, xticks=0, legend_position="top", name_position="left"):
    """
    Create interactive heatmaps using Plotly.
    ... (docstring) ...
    """

    data_length = len(lcr_data[0])
    combined_matrix = np.zeros((11, data_length))
    combined_matrix[0] = lcr_data[0]

    zipperdb_mask = zipperdb_data[0] <= zipperdb_threshold
    amylpred_mask = amylpred_data[0] > 4

    # Data mapping
    data_to_color = {
        2: 'lightgray',  # None
        3: 'blue',       # Only Amylpred
        4: 'yellow',     # Only ZipperDB
        5: 'green'       # Both
    }

    for i in range(1, 11):
        combined_matrix[i] = np.where(zipperdb_mask & amylpred_mask, 5, np.where(zipperdb_mask, 4, np.where(amylpred_mask, 3, 2)))
    fig = go.Figure(data=go.Heatmap(
        z=combined_matrix,
        zmin=0,
        zmax=5,
        colorscale=[
            [0, 'white'],       # No LCR
            [0.2, 'darkgray'],  # LCR
            [0.4, 'lightgray'], # None
            [0.6, 'blue'],      # Only Amylpred
            [0.8, 'yellow'],    # Only ZipperDB
            [1, 'green']        # Both
        ],
        showscale=False,
        hoverongaps=False,
        name='Combined'
    ))

    # Legend using plotly legend.
    legend_items = [
        {'label': 'No LCR', 'color': 'white'},
        {'label': 'LCR', 'color': 'darkgray'},
        {'label': 'None', 'color': 'lightgray'},
        {'label': 'Only Amylpred', 'color': 'blue'},
        {'label': 'Only ZipperDB', 'color': 'yellow'},
        {'label': 'Both', 'color': 'green'}
    ]

    for item in legend_items:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(color=item['color']), name=item['label']))

    fig.update_layout(
        showlegend=True,
        legend_orientation="h" if legend_position == "top" or legend_position == "bottom" else "v",
        legend=dict(x=0.5 if legend_position == "top" or legend_position == "bottom" else 1.02, y=1.1 if legend_position == "top" else 0.5, xanchor="center", yanchor="bottom" if legend_position == "top" else "middle"),
        title_text="Protein Analysis",
        title_x=0.5,
        yaxis=dict(showticklabels=False, fixedrange=True),
        margin=dict(l=150 if name_position == "left" else 50, r=50, t=100, b=50), # margin for name
        height=200, # fixed height
    )

    # Name placement.
    name_x = -0.05 if name_position == "left" else 0.5
    name_y = 0.5
    fig.add_annotation(x=name_x, y=name_y, xref='paper', yref='paper', text=name, showarrow=False, font=dict(size=14), textangle=-90 if name_position == "left" else 0)

    parameter_annotations_x = 1.25
    parameter_annotations_y = 0.2
    annotations_text = f"SEG Threshold: {seg_treshold}<br>ZipperDB Threshold: {zipperdb_threshold}"
    fig.add_annotation(x=parameter_annotations_x, y=parameter_annotations_y, 
                       xref='paper', yref='paper', text=annotations_text, 
                       showarrow=False, 
                       font=dict(size=14))


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

    return fig