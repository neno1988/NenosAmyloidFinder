import plotly.graph_objects as go
import numpy as np

def create_interactive_heatmaps(lcr_data, zipperdb_data, amylpred_data, threshold, name, xticks=0):
    """
    Create interactive heatmaps using Plotly.
    
    Args:
        lcr_data: numpy array of LCR data
        zipperdb_data: numpy array of ZipperDB data
        amylpred_data: numpy array of Amylpred data
        threshold: float for ZipperDB threshold
        name: str for protein name
        xticks: int for number of x-axis ticks (0 for auto)
    """
    # Create figure
    fig = go.Figure()

    # Create the combined matrix (11 rows x data length)
    data_length = len(lcr_data[0])
    combined_matrix = np.zeros((11, data_length))
    
    # First row: LCR data (0: no LCR, 1: LCR)
    combined_matrix[0] = lcr_data[0]
    
    # Create masks for the remaining rows
    zipperdb_mask = zipperdb_data[0] <= threshold
    amylpred_mask = amylpred_data[0] > 4
    
    # Fill the remaining 10 rows with the combined data
    for i in range(1, 11):
        combined_matrix[i] = np.where(
            zipperdb_mask & amylpred_mask, 5,  # Both (green)
            np.where(
                zipperdb_mask, 4,  # Only ZipperDB (yellow)
                np.where(
                    amylpred_mask, 3,  # Only Amylpred (blue)
                    2  # None (light gray)
                )
            )
        )

    # Add the combined heatmap
    fig.add_trace(
        go.Heatmap(
            z=combined_matrix,
            colorscale=[
                [0, 'white'],      # No LCR
                [0.2, 'darkgray'], # LCR
                [0.4, 'lightgray'],# None
                [0.6, 'blue'],     # Only Amylpred
                [0.8, 'yellow'],   # Only ZipperDB
                [1, 'green'],    # Both
                [1, 'green']
            ],
            showscale=True,
            colorbar=dict(
                title='Score',
                len=0.3,
                y=0.5,
                ticktext=['No LCR', 'LCR', 'None', 'Only Amylpred', 'Only ZipperDB', 'Both'],
                tickvals=[0.05, 0.15, 0.3, 0.5, 0.7, 0.9],
                tickmode='array'
            ),
            hoverongaps=False,
            name='Combined'
        )
    )

    # Update layout
    fig.update_layout(
        height=400,
        showlegend=False,
        title_text=f"Protein Analysis - {name}",
        title_x=0.5,
        margin=dict(t=50, b=50),
        yaxis=dict(
            showticklabels=False,
            fixedrange=True
        )
    )

    # Update axes
    fig.update_xaxes(title_text="Amino Acid Position")

    # Set x-axis ticks if specified
    if xticks > 0:
        x_ticks = np.linspace(0, data_length-1, xticks, dtype=int)
        fig.update_xaxes(tickvals=x_ticks)

    return fig 