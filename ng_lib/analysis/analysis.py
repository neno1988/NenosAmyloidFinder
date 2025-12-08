
import plotly as plt
import ng_lib as ng
import os
import json
import re
import numpy as np

def analyse_protein(output_folder, seq, description, name, threshold, SEG, xticks):
    fasta_protein = ng.utils.FastaSeq() # TODO: use official fasta struct here
    fasta_protein.description = description
    fasta_protein.seq = seq
    fasta_protein.name = name

    # Get SEG and ZipperDB data
    seg_data = ng.get_SEG_data(fasta_protein, seg=SEG)
    zipperDB_data = ng.get_zipperDB_data(fasta_protein)

    # Get Amylpred data if credentials are provided
    try:
        amylpred_data = ng.get_amylpred_data(fasta_protein.seq)
    except Exception as e:
        # TODO: reinstate warning
        # messagebox.showwarning("Warning", f"Failed to get Amylpred data: {str(e)}")
        print("Reinstate warning about Amylpred failure here")
        amylpred_data = None

    do_matplotlib_plots = False
    if do_matplotlib_plots:
        if SEG==0:
            # TODO: only return fig  in nice_heatmap_plot?
            white_cmap = ["white", "white"]
            fig = ng.heatmaps_binary_non_binary(seg_data, zipperDB_data[0], threshold=threshold, name=name, force_cmap=white_cmap, seg_bar_height=0, xticks=xticks)
            # Save the figure as a JPEG file
            file_name = fasta_protein.name+'_treshold_'+str(threshold)+'.jpg'
            fig.savefig(os.path.join(output_folder, file_name), format='jpeg', dpi=300)
            fig.show() 

        else:
            fig = ng.heatmaps_binary_non_binary(seg_data, zipperDB_data[0], threshold=threshold, name=name, xticks=xticks)
            # Save the figure as a JPEG file
            file_name = fasta_protein.name+'_treshold_'+str(threshold)+'.jpg'
            fig.savefig(os.path.join(output_folder, file_name), format='jpeg', dpi=300)
            fig.show() 
    
    # Create interactive plot
    IMAGE_HEIGHT = 200
    IMAGE_WIDTH = 1800
    fig = ng.visualization.plotly_heatmaps.create_interactive_heatmaps(
        lcr_data=seg_data.reshape(1, -1),
        seg_treshold = SEG,
        zipperdb_data=zipperDB_data,
        amylpred_data=amylpred_data if amylpred_data is not None else np.zeros_like(zipperDB_data),
        zipperdb_threshold=threshold,
        name=name,
        xticks=xticks,
        image_height=IMAGE_HEIGHT, 
        image_width=IMAGE_WIDTH 
    )
    #fig.show()

    # Save the figure as HTML
    name = re.sub(r'[^\w_. -]', '_', name)
    name_and_parameters = f"{name}_SEG{SEG}_ZDB{threshold}"
    #output_folder = os.path.dirname(output_folder)
    html_file = os.path.join(output_folder, f"{name_and_parameters}_analysis.html")
    html_file = os.path.normpath(html_file)
    png_file = os.path.join(output_folder, f"{name_and_parameters}_image.png")
    json_file = os.path.join(output_folder, f"{name_and_parameters}_data.json")
    
    fig.write_html(html_file)
    

    fig.data= [fig.data[0]]
    fig.update_layout(title_text=None)


    fig.write_image(png_file, engine='kaleido', width = IMAGE_WIDTH, height = IMAGE_HEIGHT)


    # Generate output data file
    data_file = dict()
    data_file["fasta_seq"] = fasta_protein.seq
    data_file["fasta_description"] = fasta_protein.description
    data_file["fasta_name"] = fasta_protein.name
    data_file["ZipperDB"] = zipperDB_data.tolist()
    data_file["Amylpred"] = amylpred_data.tolist()
    data_file["SEG"] = seg_data.tolist()


    # Save the JSON
    with open(json_file, 'w') as f:
        json.dump(data_file, f, indent = 4)
    
    # If Amylpred data is available, show it
    plot_amylpred_only = False
    if plot_amylpred_only and amylpred_data is not None:
        plt.figure(figsize=(10, 2))
        plt.imshow(amylpred_data, cmap='viridis', aspect='auto')
        plt.colorbar(label='Amylpred Score')
        plt.yticks([])
        plt.xticks(ticks=np.arange(len(seq)), fontsize=8)
        plt.xlabel('Amino Acid Position')
        plt.title(f'Amylpred Analysis - {name}')
        plt.show()

