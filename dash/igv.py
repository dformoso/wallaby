import dash
import dash_auth
import dash_bio as dashbio
from dash import dash_table as dt
from dash.dependencies import Input, Output
from dash import dcc
from dash import html
import os, glob, csv, time
import pandas as pd

batches = [ 
    { 'value' : 'http://159.196.33.135:8080/hpv16_rnaseq/', 'label' : 'HPV16 RNAseq Dataset' },
    { 'value' : 'http://159.196.33.135:8080/hpv18_rnaseq/', 'label' : 'HPV18 RNAseq Dataset' } 
]

donors = []
for batch in batches:
    donor = pd.read_csv(batch['value'] + 'donor_and_recipient.csv').columns[0].replace(" ", "")
    donors.append({ 'value' : donor, 'label' : donor })
recipients = []
for batch in batches:
    recipient = pd.read_csv(batch['value'] + 'donor_and_recipient.csv').columns[1].replace(" ", "")
    recipients.append({ 'value' : recipient, 'label' : recipient })


tabs_styles = {
    'height': '44px'
}
tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold'
}

tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#119DFF',
    'color': 'white',
    'padding': '6px',
    'fontWeight': 'bold'
}

app = dash.Dash(__name__)

VALID_USERNAME_PASSWORD_PAIRS = {'dash': 'r0tt3nf1sh'}
auth = dash_auth.BasicAuth( app, VALID_USERNAME_PASSWORD_PAIRS)

app.layout = html.Div([
    dcc.Tabs(
        id='tab-select', value='tab-0', style = tabs_styles,
        children=[
            dcc.Tab(
                label = 'Data Selection', value = 'tab-0', 
                style = tab_style, selected_style = tab_selected_style,
                children = [
                    html.Br(),
                    html.Div([
                        html.Div(
                            id='dataset', 
                            children = 'Please select a dataset name:',
                            style={'width': '20%', 'display': 'inline-block'}),
                        html.Div(
                            id='donor', 
                            children = 'Donor Organism:',
                            style={'width': '20%', 'display': 'inline-block'}),
                        html.Div(
                            id='recipient', 
                            children = 'Recipient Organism:',
                            style={'width': '20%', 'display': 'inline-block'}),
                        html.Div(
                            id='srr', 
                            children = 'SRR Identifier:',
                            style={'width': '20%', 'display': 'inline-block'}),
                        html.Div(
                            id='id', 
                            children = 'Overlap Locus ID:',
                            style={'width': '20%', 'display': 'inline-block'})
                    ]),
                    html.Div([
                        html.Div([
                            dcc.Dropdown(
                                id = 'batch-select',
                                options = batches, value = batches[0]['value'], 
                                clearable = False
                            )
                        ], style={'width': '20%', 'display': 'inline-block'}),
                        
                        html.Div([
                            dcc.Dropdown(
                                id = 'donor-select',
                                options = donors, value = donors[0]['value'], 
                                clearable = False, disabled = True
                            )
                        ], style={'width': '20%', 'display': 'inline-block'}),
                        
                        html.Div([
                            dcc.Dropdown(
                                id = 'recipient-select',
                                options = recipients, value = recipients[0]['value'], 
                                clearable = False, disabled = True
                            )
                        ], style={'width': '20%', 'display': 'inline-block'}),
                        
                        html.Div([
                            dcc.Dropdown(id = 'srr-select')
                        ], style={'width': '20%', 'display': 'inline-block'}),

                        html.Div([
                            dcc.Dropdown(id = 'id-select', value = '1')
                        ], style={'width': '20%', 'display': 'inline-block'})
                    ]),
                ]   
            ),
            dcc.Tab(label = 'Metrics Reports', value = 'tab-1', 
            style = tab_style, selected_style = tab_selected_style),
            dcc.Tab(label = 'Reads at Locus - Mapped to Recipient', value = 'tab-2', 
            style = tab_style, selected_style = tab_selected_style),
            dcc.Tab(label = 'Reads Mapped to Donor', value = 'tab-3', 
            style = tab_style, selected_style = tab_selected_style),
            dcc.Tab(label = 'Locus-filtered Reads', value = 'tab-4', 
            style = tab_style, selected_style = tab_selected_style),
            ]
        ),
    html.Div(id='tabs-content')
])

################
# Batch change > Update the SRRs and default SRR
@app.callback(
    Output('srr-select', 'options'),
    Output('srr-select', 'value'),
    Output('donor-select', 'value'),
    Output('recipient-select', 'value'),
    Input('batch-select', 'value')
)
def update_options(batch):
    insertion_table = pd.read_csv(batch + 'putative_insertion_table.csv', 
    dtype={'srr': str, 'id': str })
    # Extract ids for SRR and build dictionary from it
    srrs = insertion_table['srr'].drop_duplicates()

    donor = pd.read_csv(batch + 'donor_and_recipient.csv').columns[0].replace(" ", "")
    recipient = pd.read_csv(batch + 'donor_and_recipient.csv').columns[1].replace(" ", "")

    return [{ 'value' : srr, 'label' : srr } for srr in srrs], srrs[0], donor, recipient

###################
# Batch or SRR change > Update the SRR IDs
@app.callback(
    Output('id-select', 'options'),
    Input('srr-select', 'value'),
    Input('batch-select', 'value')
)
def update_options(srr, batch):
    insertion_table = pd.read_csv(batch + 'putative_insertion_table.csv', dtype={'srr': str, 'id': str })
    # Extract ids for SRR and build dictionary from it
    insertion_table_srr = insertion_table[insertion_table['srr'] == srr]
    ids = insertion_table_srr['id']
    return [{ 'value' : id, 'label' : id } for id in ids]

###########################
# Return Tab Components
@app.callback(
    Output('tabs-content', 'children'),
    Input('tab-select', 'value'),
    Input('srr-select', 'value'),
    Input('id-select', 'value'),
    Input('donor-select', 'value'),
    Input('recipient-select', 'value'),
    Input('batch-select', 'value')
)
def render_content(tab, srr, id, donor, recipient, batch):
    # Load insertion table
    insertion_table = pd.read_csv(batch + 'putative_insertion_table.csv', dtype={'srr': str, 'id': str })
    insertion_table_renamed = insertion_table.rename(columns = {
        "srr": "SRR name", "id": "Overlap locus ID", "chr": "Chromosome", "start": "Start locus",
        "stop": "End locus", "num_crossings": "Number of crossings", 
        "unique_crossings": "List of unique crossings", "num_reads": "Number of reads", 
        "gene_name": "Gene name"})
    # Defaults
    padding_left = 150
    padding_right = 150

    if tab == 'tab-0':
        return html.Div([
            html.Br(),
            html.Div(
                id='table', 
                children = 'Table of putative insertions:',
                style={'width': '20%', 'display': 'inline-block'}),
            dt.DataTable(
                id = 'tbl', data = insertion_table_renamed.to_dict('records'),
                columns = [{'name': i, 'id': i} for i in insertion_table_renamed.columns],
                style_cell = {'textAlign': 'left'}
            )
        ])

    if tab == 'tab-1':
        return html.Div([
            dcc.Tabs(style = tabs_styles,
                children = [
                dcc.Tab(label='Paired-End .FASTQ MultiQC Report',
                        style = tab_style, selected_style = tab_selected_style,
                        children=[
                    html.Iframe(
                        src = batch + srr + '_multiqc_report.html',
                        style = {"height": "1067px", "width": "100%"},
                    )
                ]),
                dcc.Tab(label='Paired-End Donor .BAM MultiQC Report',
                        style = tab_style, selected_style = tab_selected_style,
                        children=[
                    html.Iframe(
                        src = batch + srr + '-to-' + donor + '_multiqc_report.html',
                        style = {"height": "1067px", "width": "100%"}
                    )
                ]),
                dcc.Tab(label='Paired-End Recipient .BAM MultiQC Report',
                        style = tab_style, selected_style = tab_selected_style, children=[
                    html.Iframe(
                        src = batch + srr + '-to-' + recipient + '_multiqc_report.html',
                        style = {"height": "1067px", "width": "100%"}
                    )
                ]),
            ])
        ])



    elif tab == 'tab-2':
        # Extract track names
        insertion_table_srr = insertion_table[(insertion_table['srr'] == str(srr)) & (insertion_table['id'] == str(id))]
        crossings = str(insertion_table_srr['unique_crossings'].values[0]).split('|')

        # Calculate locus
        srr_id = insertion_table_srr[insertion_table_srr['id'] == id]
        chromosome = str(srr_id['chr'].values[0])
        start = str(srr_id['start'].values[0] - padding_left)
        stop = str(srr_id['stop'].values[0] + padding_right)
        locus = str(chromosome + ':' + start + '-' + stop)

        recipient_tracks = []
        # Create to-recipient alignment tracks
        for crossing in crossings:
            recipient_tracks.append(
                {
                    # Docs: 
                    # https://github.com/igvteam/igv.js/wiki/Alignment-Track
                    # https://github.com/igvteam/igv.js/blob/c6773940de86cf2938f40feec86d1905a866deba/js/bam/bamTrack.js
                    'name': crossing,
                    'url': batch + srr + '-to-' + recipient + '_' + crossing + '_filtered.bam',
                    'indexURL': batch + srr + '-to-' + recipient + '_' + crossing + '_filtered.bam.bai', 
                    'format': 'bam', 'type': 'alignment', 'showSoftClips': 'true', 'alignmentRowHeight' : '14',
                    'viewAsPairs': 'true', 'showMismatches' : 'true', 'colorBy' : 'strand', 
                    'coverageColor' : 'rgb(210, 100, 102)'
                }
            )
        
        # Return html
        return html.Div([
            html.Div([
                dashbio.Igv(
                    id = 'hg38',
                    genome = 'hg38',
                    locus = locus,
                    tracks = recipient_tracks
                )
            ], style={'width': '100%', 'display': 'inline-block'})
        ])

    elif tab == 'tab-3':
        # Extract track names
        insertion_table_srr = insertion_table[(insertion_table['srr'] == str(srr)) & (insertion_table['id'] == str(id))]
        crossings = str(insertion_table_srr['unique_crossings'].values[0]).split('|')

        # Calculate locus
        srr_id = insertion_table_srr[insertion_table_srr['id'] == id]
        chromosome = str(srr_id['chr'].values[0])
        start = str(srr_id['start'].values[0] - padding_left)
        stop = str(srr_id['stop'].values[0] + padding_right)
        locus = str(chromosome + ':' + start + '-' + stop)

        donor_tracks = []
        # Create to-donor alignment tracks
        
        for crossing in crossings:
            donor_tracks.append(
                {
                    # Docs: 
                    # https://github.com/igvteam/igv.js/wiki/Alignment-Track
                    # https://github.com/igvteam/igv.js/blob/c6773940de86cf2938f40feec86d1905a866deba/js/bam/bamTrack.js
                    'name': crossing,
                    'url': batch + srr + '-to-' + donor + '_' + crossing + '_filtered.bam',
                    'indexURL': batch + srr + '-to-' + donor + '_' + crossing + '_filtered.bam.bai', 
                    'format': 'bam', 'type': 'alignment', 'showSoftClips': 'true', 'alignmentRowHeight' : '14',
                    'viewAsPairs': 'true', 'showMismatches' : 'true', 'colorBy' : 'strand', 
                    'coverageColor' : 'rgb(210, 100, 102)'
                }
            )
        
        donor_tracks.append(
            {
                'name': 'Annotations',
                'url': batch + donor + '.gff', 'isCompressed' : 'false',
                'displayMode': 'EXPANDED',
                'nameField': 'gene'
            }
        )
        # Return html
        return html.Div([
            html.Div([
                dashbio.Igv(
                    id = donor,
                    reference = {
                        'id': donor,
                        'name': donor,
                        'fastaURL': batch + donor + '.fa',
                        'indexURL': batch + donor + '.fa.fai',
                        'tracks': donor_tracks,
                    },
                )
            ], style={'width': '100%', 'display': 'inline-block'})
        ])

    elif tab == 'tab-4':
        # Extract track names
        insertion_table_srr = insertion_table[(insertion_table['srr'] == str(srr)) & (insertion_table['id'] == str(id))]
        crossings = str(insertion_table_srr['unique_crossings'].values[0]).split('|')

        # Calculate locus
        srr_id = insertion_table_srr[insertion_table_srr['id'] == id]
        chromosome = str(srr_id['chr'].values[0])
        start = str(srr_id['start'].values[0] - padding_left)
        stop = str(srr_id['stop'].values[0] + padding_right)
        locus = str(chromosome + ':' + start + '-' + stop)

        recipient_tracks = []
        # Create to-recipient alignment tracks
        for crossing in crossings:
            recipient_tracks.append(
                {
                    # Docs: 
                    # https://github.com/igvteam/igv.js/wiki/Alignment-Track
                    # https://github.com/igvteam/igv.js/blob/c6773940de86cf2938f40feec86d1905a866deba/js/bam/bamTrack.js
                    'name': crossing,
                    'url': batch + srr + '-to-' + recipient + '_' + crossing + '_filtered_id' + id + '.bam',
                    'indexURL': batch + srr + '-to-' + recipient + '_' + crossing + '_filtered_id' + id + '.bam.bai', 
                    'format': 'bam', 'type': 'alignment', 'showSoftClips': 'true', 'alignmentRowHeight' : '14',
                    'viewAsPairs': 'true', 'showMismatches' : 'true', 'colorBy' : 'strand', 
                    'coverageColor' : 'rgb(210, 100, 102)'
                }
            )

        donor_tracks = []
        # Create to-donor alignment tracks
        for crossing in crossings:
            donor_tracks.append(
                {
                    # Docs: 
                    # https://github.com/igvteam/igv.js/wiki/Alignment-Track
                    # https://github.com/igvteam/igv.js/blob/c6773940de86cf2938f40feec86d1905a866deba/js/bam/bamTrack.js
                    'name': crossing,
                    'url': batch + srr + '-to-' + donor + '_' + crossing + '_filtered_id' + id + '.bam',
                    'indexURL': batch + srr + '-to-' + donor + '_' + crossing + '_filtered_id' + id + '.bam.bai', 
                    'format': 'bam', 'type': 'alignment', 'showSoftClips': 'true', 'alignmentRowHeight' : '14',
                    'viewAsPairs': 'true', 'showMismatches' : 'true', 'colorBy' : 'strand', 
                    'coverageColor' : 'rgb(210, 100, 102)'
                }
            )

        donor_tracks.append(
            {
                'name': 'Annotations',
                'url': batch + donor + '.gff', 'isCompressed' : 'false',
                'displayMode': 'EXPANDED',
                'nameField': 'gene'
            }
        )
        # Return html
        return html.Div([
            html.Div([
                dashbio.Igv(
                    id = 'hg38',
                    genome = 'hg38',
                    locus = locus,
                    tracks = recipient_tracks
                )
            ], style={'width': '50%', 'display': 'inline-block'}),
            html.Div([
                dashbio.Igv(
                    id = donor,
                    reference = {
                        'id': donor,
                        'name': donor,
                        'fastaURL': batch + donor + '.fa',
                        'indexURL': batch + donor + '.fa.fai',
                        'tracks': donor_tracks,
                    },
                )
            ], style={'width': '50%', 'display': 'inline-block'})
        ])

if __name__ == '__main__':
    #app.run_server(debug=True)
    app.run_server(debug=True, host='0.0.0.0')