import dash
import dash_table as dt
import dash_auth
import pandas as pd

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
import dash_bio as dashbio
import dash_core_components as dcc
import dash_html_components as html
import os, glob

donor_names = [ 
    { 'value' : 'hpv16', 'label' : 'hpv16' }, 
    { 'value' : 'hpv18', 'label' : 'hpv18' },
    { 'value' : 'hiv1', 'label' : 'hiv1' }
]
recipient_names = [ { 'value' : 'UCSChg38', 'label' : 'UCSChg38' } ]

outputs_path_hpv16 = 'http://192.168.86.200:8080/hpv16_rnaseq/'
insertion_table_hpv16 = pd.read_csv(
    outputs_path_hpv16 + 'hpv16_to_hg38_table.csv', 
    dtype={'srr': str, 'id': str }).drop(columns=['sequence'])
outputs_path_hpv18 = 'http://192.168.86.200:8080/hpv18_rnaseq/'
insertion_table_hpv18 = pd.read_csv(
    outputs_path_hpv18 + 'hpv18_to_hg38_table.csv',
     dtype={'srr': str, 'id': str }).drop(columns=['sequence'])
outputs_path_hiv1 = 'http://192.168.86.200:8080/hiv1_rnaseq/'
insertion_table_hiv1 = pd.read_csv(
    outputs_path_hiv1 + 'hiv1_to_hg38_table.csv', 
    dtype={'srr': str, 'id': str }).drop(columns=['sequence'])

VALID_USERNAME_PASSWORD_PAIRS = {
    'dash': 'r0tt3nch1ck3n'
}

app = dash.Dash(__name__)

auth = dash_auth.BasicAuth(
    app,
    VALID_USERNAME_PASSWORD_PAIRS
)

app.layout = html.Div([
    dcc.Tabs(
        id='tab-selection', 
        value='tab-1', 
        children=[
            dcc.Tab(
                label = 'Table of Available Samples', 
                value = 'tab-1',
                children = [
                    html.Div([
                        html.Div([
                            dcc.Dropdown(
                                id = 'default-igv-donor-select',
                                options = donor_names,
                                value = 'hpv16'
                            )
                        ], style={'width': '25%', 'display': 'inline-block'}),
                        
                        html.Div([
                            dcc.Dropdown(
                                id = 'default-igv-recipient-select',
                                options = recipient_names,
                                value = 'UCSChg38'
                            )
                        ], style={'width': '25%', 'display': 'inline-block'}),
                        
                        html.Div([
                            dcc.Dropdown(
                                id = 'default-igv-srr-select'
                            )
                        ], style={'width': '25%', 'display': 'inline-block'}),

                        html.Div([
                            dcc.Dropdown(
                                id = 'default-igv-id-select',
                                value = '1'
                            )
                        ], style={'width': '25%', 'display': 'inline-block'})
                    ]),
                ]   
            ),
            dcc.Tab(
                label = 'IGV Visualization', 
                value = 'tab-2'),
            ]
        ),
    
    html.Div(id='tabs-content')
])

################
# Donor change > Update the SRRs and default SRR
@app.callback(
    Output('default-igv-srr-select', 'options'),
    Output('default-igv-srr-select', 'value'),
    Input('default-igv-donor-select', 'value')
)
def update_options(donor):
    if not donor:
        raise PreventUpdate
    else:
        if donor == 'hpv16':
            insertion_table = insertion_table_hpv16
        elif donor == 'hpv18':
            insertion_table = insertion_table_hpv18
        elif donor == 'hiv1':
            insertion_table = insertion_table_hiv1
        # Extract ids for SRR and build dictionary from it
        srrs = insertion_table['srr'].drop_duplicates()
        return [{ 'value' : srr, 'label' : srr } for srr in srrs], srrs[0]

###################
# Donor or SRR change > Update the SRR IDs
@app.callback(
    Output('default-igv-id-select', 'options'),
    Input('default-igv-donor-select', 'value'),
    Input('default-igv-srr-select', 'value')
)
def update_options(donor, srr):
    if not donor or not srr:
        raise PreventUpdate
    else:
        if donor == 'hpv16':
            insertion_table = insertion_table_hpv16
        elif donor == 'hpv18':
            insertion_table = insertion_table_hpv18
        elif donor == 'hiv1':
            insertion_table = insertion_table_hiv1
        # Extract ids for SRR and build dictionary from it
        insertion_table_srr = insertion_table[insertion_table['srr'] == srr]
        ids = insertion_table_srr['id']
        return [{ 'value' : id, 'label' : id } for id in ids]

###########################
# Return the IGV components
@app.callback(
    Output('tabs-content', 'children'),
    Input('tab-selection', 'value'),
    Input('default-igv-srr-select', 'value'),
    Input('default-igv-id-select', 'value'),
    Input('default-igv-donor-select', 'value'),
    Input('default-igv-recipient-select', 'value')
)

def render_content(tab, srr, id, donor, recipient):
    if donor == 'hpv16':
        insertion_table = insertion_table_hpv16
        outputs_path = outputs_path_hpv16
    elif donor == 'hpv18':
        insertion_table = insertion_table_hpv18
        outputs_path = outputs_path_hpv18
    elif donor == 'hiv1':
        insertion_table = insertion_table_hiv1
        outputs_path = outputs_path_hiv1

    if tab == 'tab-1':
        return html.Div([
            dt.DataTable(
                id = 'tbl', data = insertion_table.to_dict('records'),
                columns = [{'name': i, 'id': i} for i in insertion_table.columns],
                style_cell = {'textAlign': 'left'}
            )
        ])
    elif tab == 'tab-2':

        # Defaults
        padding_left = 150
        padding_right = 150

        # Extract track names
        insertion_table_srr = insertion_table[(insertion_table['srr'] == str(srr)) & (insertion_table['id'] == str(id))]
        crossings = str(insertion_table_srr['unique_crossings'].values[0]).split('|')

        # Calculate locus
        insertion_table_srr_id = insertion_table_srr[insertion_table_srr['id'] == id]
        chromosome = str(insertion_table_srr_id['chr'].values[0])
        start = str(insertion_table_srr_id['start'].values[0] - padding_left)
        stop = str(insertion_table_srr_id['stop'].values[0] + padding_right)
        locus = str(chromosome + ':' + start + '-' + stop)

        repicient_tracks = []
        # Create to-recipient alignment tracks
        for crossing in crossings:
            repicient_tracks.append(
                {
                    # Docs: 
                    # https://github.com/igvteam/igv.js/wiki/Alignment-Track
                    # https://github.com/igvteam/igv.js/blob/c6773940de86cf2938f40feec86d1905a866deba/js/bam/bamTrack.js
                    'name': crossing,
                    'url': outputs_path + srr + '-to-' + recipient + '_' + crossing + '_filtered.bam',
                    'indexURL': outputs_path + srr + '-to-' + recipient + '_' + crossing + '_filtered.bam.bai', 
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
                    'url': outputs_path + srr + '-to-' + donor + '_' + crossing + '_filtered.bam',
                    'indexURL': outputs_path + srr + '-to-' + donor + '_' + crossing + '_filtered.bam.bai', 
                    'format': 'bam', 'type': 'alignment', 'showSoftClips': 'true', 'alignmentRowHeight' : '14',
                    'viewAsPairs': 'true', 'showMismatches' : 'true', 'colorBy' : 'strand', 
                    'coverageColor' : 'rgb(210, 100, 102)'
                }
            )
        
        if donor == 'hpv16':
            donor_tracks.append(
                {
                    'name': 'Annotations',
                    'url': outputs_path + 'hpv16.gff', 'isCompressed' : 'false',
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
                    tracks = repicient_tracks
                )
            ], style={'width': '50%', 'display': 'inline-block'}),
            html.Div([
                dashbio.Igv(
                    id = donor,
                    reference = {
                        'id': donor,
                        'name': donor,
                        'fastaURL': outputs_path + donor + '.fa',
                        'indexURL': outputs_path + donor + '.fa.fai',
                        'tracks': donor_tracks,
                    },
                )
            ], style={'width': '50%', 'display': 'inline-block'})
        ])

if __name__ == '__main__':
    #app.run_server(debug=True)
    app.run_server(debug=True, host='0.0.0.0')