import networkx as nx
import numpy as np
#import matplotlib.pyplot as plt

#import matplotlib.pyplot as plt
#import matplotlib.cm as cmx
from math import ceil, floor
#import matplotlib.pyplot as plt
import numpy as np
# import matplotlib.colors as colors 
import time
import random
#import sklearn.metrics as mt
import re
import networkx.algorithms.community as nx_comm



import dash_cytoscape as cyto  # pip install dash-cytoscape==0.2.0 or higher
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Output, Input,State
import pandas as pd  # pip install pandas
import plotly.express as px
import math
from dash import no_update
from dash.exceptions import PreventUpdate
#from utils import *
import dash
import json
from networkx.readwrite import json_graph
import dash_daq as daq

import base64
import datetime
import io

import dash_table
from dash.dependencies import Input, Output, State


def create_graph(nodes,edges,node_emphasize=True):
    node_label={}
    G = nx.MultiGraph()
    if node_emphasize==True:
        #G.add_nodes_from(list(nodes["#BIOGRID ID"]))
        #for nodeid in range(nodes.shape[0]):
        for row in nodes.iterrows():    
            node=row[1].loc["#BIOGRID ID"]
            OFFICIAL_SYMBOL=row[1].loc["OFFICIAL SYMBOL"]
            node_Entrez=row[1].loc["ENTREZ GENE ID"]
            G.add_node(node)
            G.nodes[node]["Official Symbol"]=OFFICIAL_SYMBOL
            G.nodes[node]["Entrez Gene ID"]=node_Entrez
            node_label[OFFICIAL_SYMBOL.lower()]=node
        for row in edges.iterrows():
            nodeA_bioGrid=row[1].loc["BioGRID ID Interactor A"]
            
            nodeB_bioGrid=row[1].loc["BioGRID ID Interactor B"]

            if nodeA_bioGrid in G and nodeB_bioGrid in G:
                #if nodeA_bioGrid not in G[nodeB_bioGrid] and nodeB_bioGrid not in G[nodeA_bioGrid]
                # Create edge attributes
                #G.add_edge(1, 2, weight=4.7)
                edge_bioGrid=row[1].loc["#BioGRID Interaction ID"]
                edge_Throughput=row[1].loc["Throughput"]
                edge_OTC=row[1].loc["Ontology Term Categories"]  
                edge_OTN=row[1].loc["Ontology Term Names"] 
                edge_OTQN=row[1].loc["Ontology Term Qualifier Names"] 
                G.add_edge(nodeA_bioGrid,nodeB_bioGrid, edge_bioGrid=edge_bioGrid,edge_Throughput=edge_Throughput,edge_OTC=edge_OTC,edge_OTN=edge_OTN,edge_OTQN=edge_OTQN)    
            
    else:    
        for row in edges.iterrows():
            nodeA_bioGrid=row[1].loc["BioGRID ID Interactor A"]
            nodeA_Entrez=row[1].loc["Entrez Gene Interactor A"]
            nodeA_Symbol=row[1].loc["Official Symbol Interactor A"]

            nodeB_bioGrid=row[1].loc["BioGRID ID Interactor B"]
            nodeB_Entrze=row[1].loc["Entrez Gene Interactor B"]
            nodeB_Symbol=row[1].loc["Official Symbol Interactor B"]

            if nodeA_bioGrid not in G:
                G.add_node(nodeA_bioGrid)
                G.nodes[nodeA_bioGrid]["Entrez Gene ID"]=nodeA_Entrez
                G.nodes[nodeA_bioGrid]["Official Symbol"]=nodeA_Symbol
            if nodeB_bioGrid not in G:
                G.add_node(nodeB_bioGrid)
                G.nodes[nodeB_bioGrid]["Entrez Gene ID"]=nodeB_Entrze
                G.nodes[nodeB_bioGrid]["Official Symbol"]=nodeB_Symbol

            #if nodeA_bioGrid not in G[nodeB_bioGrid] and nodeB_bioGrid not in G[nodeA_bioGrid]
            # Create edge attributes
            #G.add_edge(1, 2, weight=4.7)

            edge_bioGrid=row[1].loc["#BioGRID Interaction ID"]
            edge_Throughput=row[1].loc["Throughput"]
            edge_OTC=row[1].loc["Ontology Term Categories"]  
            edge_OTN=row[1].loc["Ontology Term Names"] 
            edge_OTQN=row[1].loc["Ontology Term Qualifier Names"] 

            G.add_edge(nodeA_bioGrid,nodeB_bioGrid, edge_bioGrid=edge_bioGrid,edge_Throughput=edge_Throughput,edge_OTC=edge_OTC,edge_OTN=edge_OTN,edge_OTQN=edge_OTQN)
            
    return G,node_label

def get_communities(G):
    c=0
    communities={}
    com_found=nx.community.greedy_modularity_communities(G)
    for com in com_found:
        c+=1
        for node in com:
            communities[node]=c
    nx_comm.modularity(G,com_found)
    return communities, nx_comm.modularity(G,com_found)

def create_elements(G,node_class,G_org):
    # node_class is a dictionary with key=node and value that is the class of node
    Cy_nodes=[
    # Nodes elements
    {'data': {'id': str(node), 'label': G.nodes[node]['Official Symbol'],'degree':G_org.degree(node)},
     'selectable': True,'locked': False,'grabbable': True,'classes': str(node_class[node])} for node in G
    ]

    Cy_edges=[
        {'data': {'source': str(edge[0]), 'target': str(edge[1]), 'label': str(edge[0])+">"+str(edge[1])}} for edge in G.edges()
    ]

    nodes_edges_Cy_with_label=Cy_nodes+Cy_edges
    Cy_nodes_without_label=[
    # Nodes elements
    {'data': {'id': str(node)},
     'selectable': True,'locked': False,'grabbable': True,'classes': str(node_class[node])} for node in G
    ]
    
    #Cy_edges_without_label=[
    #    {'data': {'source': str(edge[0]), 'target': str(edge[1]), 'label': str(edge[0])+">"+str(edge[1])}} for edge in G.edges()
    #]
    nodes_edges_Cy_without_label=Cy_nodes_without_label+Cy_edges
    
    return nodes_edges_Cy_with_label,nodes_edges_Cy_without_label

def create_unique_class(G):
    unique_class_node={}
    for node in G:
        unique_class_node[node]=0
    return unique_class_node

def best_classes(classes):
    best_comunities={}
    for i in list(classes.values()):
        best_comunities[i]=best_comunities.get(i,0)+1 
    best=[i[0] for i in sorted(best_comunities.items())[:23]]
    # Return only best 23 classes
    return best

def coloring (cids):
    # cids is a list of class number that should be colored. Note that the list should contain 23 class ids as the last color is meant for the rest (There are at most 24 colors) 
    col_swatch = px.colors.qualitative.Dark24
    style=[
            {
                'selector': '.'+str(cid),
                'style': {
                    'background-color': col_swatch[cid],
                }
            } for cid in cids
    ]
    return style

def degree_dist_df(G,q,selected_deg):
    # q is quartile
    degrees=np.array([G.degree(node) for node in G])
    limit=int(np.round(np.quantile(degrees,q))+2)
    degrees[degrees>limit]=limit
    hist=np.histogram(degrees,bins=limit)
    freq=hist[0]/np.sum(hist[0])
    degrees=hist[1]
    degrees=degrees[:-1]
    degrees=[np.round(i) for i in degrees]
    df=pd.DataFrame(None,columns=["Degree", "Freq","Selected"])
    df["Degree"]=degrees
    df["Freq"]=freq
    df["Selected"]=["No" for i in range(len(degrees))]
    if selected_deg is not None:
        df.iloc[selected_deg,2]="Yes"
    #fig = px.bar(df,x="Degree",y="Freq",color="Selected")
    return df

def get_centrality(G,centrality):
    res_normalized={} 
    if centrality=="betweenness":
        res=nx.centrality.betweenness_centrality(G)
    elif centrality=="eigenvector":
        res=nx.centrality.eigenvector_centrality(G)
    elif centrality=="closeness":
        res=nx.centrality.closeness.closeness_centrality(G)
    elif centrality=="clustering":
        res=nx.clustering(G)
    mn=np.min(list(res.values()))
    mx=np.max(list(res.values()))
    rg=mx-mn
    for node in res:
        res_normalized[node]=(res[node]-mn)/rg
    return res,res_normalized

def shortest_path(G,source,target):
    shortest_path_res=None
    edges=[]
    try:
        shortest_path_res=nx.shortest_paths.shortest_path(nx.Graph(G),source,target)
    except:
        shortest_path_res=None
    
    if shortest_path_res is not None:
        for i in range(len(shortest_path_res)-1):
            left=shortest_path_res[i]
            right=shortest_path_res[i+1]
            edges.append((left,right))  
    else:
        edges=None
    return edges

from networkx.readwrite import json_graph

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
#app = JupyterDash(__name__, external_stylesheets=external_stylesheets)
#app=dash.Dash(__name__, external_stylesheets=external_stylesheets,suppress_callback_exceptions=True)
#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__,external_stylesheets=[dbc.themes.SLATE])
#app = dash.Dash(__name__,external_stylesheets=external_stylesheets)
cyto.load_extra_layouts()

Msg="In the modern era when data exist everywhere, Information visualization is considered as an important analytical tool to display and process data for different aspects and it is quickly being recognized as an essential part of effective research communication. In this project for 'Data Visualization' course taught by Prof. John A. Lee at UClouvain, we were asked to provide an Information visualization dashboard to visualize a glioblastoma gene interaction dataset that includes both edges and nodes dataset separately. Furthermore, this interactional dashboard should have designed so that it can satisfy user's needs such as uploading new datasets with a specific format, modifying the color and size of the nodes and edges according to some metrics, depicting multiple views of the graph structure, and so on. In order to achieve this goal, we developed a user-friendly dashboard (called BioVisualizer) with Python 3 and Dash framework for building GUI. Dash is a powerful framework built specially for creating interactive data visualization apps. In order to create a Network visualization app using Dash, we utilized Dash Cytoscapethat is a network visualization component for Dash. Moreover, as we required to apply some network algorithms such as Community Detection, Shortest-path, we used python NetworkX library. This Dashbord was designed by Nima Farnoodian [nima.farnoodian@student.uclouvain.be] and Atefeh Bahrami [atefeh.bahrami@student.uclouvain.be] at EPL, Université catholique de Louvain, December 2021."

node_style= [{
        "selector": 'node',
        'style': {
            "content": "data(label)",
            "opacity": 0.5,
            "width": "data(size)",
            "height": "data(size)"
        }
    }]

edge_style= [{
        "selector": 'edge',
        'style': {
            "curve-style": "bezier",
            "opacity": 0.2,
            'width': 1
        }
    }]

#my_style_sheet=styles_color+edge_style+node_style
my_style_sheet=edge_style+node_style
def GCC_button():
    return dbc.Button("Plot all Connected Components", id="GCC", className="operation", color="light",size="sm")
def Compute_degree():
    return dbc.Button("Compute Degree", id="compute_degree", className="operation", color="link")
def Detect_community():
    return dbc.Button("Detect Comunities", id="detect_communities", className="operation" , color="link")
def Compute_betweenness():
    return dbc.Button("Compute Betweenness", id="compute_betweenness", className="operation", color="link")
def Compute_closeness():
    return dbc.Button("Compute Closeness", id="compute_closeness", className="operation", color="link")
def Compute_eigenvector():
    return dbc.Button("Compute Eigenvector", id="compute_eigenvector", className="operation", color="link")
def Compute_clustering():
    return dbc.Button("Compute Clustering", id="compute_clustering", className="operation", color="link")
def apply_coloring():
    return dbc.Button("Apply Coloring", id="apply_coloring", className="operation",size="sm")

def Compute_shortest_path():
    return dbc.Button("Compute Shortest_path", id="compute_shortest", className="operation" , color="link")

def build_button():
    return dbc.Button("Build Graph", id="build", className="operation",outline=True, color="primary")


def upload_edge():
    return dcc.Upload(
        id='upload-edge',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Edge File')
        ]),
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True
    )

def upload_node():
    return dcc.Upload(
        id='upload-node',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Node File')
        ]),
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True
    )

def slide_size_edge():
    return daq.Knob(id="size_slide_edge",
    label="Edge Size",
    value=1,
    #color={"gradient":True,"ranges":{"green":[0,5],"yellow":[5,9],"red":[9,10]}},
    max=5,
    scale={'start':0, 'labelInterval': 1, 'interval': 1}
    )

def slide_size_node():
    return daq.Knob(id="size_slide",
    label="Node Size",
    value=5,
    #color={"gradient":True,"ranges":{"green":[0,5],"yellow":[5,9],"red":[9,10]}},
    max=15,
    scale={'start':0, 'labelInterval': 2, 'interval': 2}
    )

def node_centrality():
    return dcc.Dropdown(
            id='centrality',
            value='None',
            clearable=False,
            options=[
                {'label': "None", 'value': "None"}           
            ])

def show_node_label():
    return dcc.Checklist(
        id="show_label",
        options=[
            {'label': "Node label", 'value': 'show'},
        ],
        value=['show'],
        labelStyle={'display': 'inline-block'}
    )

colorpick_style = {
    #'color':'white',
    'backgroundColor': '#272b30',
    'borderColor':'#272b30',
    'color':'#272b30',
    'theme':'#272b30'
}

tab_style = {
    #'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold',
    'color':'white',
    'backgroundColor': '#272b30',
    'borderColor':'#272b30'
}

tab_selected_style = {
    #'borderTop': '1px solid #d6d6d6',
    #'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#272b30',
    'color': 'white',
    'padding': '6px'
}

def node_color():
    return daq.ColorPicker(
    label='Node Color',
    value=dict(hex= '#8408E6') 
)

def node_label_color():
    return daq.ColorPicker(
    label='Node Label Color',
    value=dict(hex= '#1291ea')
)


def edge_color():
    return daq.ColorPicker(
    label='Edge color',
    value=dict(hex= '#ccd0d0')
)

def background_color():
    return daq.ColorPicker(
    label='Background color',
    value=dict(hex= '#272b30')
)

tabs_styles = {
    'height': '44px' , 
    #'width':'50px',
    'backgroundColor': '#272b30'
}
def color_tab():
    return dcc.Tabs(id="color_tab",style=tabs_styles,
        children=[
        dcc.Tab(label='Node', children=[node_color()],style=tab_style, selected_style=tab_selected_style),
        dcc.Tab(label='Lable', children=[node_label_color()],style=tab_style, selected_style=tab_selected_style), 
        dcc.Tab(label='Edge', children=[edge_color()],style=tab_style, selected_style=tab_selected_style),
        dcc.Tab(label='BG', children=[background_color()],style=tab_style, selected_style=tab_selected_style),
    ])

def build_error():
    return dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("Error while buiding graph")),
                dbc.ModalBody("An error appear occurred during building the graph. You may have uploaded the edge or node file incorrectly or one of them was not uploaded."),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close", id="close", className="ms-auto", n_clicks=0
                    )
                ),
            ],
            id="error_modal",
            is_open=False,
        )


######################################### nima layout#########################



# layouts: preset, random, cose, circular, grid, breadthfirst, concentric, external layouts
############################Draw figure##############################
def drawGraph():
    return  html.Div([
        dbc.Card(
            dbc.CardBody([
                cyto.Cytoscape(
                id='org-chart',
                autoungrabify=False,
                minZoom=0.05,
                maxZoom=3,
                layout={'name': 'random'},
                style={'width': '120%', 'height': '500px'},
                #elements=nodes_edges_Cy_wo_lbl,
                stylesheet=my_style_sheet
                )
            ])
        ),  
    ])
def degree_graph():
    return dcc.Graph(id='deg_dist_bar')

      
def slide_edge_opacity():
    return dcc.Slider(
                id="edge_opacity",
                min=0,
                max=1,
                step=0.05,
                marks={0: "0", 1: "1",},
                value=.5,
            )
def slide_node_opacity():
    return dcc.Slider(
                id="node_opacity",
                min=0,
                max=1,
                step=0.05,
                marks={0: "0", 1: "1",},
                value=.8,
            )
def quartile_dist():
    return dcc.Slider(
                id="quartile_dist",
                min=0,
                max=1,
                step=0.02,
                marks={0: "0", 1: "1",},
                value=.9,
            )

def layout_dropDown():
    return dcc.Dropdown(
                id='dpdn',
                value='random',
                clearable=False,
                options=[
                    {'label': name.capitalize(), 'value': name}
                    for name in  ['random','grid','circle','concentric','breadthfirst','cose','cose-bilkent','dagre','cola','klay','spread','euler']
                ])
#############1

app.layout = html.Div([
    dcc.Store(id='graph_type', data="allC", storage_type='memory'),
    dcc.Store(id='degree_dist_df', data=None, storage_type='memory'),
    dcc.Store(id='previously_selected_node', data=None, storage_type='memory'),
    dcc.Store(id='selected_node', data=None, storage_type='memory'),
    dcc.Store(id='betweenness_data', data=None, storage_type='memory'),
    dcc.Store(id='closeness_data', data=None, storage_type='memory'),
    dcc.Store(id='eigenvector_data', data=None, storage_type='memory'),
    dcc.Store(id='clustering_data', data=None, storage_type='memory'),
    dcc.Store(id='community_data', data=None, storage_type='memory'),
    dcc.Store(id='is_build', data=False),
    dbc.Carousel(
    items=[
        {"key": "1", "src": "https://s4.uupload.ir/files/s1_1e9p.jpg"},
        {"key": "2", "src": "https://s4.uupload.ir/files/s2_yo5.jpg"},
        {"key": "3", "src": "https://s4.uupload.ir/files/s3_ucxa.jpg"},
    ],
    controls=True,
    indicators=False,
    ),
    
    #html.H1("Bio Visualizer", style={'text-align': 'center', 'color' : 'white' }),
    #dbc.CardImg(src="https://s4.uupload.ir/files/h4_3sg4.jpg", top=True),
    dbc.Card(
        dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    
                    dbc.Row([
                        dbc.CardImg(src="https://s4.uupload.ir/files/logo_mi9z.png", top=True),
                    ],style={"width": "20rem", "align":'center'}),
                    
                    html.Br(),
                    
                    dbc.Row([
                    html.H1("BioVisualizer", style={'text-align': 'center', 'color' : '#B052FA' }),
                    html.H5("Glioblastoma Gene Interaction Visualization Dashboard", style={'text-align': 'center', 'color' : 'white' }),
                   
                    ],style={"width": "20rem"}),
                    html.Br(),
                    html.Br(),
                    dbc.Row([
                            html.P("Please upload your node and edge datasets here.", style={'text-align': 'center', 'color' : 'white' }),
                            upload_node(),
                            html.Div(id='output-node-data'),
                            upload_edge(),
                            html.Div(id='output-edge-data'),
                    
                        
                     ],style={"width": "20rem"}),
                    
                    dbc.Row([
                            build_button(),
                            build_error(),
                    ],style={"width": "8rem" ,'vertical-align': 'middle', "margin-left": "90px"}),
                    
                   
                    dbc.Row([
                            html.Div(id='print_info'),
                    ],style={"width": "25rem" ,'vertical-align': 'middle', "margin-left": "10px"}),
                    
                    html.Br(),
                    dbc.Row([
                      html.Hr(style={"size":150, "width":"200%" ,"color":"white", "margin-left": "90px"}) ,
                        
                     ],align='center' ,style={"width": "10rem"}),
                    
                    html.Br(),
                    
                    
                    html.Br(),
                    
                    
                    
                    dbc.Row([
                        
                    ]),
                    
                    html.Br(),
                    
                    dbc.Row([
                        
                        
                    ],style={"width": "15rem"}),
                    
                    
                        
                    dbc.Row([

                        #dbc.Button("Upload your data", color="primary"),
                    ],style={"width": "15rem"}),

                    #########button
                    




                    ###############
                     

                    
                       
                    
            ],width=2),
                                      
                
                dbc.Col([
                    dbc.Row([
                        dbc.Col([
                                html.P("Please Draw Network"),
                                GCC_button(),
                                                        
                        ],width=5),
                        
                        
                        dbc.Col([
                            html.P("Please Select layout"),
                            layout_dropDown(),
                                                        
                        ],width=3),
                        
                        dbc.Col([
                            html.P("Metrics Visualization mode"),
                            node_centrality(),
    
                       ],width=3),
                    ], style={"width": "73rem" ,'vertical-align': 'right', "margin-left": "0px"} ),
                       
                     dbc.Row([
                        
                        drawGraph(),
                        show_node_label(),
                     
                    ]),    
                    
                    
                                       
                ],width=7),    
                
                
                dbc.Col([
                    dbc.Row([
                        degree_graph(),
                        quartile_dist(),
                    dbc.Row([
                        #'/files/report-tutorial.pdf'
                        dbc.Offcanvas(
                            html.P(
                                [Msg,
                            html.Br(),
                            html.A('Click here to see the tutorial.', href='/files/report-tutorial.pdf'),
                            html.Br(),
                            html.A('Click here to download sample Nodes dataset.', href='/files/nodes.txt'),
                            html.Br(),
                            html.A('Click here to download sample Edges dataset.', href='/files/edges.txt')
                            ]
                            
                                
                            ),
                            
                            id="offcanvas",
                            title="BioVisualizer",
                            is_open=False,
                        ),
                    ],style={"width": "15rem"}),
                        html.Div(children='Node Info (Selected):'),
                        html.Div(id='node_info', children='Nothing selected'),
                        html.Br(),
                        html.Div(id="edge_selected",children='Edge Info:'),
                        html.Div(id='edge_info', children='Nothing selected'),
                        html.Br(),
                        html.Div(id='empty-div', children=''),
                        html.Br(),
                        html.Br(),
                        dbc.Button("About Us", id="open-offcanvas", n_clicks=0),

                        #check_color()
                        #drawFigure()
                    ]),
                    
                    
                ],width=3),  
                
                ], align='center'),
            
                dbc.Row([
                    
                    dbc.Col([
                        dbc.Row([
                        dbc.Accordion(
                            [  
                                                    
                            dbc.AccordionItem(
                                [
                                    
                                    Compute_degree(),
                                    Compute_betweenness(),
                                    Compute_closeness(),
                                    Compute_eigenvector(),
                                    Compute_clustering(),
                                ],
                                title="Compute metrics",
                            ),
                             dbc.AccordionItem(
                                [
                                    Detect_community(),
                                    dbc.Row([
                                    html.Div(id='print_info_community'),
                                    ],style={"width": "25rem" ,'vertical-align': 'middle', "margin-left": "1px"}),
                                    daq.PowerButton(
                                        id='community_show',
                                        label='Show',
                                        labelPosition='top',
                                        on=False,
                                        theme='dark',
                                        color='#FF5E5E'
                                    )
                                ],
                                title="Community Detection",
                            ),
                            dbc.AccordionItem(
                                [html.Div(children='Shortest path Based on Official Symbol:'),
                                dcc.Input(id="source",type="text",placeholder="Source Node"),
                                dcc.Input(id="target",type="text",placeholder="Target Node"),
                                Compute_shortest_path(),
                                daq.PowerButton(
                                        id='sp_show',
                                        label='Show',
                                        labelPosition='top',
                                        on=False,
                                        theme='dark',
                                        color='#FF5E5E'
                                    )
                                ],
                                title="Shortest-path")
                            ]
                                ),
 
                    ],style={"width": "20rem"}, align='center'),
                    
                    ],width=3),
                    
                    dbc.Col([
                        dbc.Row([
                        dbc.Col([
                             html.P("Edge Opacity", style={'textAlign': 'center'}),
                             slide_edge_opacity(),
                             slide_size_edge()
                        ],width=3),
                        
                        dbc.Col([
                             html.P("Node Opacity", style={'textAlign': 'center'}),
                            slide_node_opacity(),
                            slide_size_node()
                        ],width=3),
                            
                            dbc.Col([
                                 color_tab(),
                        ],width=3),
                            
                            dbc.Col([
                                html.Br(),
                                html.Br(),
                                html.Br(),
                                apply_coloring(),
                        ],width=2),
                           
                           

                    ]),
                    
                        
                    ],width=7),
                    
                    dbc.Col([
                    ],width=3),
                    
                ]),
    
   
        ]), color = 'dark'
    )
])




@app.callback(
    Output("offcanvas", "is_open"),
    Input("open-offcanvas", "n_clicks"),
    [State("offcanvas", "is_open")],
)
def toggle_offcanvas(n1, is_open):
    if n1:
        return not is_open
    return is_open




#########################################################




def parse_contents_edge(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' or 'txt' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')),header=0,sep='\t')
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        html.Hr(),  # horizontal line
        dcc.Store(id='edge-df', data=df.to_dict('records')),
        # For debugging, display the raw contents provided by the web browser
    ])


def parse_contents_node(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' or 'txt' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')),header=0,sep='\t')
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        html.Hr(),  # horizontal line
        dcc.Store(id='node-df', data=df.to_dict('records')),
        # For debugging, display the raw contents provided by the web browser
    ])


@app.callback(Output('upload-edge', 'children'),
              Input('upload-edge', 'contents'),
              State('upload-edge', 'filename'),
              State('upload-edge', 'last_modified'),
              prevent_initial_call=True)
def upload_edge(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents_edge(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children
    
    
@app.callback(Output('upload-node', 'children'),
              Input('upload-node', 'contents'),
              State('upload-node', 'filename'),
              State('upload-node', 'last_modified'),
              prevent_initial_call=True)
def upload_node(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents_node(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
                
        return children
     
    
@app.callback(Output('is_build', 'data'),
              Output('print_info', 'children'),
                Input('build', 'n_clicks'),
                State('node-df','data'),
                State('edge-df','data'),
                Input("close", "n_clicks"))

def build_graph(n, node_df, edge_df,close):
    print("I am here in build graph")
    if n is None:
        return no_update,no_update
    else:
        print("I am also here in build graph")
        global G
        global mn
        global mx
        global G_subgraph_multi
        global G_simple
        print("clicked")
        global node_label
        
        node_pd=pd.DataFrame(node_df)
        edge_pd=pd.DataFrame(edge_df)
        G,nb=create_graph(node_pd,edge_pd,node_emphasize=True)
        node_label=nb
        mn=np.mean([i[1] for i in list(G.degree())])
        mx=np.max([i[1] for i in list(G.degree())])

        # Creating G_subgraph_simple for future analysis
        largest_cc = max(nx.connected_components(G), key=len)
        G_subgraph_multi=G.subgraph(largest_cc)
        G_subgraph_multi=nx.MultiGraph(G_subgraph_multi)

        #global G_simple
        G_simple=nx.Graph(G)
        
        print(G.number_of_edges())
        print(G.number_of_nodes())
        #info="No of Nodes: " + str(G.number_of_nodes())+"|No of Edges: "+str(G.number_of_edges())
        
        print_info=[html.Br(),
                    html.Div(children="No of Nodes: " + str(G.number_of_nodes())),
                    html.Div(children="No of Edges: "+str(G.number_of_edges())),
                    html.Div(children="Average Degree: "+str(np.round((2*G.number_of_edges())/G.number_of_nodes(),4))),
                    html.Div(children="Is Graph Multiple? "+str(G.is_multigraph()))]
        
        return True,print_info
    
@app.callback(
    Output('edge_info', 'children'),
    Output('edge_selected','children'),
    Input('org-chart','tapEdgeData')
)
def show_edge_info(tap_edge):
    if tap_edge is None:
        raise PreventUpdate
    else:
        print("Tapped Edge: {}".format(tap_edge))
        source=int(tap_edge['source'])
        target=int(tap_edge['target'])
        print("test,",source,target)
        print("------------------------------------------------------------")
        print("------------------------------------------------------------")
        results=[html.Br()]
        edges=G[source][target]
        col1=[]
        col2=[]
        col3=[]
        col4=[]
        col5=[]
        edges=G[source][target]
        sourceOff=G.nodes[source]['Official Symbol']
        targetOff=G.nodes[target]['Official Symbol']
        for edge in edges:
            edge_bioGrid=edges[edge]['edge_bioGrid']
            col1.append(edge_bioGrid)
            edge_Throughput=edges[edge]['edge_Throughput']
            col2.append(edge_Throughput)
            edge_OTC=edges[edge]['edge_OTC']
            col3.append(edge_OTC)
            edge_OTN=edges[edge]['edge_OTN']
            col4.append(edge_OTN)
            edge_OTQN=edges[edge]['edge_OTQN']
            col5.append(edge_OTQN)
        data={"Bigrid ID":col1,"Throughput":col2,"Ontology Term Category":col3,"Ontology Term Name":col4,"Ontology Term Qualifier Name":col5}
        df=pd.DataFrame.from_dict(data)
        tbl=dash_table.DataTable(
            id='edge_information',
            columns=[{"name": i, "id": i} for i in df.columns],page_size=1,
            data=df.to_dict('records'),style_table={'overflowX': 'auto'},
                style_header={
                'backgroundColor': 'rgb(30, 30, 30)',
                'color': 'white'
            },
            style_data={
                'backgroundColor': 'rgb(50, 50, 50)',
                'color': 'white'
            },
        )
        edge_selected='Edge Info:' + "("+sourceOff+"-"+targetOff+")"
        return tbl,edge_selected



@app.callback(
    Output('org-chart', "stylesheet"),
    Output('org-chart', "style"),
    Input("edge_opacity", "value"),
    Input("size_slide_edge", "value"),
    Input("node_opacity", "value"),
    Input("org-chart", "stylesheet"),
    Input(component_id='show_label', component_property='value'),
    Input("color_tab", "children"),
    Input("apply_coloring", "n_clicks"),
    Input("community_show", "on"),
    State("community_data", "data"),
    Input("sp_show", "on")
)
def update_style(edge_value,size_slide_edge,node_value,whole_style_sheet,show,color_tab,color_click,com_on,community_data,sp_show):    
    if color_click is not None:
        node_color=color_tab[0]['props']['children'][0]['props']['value']['hex']
        node_label_color=color_tab[1]['props']['children'][0]['props']['value']['hex']
        edge_color=color_tab[2]['props']['children'][0]['props']['value']['hex']
        back_color=color_tab[3]['props']['children'][0]['props']['value']['hex']
    else:
        node_color='#1291ea'
        edge_color='#ccd0d0'
        back_color='#494A4A'
        node_label_color="#FFFFFF"
    if "show" in show:   
        node_style= [{
            "selector": 'node',
            'style': {
                "content": "data(label)",
                "opacity": node_value,
                "width": "data(size)",
                "height": "data(size)",
                "color":node_label_color,
                "background-color":node_color
            }
        }]
    else:
        node_style= [{
            "selector": 'node',
            'style': {
                "opacity": node_value,
            "width": "data(size)",
            "height": "data(size)",
            "color":node_label_color,
            "background-color":node_color
            }
        }]
    print(node_color)
    #############################     
    if sp_show==True:
        shortest_path_style=[{
            'selector': '.sp',
            'style': {
                'background-color': '#D51038',
                'line-color': '#D51038',
                "opacity": 1
                    }
                }]
        
        edge_style= [{
            "selector": 'edge',
            'style': {
                "curve-style": "bezier",
                "opacity": edge_value,
                'width': size_slide_edge*2}
            }] 
        my_new_style_sheet=edge_style+node_style+shortest_path_style
        
    else:
        edge_style= [{
            "selector": 'edge',
            'style': {
                "curve-style": "bezier",
                "opacity": edge_value,
                'width': size_slide_edge*2,
                "line-color":edge_color
            }
        }] 
        my_new_style_sheet=edge_style+node_style
        
    print(edge_color)
    #styles_color=[whole_style_sheet[0]]
    #my_new_style_sheet=styles_color+edge_style+node_style
    #my_new_style_sheet=edge_style+node_style
    
    if com_on:
        if community_data is not None:
            communities=json.loads(community_data)
            best_communities=best_classes(communities)
            com_color_style=coloring(best_communities)
            my_new_style_sheet=my_new_style_sheet+com_color_style
    style={'width': '100%', 'height': '500px', "background-color":back_color}
    
    return my_new_style_sheet,style

@app.callback(
    Output("org-chart", "elements"),
    Output('org-chart', 'layout'),
    Output('dpdn', 'value'),
    Output("GCC", "children"),
    Output("graph_type", "data"),
    Output("print_info_community", "children"),
    Output("community_data", "data"),
    Input("GCC", "n_clicks"),
    Input('dpdn', 'value'),
    Input("graph_type", "data"),
    State("org-chart", "elements"),
    Input("centrality","value"),
    Input("size_slide","value"),
    Input('betweenness_data','data'),
    Input('closeness_data','data'),
    Input('eigenvector_data','data'),
    Input('clustering_data','data'),
    State('is_build', 'data'),
    Input("detect_communities","n_clicks"),
    State("source","value"),
    State("target","value"),
    Input("compute_shortest","n_clicks"),
    prevent_initial_call=True
)

def Change_element(n_GCC,layout_value,g_type,elements,centrality_option,scale,betweenness_data,closeness_data,eigenvector_data,clustering_data,is_build,n_community,source,target,path_n):
    print("I am here in Change_element")
    if is_build==False:
        raise PreventUpdate
    else:
        print("I am also here in Change_element")
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        else:
            lyout={
                'name': "random",
                'animate': True
            }
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            print("item_id:",button_id)
            if button_id=="betweenness_data" or button_id=="closeness_data" or button_id=="eigenvector_data" or button_id=="clustering_data":
                 raise PreventUpdate 
            elif button_id=="GCC":
                inp=int(n_GCC%2) 
                if inp==0:
                    print(button_id,"Clicked")
                    print("I am also here")
                    #G_final_2=G_final.copy()
                    largest_cc = max(nx.connected_components(G_simple), key=len)
                    G_subgraph=G.subgraph(largest_cc)
                    G_subgraph_simple=nx.Graph(G_subgraph)
                    classes=create_unique_class(G_subgraph_simple)
                    nodes_edges_Cy_lbl_new,nodes_edges_Cy_wo_lbl_new=create_elements(G_subgraph_simple,classes,G_subgraph_multi)
                    #G_final=G_subgraph_simple
                    caption="Plot all Connected Components"
                    graph_type="GCC"
                elif inp==1:            
                    #global G_final
                    print(button_id,"Clicked")
                    #G_final_3=G_simple.copy()
                    classes_simple=create_unique_class(G_simple)
                    nodes_edges_Cy_lbl_new,nodes_edges_Cy_wo_lbl_new=create_elements(G_simple,classes_simple,G)
                    #G_final=G_simple
                    caption="Plot only Giant Connected Component"
                    graph_type="allC"
                return nodes_edges_Cy_lbl_new,lyout,"random",caption,graph_type,no_update,no_update
            elif button_id=='dpdn':
                lyt={
                'name': layout_value,
                'animate': True
                }
                if layout_value == 'random':
                    return no_update,lyt,layout_value,no_update,no_update,no_update,no_update
                else:
                    return no_update,lyt,layout_value,no_update,no_update,no_update,no_update
                
            ######Communitiy Detection   
            elif button_id=="detect_communities":
                if g_type=="GCC":
                    no=G_subgraph_multi.number_of_nodes()
                    my_g_simple=nx.Graph(G_subgraph_multi) 
                    my_g_org=G_subgraph_multi
                elif g_type=="allC":
                    no=G.number_of_nodes()
                    my_g_simple=G_simple
                    my_g_org=G

                comunities,modularity=get_communities(my_g_org)
                
                nodes_edges_Cy_lbl_new,_=create_elements(my_g_simple,comunities,my_g_org)
                print("I am done with community detection")
                
                no_com=len(np.unique(list(comunities.values())))
                print_info=[html.Br(),
                    html.Div(children="Method: Louvain Algorithm"),
                    html.Div(children="No of Communities: " + str(no_com)),
                    html.Div(children="Modularity Score "+str(np.round(modularity,4)))]
                community_data = json.dumps(comunities)
                                             
                return nodes_edges_Cy_lbl_new,no_update,no_update,no_update,no_update,print_info,community_data
                
            ######End           
            #############Shortest_path 
            elif button_id=="compute_shortest":
                if source is not None and target is not None:
                    if g_type=="GCC":
                        no=G_subgraph_multi.number_of_nodes()
                        my_g_simple=nx.Graph(G_subgraph_multi) 
                    elif g_type=="allC":
                        no=G.number_of_nodes()
                        my_g_simple=G_simple
                        
                    source=node_label.get(source.lower(),False)
                    target=node_label.get(target.lower(),False)
                    if source and target:
                        s_path=shortest_path(my_g_simple,source,target)
                        print("Hey I am in shortest-path")
                        print(s_path)
                        if s_path is not None:
                            for el_id in range(len(elements)):
                                el=elements[el_id]
                                if 'source' in el['data'] and 'target' in el['data']:
                                    s=int(el['data']['source'])
                                    t=int(el['data']['target'])
                                    if 'classes' in el:
                                        el['classes']="null"
                                        #del[el['classes']]
                                    if (s,t) in s_path or (t,s) in s_path:
                                        print("found (s,t)", (s,t))
                                        print("found (t,s)", (t,s))
                                        el['classes']='sp'
                                        print(el)
                                        elements[el_id]=el
                            elements_new=elements
                        else:
                            elements_new=no_update
                    else:
                        elements_new=no_update
                return elements_new,no_update,no_update,no_update,no_update,no_update,no_update
            ##########################
            elif button_id=="centrality":
                if g_type=="GCC":
                    no=G_subgraph_multi.number_of_nodes()
                    my_g=G_subgraph_multi
                elif g_type=="allC":
                    no=G.number_of_nodes()
                    my_g=G
                if centrality_option != "None":
                    print("centrality_option",centrality_option)
                    if centrality_option=="degree":
                        for i in range(no):    
                            deg=my_g.degree(int(elements[i]['data']["id"])) 
                            de_norm=(deg-mn)/mx
                            elements[i]['data']['size']= int((1 + de_norm*10 )* scale) 
                        return elements,no_update,no_update,no_update,no_update,no_update,no_update

                    if centrality_option=="betweenness":
                        if betweenness_data is not None:
                            r_normalized=json.loads(betweenness_data)["r_normalized"]
                            #r,r_normalized=get_centrality(G_simple,"betweenness")
                            for i in range(no):    
                                #deg=my_g.degree(int(elements[i]['data']["id"]))
                                bet=r_normalized[elements[i]['data']["id"]]
                                elements[i]['data']['size']= int((1 + bet*10 )* scale) 
                            return elements,no_update,no_update,no_update,no_update,no_update,no_update
                        else:
                            return no_update,no_update,no_update,no_update,no_update,no_update,no_update

                    if centrality_option=="closeness":
                        if closeness_data is not None:
                            r_normalized=json.loads(closeness_data)["r_normalized"]
                            #r,r_normalized=get_centrality(G_simple,"betweenness")
                            for i in range(no):    
                                #deg=my_g.degree(int(elements[i]['data']["id"]))
                                bet=r_normalized[elements[i]['data']["id"]]
                                elements[i]['data']['size']= int((1 + bet*10 )* scale) 
                            return elements,no_update,no_update,no_update,no_update,no_update,no_update
                        else:
                            return no_update,no_update,no_update,no_update,no_update,no_update,no_update

                    if centrality_option=="eigenvector":
                        if eigenvector_data is not None:
                            r_normalized=json.loads(eigenvector_data)["r_normalized"]
                            #r,r_normalized=get_centrality(G_simple,"betweenness")
                            for i in range(no):    
                                #deg=my_g.degree(int(elements[i]['data']["id"]))
                                bet=r_normalized[elements[i]['data']["id"]]
                                elements[i]['data']['size']= int((1 + bet*10 )* scale) 
                            return elements,no_update,no_update,no_update,no_update,no_update,no_update
                        else:
                            return no_update,no_update,no_update,no_update,no_update,no_update,no_update

                    if centrality_option=="clustering":
                        if clustering_data is not None:
                            r_normalized=json.loads(clustering_data)["r_normalized"]
                            #r,r_normalized=get_centrality(G_simple,"betweenness")
                            for i in range(no):    
                                #deg=my_g.degree(int(elements[i]['data']["id"]))
                                bet=r_normalized[elements[i]['data']["id"]]
                                elements[i]['data']['size']= int((1 + bet*10 )* scale) 
                            return elements,no_update,no_update,no_update,no_update,no_update,no_update
                        else:
                            return no_update,no_update,no_update,no_update,no_update,no_update,no_update

                if centrality_option == "None":
                    for i in range(no):    
                        elements[i]['data']['size']=int(5 * scale)
                    return elements,no_update,no_update,no_update,no_update,no_update,no_update 
                return no_update,no_update,no_update,no_update,no_update,no_update,no_update 

            elif button_id=="size_slide":
                if g_type=="GCC":
                    no=G_subgraph_multi.number_of_nodes()
                    my_g=G_subgraph_multi
                elif g_type=="allC":
                    no=G.number_of_nodes()
                    my_g=G
                if centrality_option != "None":
                    if centrality_option=="degree":
                        for i in range(no):    
                            deg=my_g.degree(int(elements[i]['data']["id"])) 
                            de_norm=(deg-mn)/mx
                            elements[i]['data']['size']= int((1 + de_norm*10 )* scale)  
                        return elements,no_update,no_update,no_update,no_update,no_update,no_update 

                    if centrality_option=="betweenness":
                        if betweenness_data is not None:
                            r_normalized=json.loads(betweenness_data)["r_normalized"]
                            #r,r_normalized=get_centrality(G_simple,"betweenness")
                            for i in range(no):    
                                #deg=my_g.degree(int(elements[i]['data']["id"]))
                                bet=r_normalized[elements[i]['data']["id"]]
                                elements[i]['data']['size']= int((1 + bet*10 )* scale) 
                            return elements,no_update,no_update,no_update,no_update,no_update,no_update 
                        else:
                            return no_update,no_update,no_update,no_update,no_update,no_update,no_update 

                    if centrality_option=="closeness":
                        if closeness_data is not None:
                            r_normalized=json.loads(closeness_data)["r_normalized"]
                            #r,r_normalized=get_centrality(G_simple,"betweenness")
                            for i in range(no):    
                                #deg=my_g.degree(int(elements[i]['data']["id"]))
                                bet=r_normalized[elements[i]['data']["id"]]
                                elements[i]['data']['size']= int((1 + bet*10 )* scale) 
                            return elements,no_update,no_update,no_update,no_update,no_update,no_update 
                        else:
                            return no_update,no_update,no_update,no_update,no_update,no_update,no_update 

                    if centrality_option=="eigenvector":
                        if eigenvector_data is not None:
                            r_normalized=json.loads(eigenvector_data)["r_normalized"]
                            #r,r_normalized=get_centrality(G_simple,"betweenness")
                            for i in range(no):    
                                #deg=my_g.degree(int(elements[i]['data']["id"]))
                                bet=r_normalized[elements[i]['data']["id"]]
                                elements[i]['data']['size']= int((1 + bet*10 )* scale) 
                            return elements,no_update,no_update,no_update,no_update,no_update,no_update 
                        else:
                            return no_update,no_update,no_update,no_update,no_update,no_update,no_update 

                    if centrality_option=="clustering":
                        if clustering_data is not None:
                            r_normalized=json.loads(clustering_data)["r_normalized"]
                            #r,r_normalized=get_centrality(G_simple,"betweenness")
                            for i in range(no):    
                                #deg=my_g.degree(int(elements[i]['data']["id"]))
                                bet=r_normalized[elements[i]['data']["id"]]
                                elements[i]['data']['size']= int((1 + bet*10 )* scale) 
                            return elements,no_update,no_update,no_update,no_update,no_update,no_update 
                        else:
                            return no_update,no_update,no_update,no_update,no_update,no_update,no_update 
                if centrality_option == "None":
                    for i in range(no):    
                        elements[i]['data']['size']=int(5 * scale)
                    return elements,no_update,no_update,no_update,no_update,no_update,no_update 
                return no_update,no_update,no_update,no_update,no_update,no_update,no_update 


            
@app.callback(
    Output('deg_dist_bar','figure'),
    Output('node_info', 'children'),
    Output('selected_node','data'),
    Output("centrality","options"),
    Output('betweenness_data','data'),
    Output('closeness_data','data'),
    Output('eigenvector_data','data'),
    Output('clustering_data','data'),
    Input("compute_degree", "n_clicks"),
    Input("compute_betweenness", "n_clicks"),
    Input("compute_closeness", "n_clicks"),
    Input("compute_eigenvector", "n_clicks"),
    Input("compute_clustering", "n_clicks"),
    Input("graph_type","data"),
    Input('org-chart','tapNodeData'),
    Input('quartile_dist','value'),
    Input('selected_node','data'),
    Input("centrality","options"),
    State('is_build', 'data'),
    State('betweenness_data', "data"),
    State('closeness_data', "data"),
    State('eigenvector_data', "data"),
    State('clustering_data', "data"),
    State('community_data',"data"),
    prevent_initial_call=True
)
def compute_centrality(n,n_between,n_closness,n_eigen,n_clustering,g_type,tap_node,q,selected_node,options,is_build,bt,cl,ei,clu,cm):
    if is_build==False:
        raise PreventUpdate
    global G
    global mn
    global mx
    global G_subgraph_multi
    global G_simple
    ctx = dash.callback_context
    if not ctx.triggered:
        print("test")
        raise PreventUpdate
    else:
        print(q)
        if g_type=="GCC":
            my_G=G_subgraph_multi
        elif g_type=="allC":
            my_G=G
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if button_id=="compute_degree":
            df=degree_dist_df(my_G,q,None)
            fig=px.bar(df,x="Degree",y="Freq",color="Selected")
            if {'label': "Degree Centrality", 'value': "degree"} not in options:
                options.append({'label': "Degree Centrality", 'value': "degree"})
            return fig.update_layout(template='plotly_dark',plot_bgcolor= 'rgba(0, 0, 0, 0)',paper_bgcolor= 'rgba(0, 0, 0, 0)'),no_update,no_update,options,no_update,no_update,no_update,no_update
        
        elif button_id=="compute_betweenness":
            r,r_normalized=get_centrality(G_simple,"betweenness")
            if {'label': "Betweenness Centrality", 'value': "betweenness"} not in options:
                options.append({'label': "Betweenness Centrality", 'value': "betweenness"})
                result={"r":r,"r_normalized":r_normalized}
                json_betweenness = json.dumps(result)
            else:
                json_betweenness=no_update
            return no_update,no_update,no_update,options,json_betweenness,no_update,no_update,no_update
        
        elif button_id=="compute_closeness":
            r,r_normalized=get_centrality(G_simple,"closeness")
            if {'label': "Closeness Centrality", 'value': "closeness"} not in options:
                options.append({'label': "Closeness Centrality", 'value': "closeness"})
                result={"r":r,"r_normalized":r_normalized}
                closeness_data = json.dumps(result)
            else:
                closeness_data=no_update
            return no_update,no_update,no_update,options,no_update,closeness_data,no_update,no_update
        
        elif button_id=="compute_eigenvector":
            r,r_normalized=get_centrality(G_simple,"eigenvector")
            if {'label': "Eigenvector Centrality", 'value': "eigenvector"} not in options:
                options.append({'label': "Eigenvector Centrality", 'value': "eigenvector"})
                result={"r":r,"r_normalized":r_normalized}
                eigenvector_data = json.dumps(result)
            else:
                eigenvector_data=no_update
            return no_update,no_update,no_update,options,no_update,no_update,eigenvector_data,no_update
        
        
        elif button_id=="compute_clustering":
            r,r_normalized=get_centrality(G_simple,"clustering")
            if {'label': "Clustering", 'value': "clustering"} not in options:
                options.append({'label': "Clustering", 'value': "clustering"})
                result={"r":r,"r_normalized":r_normalized}
                clustering_data = json.dumps(result)
            else:
                clustering_data=no_update
            return no_update,no_update,no_update,options,no_update,no_update,no_update,clustering_data
        
        elif button_id=='org-chart':
            dtTable={}
            node=int(tap_node["id"])
            Official_Symbol=my_G.nodes[node]['Official Symbol']
            dtTable['Official Symbol']=Official_Symbol
            Entrez_Gene_ID=my_G.nodes[node]['Entrez Gene ID']
            dtTable['Entrez Gene ID']=Entrez_Gene_ID
            node_deg=my_G.degree(node)
            dtTable['Degree']=node_deg
            if bt is not None:
                #print(json.loads(bt)["r"])
                dtTable['Betweenness']=np.round(json.loads(bt)["r"][str(node)],5)
            if cl is not None:
                dtTable['Closeness']=np.round(json.loads(cl)["r"][str(node)],5)
            if ei is not None:
                dtTable['Eigenvector']=np.round(json.loads(ei)["r"][str(node)],5)
            if clu is not None:
                dtTable['Clustering Coef.']=np.round(json.loads(clu)["r"][str(node)],5)
            if cm is not None:
                dtTable['Community ID']=json.loads(cm)[str(node)]
                
            tbl=pd.DataFrame.from_records([dtTable])
            ds_tbl=dash_table.DataTable(
                id='node_information',
                columns=[{"name": i, "id": i} for i in tbl.columns],page_size=1,
                data=tbl.to_dict('records'),style_table={'overflowX': 'auto'},
                    style_header={
                    'backgroundColor': 'rgb(30, 30, 30)',
                    'color': 'white'
                },
                style_data={
                    'backgroundColor': 'rgb(50, 50, 50)',
                    'color': 'white'
                },
            )
            df=degree_dist_df(my_G,q,None)
            if node_deg in df["Degree"]:
                df.loc[df.Degree == node_deg, 'Selected'] = "Yes"
            else:
                df.iloc[df.shape[0]-1,2] = "Yes"
            fig=px.bar(df,x="Degree",y="Freq",color="Selected")
            result="Official Symbol: "+ Official_Symbol+ ", Entrez Gene ID: "+ str(Entrez_Gene_ID)

            print("tapped Node: {}".format(tap_node))

            print(Official_Symbol,Entrez_Gene_ID)
            print("------------------------------------------------------------")
            print("------------------------------------------------------------")
            print(node)

            return fig.update_layout(template='plotly_dark',plot_bgcolor= 'rgba(0, 0, 0, 0)',paper_bgcolor= 'rgba(0, 0, 0, 0)'),ds_tbl,node,no_update,no_update,no_update,no_update,no_update
        
        
        elif button_id=="graph_type":
            df=degree_dist_df(my_G,q,None)
            fig=px.bar(df,x="Degree",y="Freq",color="Selected")
            return fig.update_layout(template='plotly_dark',plot_bgcolor= 'rgba(0, 0, 0, 0)',paper_bgcolor= 'rgba(0, 0, 0, 0)'),no_update,no_update,no_update,no_update,no_update,no_update,no_update
        
        elif button_id=="selected_node":
            return no_update,no_update,no_update,no_update,no_update,no_update,no_update,no_update
        
        elif button_id=="quartile_dist":
            if selected_node is not None:
                print("selected_node",selected_node)
                node_deg=my_G.degree(selected_node)
                df=degree_dist_df(my_G,q,None)
                if node_deg in df["Degree"]:
                    df.loc[df.Degree == node_deg, 'Selected'] = "Yes"
                else:
                    df.iloc[df.shape[0]-1,2] = "Yes"
                fig=px.bar(df,x="Degree",y="Freq",color="Selected")
            else:
                df=degree_dist_df(my_G,q,None)
                fig=px.bar(df,x="Degree",y="Freq",color="Selected")
            return fig.update_layout(template='plotly_dark',plot_bgcolor= 'rgba(0, 0, 0, 0)',paper_bgcolor= 'rgba(0, 0, 0, 0)'),no_update,selected_node,no_update,no_update,no_update,no_update,no_update


if __name__ == '__main__':
    app.run_server(port=8000, host='127.0.0.1')

#app.run_server(mode='jupyterlab')