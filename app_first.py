from dash import Dash, html, dash_table, dcc, callback, Output, Input
import pandas as pd
import plotly.express as px
#different options to style
#import dash_design_kit as ddk
import dash_bootstrap_components as dbc

#dcc is dash core components and used for render interactive graphs. Plotly express is used to build the graph
#callback, output,input are used for user interactivity
#example app for connecting to data

#read in data
df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/gapminder2007.csv')

app = Dash()
#now using dbc 
external_stylesheets = [dbc.themes.CERULEAN]
app = Dash(__name__, external_stylesheets=external_stylesheets)

#use the list format, dash_table is used for tables
#and we can add a graph as well
# app.layout = [
#     html.Div(children= 'My first app with Dash'),
#     dash_table.DataTable(data =df.to_dict('records'), page_size = 10),
#     dcc.Graph(figure = px.histogram(df, x = 'continent', y = 'lifeExp', histfunc = 'avg'))
# ]

#make a different type of graph with callbacks/controls
#radioitem and graph are given IDs, which are used in the callback
# app.layout = [
#     html.Div(children= 'My first app with Data, Graph, and Controls'),
#     html.Hr(),
#     dcc.RadioItems(options = ['pop', 'lifeExp','gdpPercap'], value = 'lifeExp', id = 'controls-and-radio-item'),
#     dash_table.DataTable(data =df.to_dict('records'), page_size = 6),
#     dcc.Graph(figure = {}, id ='controls-and-graph')
# ]

#using dbc
#uses a container to wrap everything, similar to Rshiny when organizing the structure of where to put graphs and tables
#row for title, then row for radio buttons, the next row has two columns - first is the table and second is the graph
app.layout = dbc.Container([
    dbc.Row([
        html.Div('My First App with Data, Graph, and Controls', className="text-primary text-center fs-3")
    ]),

    dbc.Row([
        dbc.RadioItems(options=[{"label": x, "value": x} for x in ['pop', 'lifeExp', 'gdpPercap']],
                       value='lifeExp',
                       inline=True,
                       id='radio-buttons-final')
    ]),

    dbc.Row([
        dbc.Col([
            dash_table.DataTable(data=df.to_dict('records'), page_size=12, style_table={'overflowX': 'auto'})
        ], width=6),

        dbc.Col([
            dcc.Graph(figure={}, id='fig_1')
        ], width=6),
    ]),
    dbc.Row([
        dcc.Graph(figure = {}, id = "fig_2")
    ])

], fluid=True)

#Controls and callbacks are needed to give user options
@callback(
    [Output(component_id='fig_1', component_property='figure'),Output(component_id='fig_2', component_property='figure')],
    Input(component_id='radio-buttons-final', component_property='value')
)
# @callback(
# Output(component_id='fig_2', component_property='figure'),
# Input(component_id='radio-buttons-final', component_property='value')
# )


def update_graph(col_chosen):
    fig1 = px.histogram(df, x= 'continent', y = col_chosen, histfunc= 'avg')
    fig2 = px.histogram(df, x= 'country', y = col_chosen, histfunc= 'avg')
    return [fig1,fig2]



if __name__ == '__main__':
    app.run(debug = True)