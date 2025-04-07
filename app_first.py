from dash import Dash, html, dash_table, dcc
import pandas as pd
import plotly.express as px

#dcc is dash core components and used for render interactive graphs. Plotly express is used to build the graph
#example app for connecting to data

#read in data
df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/gapminder2007.csv')

app = Dash()

#use the list format, dash_table is used for tables
#and we can add a graph as well
app.layout = [
    html.Div(children= 'My first app with Dash'),
    dash_table.DataTable(data =df.to_dict('records'), page_size = 10),
    dcc.Graph(figure = px.histogram(df, x = 'continent', y = 'lifeExp', histfunc = 'avg'))
]

#Controls and callbacks are needed to give user options

if __name__ == '__main__':
    app.run(debug = True)