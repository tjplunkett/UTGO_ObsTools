"""
app.py

Author: Thomas Plunkett

Date: 10/06/23

Purpose:

Creates an interactive dashboard on the local host to help with observation planning using Dash
"""
# Import necessary packages 
from dash import Dash, html, dcc, callback, Output, Input, dash_table, State
from dash.exceptions import PreventUpdate
import plotly.express as px
from datetime import datetime, timedelta, date
from functions import *
from skymap import *
warnings.simplefilter("ignore", UserWarning)

# Calculate and define all necesary information for tables and dropdowns
default_date = datetime.combine(date.today(), datetime.min.time()) + timedelta(days=1) # Tomorrow midnight
m_phase, m_illum = moon_phase(default_date)
meridian_transit = calc_meridian_transit(219.874, -60.83, default_date)
moon_sep = calc_moon_sep(219.874, -60.83, default_date)
l, b = convert2gal(219.874, -60.83)
night_time, morning_time = find_night(default_date)
jd = get_jd(default_date)
tele_df = pd.DataFrame(columns = ['Telescopes'], data=['50cm', '1.3m'])
filter_df = pd.DataFrame(columns = ['Filters'], data=['i', 'r', 'g', 'V', 'B'])
seeing_df = pd.DataFrame(columns=['Seeing'], data = ['Good (1.5")', 'Average (2")', 'Poor (2.5")'])

# Create a dataframe with this information
night_dict = {'Quantity': ['Latitude (deg)', 'Longitude (deg)', 'Height (m)','Julian Day','RA (deg)', 'RA (hms)', 'DEC (deg)', 'DEC (dms)', 'l (deg)', 'b (deg)', 'Constellation', 'Ast. Twilight Ends', 'Ast. Morning Starts', 'Meridian Transit', 'Moon Phase (deg)', 'Perc. Illuminated', 'Moon Seperation (deg)'], 'Value': [-42.43,147.3, 646, jd, 219.874,'14:39:29.72',-60.83,'-60:49:56.00', round(l,3), round(b,3), 'Centaurus', str(night_time.iso), str(morning_time.iso), str(meridian_transit.iso), round(m_phase,3), round(m_illum,3), str(round(moon_sep,3))]}
night_df = pd.DataFrame(data=night_dict)

# Load the figures
sky_chart = create_star_chart(str(default_date.date())+" 00:00", 12, 200, 219.874, -60.93, 0)
airmass = plot_airmass(219.874, -60.83, default_date) # Alpha cent by default...
snrvsexp = plot_SNRvsTime(14, 'r', get_npix('Average (2")'), default_date)
magvstime = plot_MagLimvsTime(10, 'r', get_npix('Average (2")'), default_date)
snrvsmag = plot_SNRvsMag(30, 'r', get_npix('Average (2")'), default_date)

# Start the web app and define the layout

app = Dash(__name__)

app.layout = html.Div([
    #Title
    html.H1(children='Greenhill Observatory Observation Tools', style={'textAlign':'center'}), html.Br(),
    
    # Choose target and date
    html.Div([html.Label('Please input the target RA (in deg or sexigesimal): ', style={'width':'10vw'}), dcc.Input(id = 'RA', placeholder='14:39:29.72', type='text', style={'margin-right':'10px', 'width':'8vw'}),
    html.Label('Please input the target DEC (in deg or sexigesimal): ', style={'width':'10vw'}), dcc.Input(id = 'DEC', placeholder='-60:49:56.00', type='text', style={'margin-right':'10px', 'width':'8vw'}),
    html.Label('Please input the date of observation: '), dcc.DatePickerSingle(id='ObsDate', min_date_allowed=date(2022, 12, 1), max_date_allowed=date(2050, 12, 1), initial_visible_month=default_date, date=default_date)]), html.Br(),
    
    # Update button
    html.Button('Update', id='btn', style={'width':'75px', 'height':'30px', 'margin-left':'25px', 'textAlign':'center'}), html.Br(),
    
    # Sky map, coordinate conversion and moon phase
    html.Div([html.Img(id = 'skychart',src=sky_chart, width='600px', height='600px', style={'margin-left':'50px', 'margin-right':'25px'}), dash_table.DataTable(id = 'skytable', data=night_df.to_dict('records'), style_header={'fontWeight': 'bold', 'color':'#FFFFFF', 'backgroundColor':'#008B8B'}, style_table={'width':'650px', 'height':'600px', 'padding':'25px'}, style_cell={'font_size': '20px','text_align': 'center'})], style = {'display':'flex'}), html.Br(),
    
    # Airmass plotter
    html.H2(children = 'Airmass Calculator', style = {'textAlign':'left'}), dcc.Graph(id='Airmass', figure=airmass, style={'width':'96vw'}), html.Br(),
    
    # Exposure time calculator
    html.H2(children = 'Exposure Time Calculator', style = {'textAlign':'left'}), html.Label('Please select a telescope: ', style={'margin-right':'85px'}), html.Label('Please select a filter: ', style = {'margin-right':'110px'}), html.Label('What is the seeing like?'),
html.Div([dcc.Dropdown(tele_df['Telescopes'], '50cm', id='TeleSelector', style={'width':'200px', 'margin-right':'50px'}), dcc.Dropdown(filter_df['Filters'], 'r', id='FilterSelector', style={'width':'200px', 'margin-right':'50px'}), dcc.Dropdown(seeing_df['Seeing'], 'Average (2")', id='SeeingSelector', style={'width':'200px', 'margin-right':'50px'})], style = {'display':'flex'}), html.Br(),
    
html.Div([html.Div([html.Label('Magnitude: '), dcc.Input(id='sve_mag', value = 14), dcc.Graph(id='snrvexp_fig', figure=snrvsexp, style = {'width':'32vw'})]), html.Div([html.Label('SNR: '), dcc.Input(id='mvt_snr', value = 10), dcc.Graph(id='magvstime_fig', figure=magvstime, style = {'width':'32vw'})]), html.Div([html.Label('Exposure (s): '), dcc.Input(id='svm_exp', value = 30), dcc.Graph(id='snrvsmag_fig', figure=snrvsmag, style = {'width':'32vw'})])], style = {'display':'flex'})
        
])

# Define the interactivity 

@app.callback(
    Output(component_id='skytable', component_property='data'),
    Input(component_id='btn', component_property='n_clicks'),
    State('RA', 'value'),
    State('DEC', 'value'),
    State('ObsDate', 'date'))
def update_table(n_clicks, value1, value2, date_value):
    date_object = datetime.fromisoformat(date_value)
    if n_clicks is None:
        raise PreventUpdate
        
    # Recalculate values for the tables with new inputs
    ra_deg, dec_deg = convert2Deg(value1, value2)
    ra_sex, dec_sex = convert2Sex(value1, value2)
    l_new, b_new = convert2gal(ra_deg, dec_deg)
    m_phase, m_illum = moon_phase(date_object)
    meridian_transit = calc_meridian_transit(ra_deg, dec_deg, date_object)
    moon_sep = calc_moon_sep(ra_deg, dec_deg, date_object)
    constel = find_constellation(ra_deg, dec_deg)
    night_time_new, morning_time_new = find_night(date_object)
    jd = get_jd(date_object)
    
    # Put into dataframe
    night_dict = {'Quantity': ['Latitude (deg)', 'Longitude (deg)', 'Height (m)','Julian Day','RA (deg)', 'RA (hms)', 'DEC (deg)', 'DEC (dms)', 'l (deg)', 'b (deg)', 'Constellation', 'Ast. Twilight Ends', 'Ast. Morning Starts', 'Meridian Transit', 'Moon Phase (deg)', 'Perc. Illuminated', 'Moon Seperation (deg)'], 'Value': [-42.43, 147.3, 646, jd, round(ra_deg, 3), str(ra_sex), round(dec_deg, 3), str(dec_sex), round(l_new,3), round(b_new,3), constel, str(night_time_new.iso), str(morning_time_new.iso), str(meridian_transit.iso), round(m_phase, 3), round(m_illum, 3), str(round(moon_sep, 3))]}
    night_df = pd.DataFrame(data=night_dict)
    
    return night_df.to_dict('records')
    

@app.callback(
    Output(component_id='Airmass', component_property='figure'),
    Input(component_id='btn', component_property='n_clicks'),
    State('RA', 'value'),
    State('DEC', 'value'),
    State('ObsDate', 'date'))
def update_airmass(n_clicks, value1, value2, date_value):
    date_object = datetime.fromisoformat(date_value)
    if n_clicks is None:
        raise PreventUpdate
    ra_deg, dec_deg = convert2Deg(value1, value2)
    airmass = plot_airmass(float(ra_deg), float(dec_deg), datetime.combine(date_object, datetime.min.time()))
    return airmass

@app.callback(
    Output(component_id='skychart', component_property='src'),
    Input(component_id='btn', component_property='n_clicks'),
    State('RA', 'value'),
    State('DEC', 'value'),
    State('ObsDate', 'date'))
def update_skychart(n_clicks, value1, value2, date_value):
    date_object = datetime.fromisoformat(date_value)
    if n_clicks is None:
        raise PreventUpdate
    ra_deg, dec_deg = convert2Deg(value1, value2)
    skychart = create_star_chart(str(date_object.date())+" 00:00", 12, 200, ra_deg, dec_deg, n_clicks)
    return skychart

@app.callback(
    Output(component_id='snrvexp_fig', component_property='figure'),
    Input(component_id = 'FilterSelector', component_property = 'value'),
    Input(component_id='sve_mag', component_property='value'),
    Input(component_id = 'SeeingSelector', component_property = 'value'),
    State('RA', 'value'),
    State('DEC', 'value'),
    State('ObsDate', 'date'))
def update_snrvsexp(value1, value2, value3, value4, value5, date_value):
    date_object = datetime.fromisoformat(date_value)
    if value2 == '':
        raise PreventUpdate
    else:
        n_pix = get_npix(str(value3))
        snrvsexp = plot_SNRvsTime(float(value2), value1, n_pix, date_object)
        
    return snrvsexp

@app.callback(
    Output(component_id='magvstime_fig', component_property='figure'),
    Input(component_id = 'FilterSelector', component_property = 'value'),
    Input(component_id='mvt_snr', component_property='value'),
    Input(component_id = 'SeeingSelector', component_property = 'value'),
    State('RA', 'value'),
    State('DEC', 'value'),
    State('ObsDate', 'date'))
def update_snrvsexp(value1, value2, value3, value4, value5, date_value):
    date_object = datetime.fromisoformat(date_value)
    if value2 == '':
        raise PreventUpdate
    else:
        n_pix = get_npix(str(value3))
        magvstime = plot_MagLimvsTime(float(value2), value1, n_pix, date_object)
        
    return magvstime

@app.callback(
    Output(component_id='snrvsmag_fig', component_property='figure'),
    Input(component_id = 'FilterSelector', component_property = 'value'),
    Input(component_id='svm_exp', component_property='value'),
    Input(component_id = 'SeeingSelector', component_property = 'value'),
    State('RA', 'value'),
    State('DEC', 'value'),
    State('ObsDate', 'date'))
def update_snrvsexp(value1, value2, value3, value4, value5, date_value):
    date_object = datetime.fromisoformat(date_value)
    if value2 == '':
        raise PreventUpdate
    else:
        n_pix = get_npix(str(value3))
        snrvsmag = plot_SNRvsMag(float(value2), value1, n_pix, date_object)
        
    return snrvsmag
    
if __name__ == '__main__':
    app.run_server(debug=True, dev_tools_hot_reload=False)