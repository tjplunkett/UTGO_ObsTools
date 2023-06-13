"""
skymap.py

Author: Viyaleta Apgar (modified by Thomas Plunkett)

Purpose: Generates a sky map for Greenhill Observatory, given a date of observation. This code defines methods that are called
in the main dashboard app (app.py).

Date: 10/06/23

"""

# Import necessary packages 
from datetime import datetime
from geopy import Nominatim
from tzwhere import tzwhere
from pytz import timezone, utc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from skyfield.api import Star, load, load_file, wgs84
from skyfield.data import hipparcos, mpc, stellarium
from skyfield.projections import build_stereographic_projection
from skyfield.constants import GM_SUN_Pitjeva_2005_km3_s2 as GM_SUN
from astropy.coordinates import SkyCoord
import astropy.units as u

import warnings
warnings.simplefilter("ignore", UserWarning)

def load_data():
    """
    A function to load the ephermeris (de421), hipparcos star catalog and
    constellations needed to generate sky map
    """
    # de421 shows position of earth and sun in space
    eph = load_file("de421.bsp")

    # hipparcos dataset
    with load.open("hip_main.dat") as f:
        stars = hipparcos.load_dataframe(f)

    with load.open("constellationship.fab") as f:
        constellations = stellarium.parse_constellations(f)
        
    return eph, stars, constellations
        
# load celestial data
eph, stars, constellations = load_data()

def collect_celestial_data(when, RA, DEC):
    """
    A function to generate the celestial data needed for sky map and then
    project into a stereographic representation.
    
    params:
    when -  The date of the observation (datetime)
    RA - The right ascension of the desired target (J2000, in degrees)
    DEC - The declination of the desired target (J2000, in degrees)
    
    return:
    stars - The star positions in the sterographic projection
    edge_star1, edge_star2 - Constellation edges 
    target_x, target_y - The position of the target in the stereographic projection
    """
    # Get coordinates of observatory
    lat, long = -42.43, 147.3
    
    # Defined a sky coord and convert to hours and degrees for target
    c = SkyCoord(ra=RA, dec=DEC, unit=u.deg)
    RA = c.ra.hour
    DEC = c.dec.degree
    
    # Convert date string into datetime object
    dt = datetime.strptime(when, '%Y-%m-%d %H:%M')

    # Define datetime and convert to UTC based on location coordinates
    timezone_str = 'Australia/Tasmania'
    local = timezone(timezone_str)
    utc_dt = local.localize(dt, is_dst=None).astimezone(utc)

    # Define observer using location coordinates and UTC time
    t = load.timescale().from_datetime(utc_dt)
    observer = wgs84.latlon(latitude_degrees=lat, longitude_degrees=long).at(t)

    # An ephemeris on Sun and Earth positions.
    sun = eph['sun']
    earth = eph['earth']
    
    # And the constellation outlines list.
    edges = [edge for name, edges in constellations for edge in edges]
    edges_star1 = [star1 for star1, star2 in edges]
    edges_star2 = [star2 for star1, star2 in edges]

    
    # Define the angle and center the observation location by the angle
    position = observer.from_altaz(alt_degrees=90, az_degrees=0)
    ra, dec, distance = observer.radec()
    center_object = Star(ra=ra, dec=dec)

    # Build the stereographic projection
    center = earth.at(t).observe(center_object)
    projection = build_stereographic_projection(center)
    field_of_view_degrees = 180.0

    # Compute the x and y coordinates based on the projection
    star_positions = earth.at(t).observe(Star.from_dataframe(stars))
    stars['x'], stars['y'] = projection(star_positions)
    target_x, target_y = projection(earth.at(t).observe(Star(ra_hours=RA, dec_degrees=DEC)))
    
    return stars, edges_star1, edges_star2, target_x, target_y
    
    
def create_star_chart(when, chart_size, max_star_size, RA, DEC, n_clicks):
    """
    Function to generate the star chart/sky map for GHO
    
    when - The datetime of the observation
    chart_size - Size of the chart
    max_star_size - Maximum size of stars to be plotted
    RA, DEC - The coordinates of the target in degrees
    n_clicks - The number of times the update button has been clicked in the dashboard
    """
    stars, edges_star1, edges_star2, tar_x, tar_y = collect_celestial_data(when, RA, DEC)
    
    # Define the number of stars and brightness of stars to include
    limiting_magnitude = 10
    bright_stars = (stars.magnitude <= limiting_magnitude)
    magnitude = stars['magnitude'][bright_stars]
    marker_size = max_star_size * 10 ** (magnitude / -2.5)
    
    # Calculate the constellation lines
    xy1 = stars[['x', 'y']].loc[edges_star1].values
    xy2 = stars[['x', 'y']].loc[edges_star2].values
    lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)
    
    fig, ax = plt.subplots(figsize=(chart_size, chart_size))
    
    # Add circle patch
    border = plt.Circle((0, 0), 1, color='#041A40', fill=True)
    ax.add_patch(border)

    # Draw the stars and constellation
    ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars],
           s=marker_size, color='white', marker='.', linewidths=0,
           zorder=2)
    ax.scatter(tar_x, tar_y, s=200, color = 'r', marker='x')
    ax.add_collection(LineCollection(lines_xy, colors='#ffff', linewidths=0.15))
    
    # Remove points outside the circle
    horizon = Circle((0, 0), radius=1, transform=ax.transData)
    for col in ax.collections:
        col.set_clip_path(horizon)

    # Finally, add other settings
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    plt.axis('off')
    when_datetime = datetime.strptime(when, '%Y-%m-%d %H:%M')
    plt.title(f"Observation Location: Bisdee Tier, Time: {when_datetime.strftime('%Y-%m-%d %H:%M')}", loc='right',color = 'white', fontsize=10)
    plt.tight_layout()
    
    # Check the n_clicks on web to allow for map to update (this is a hacky solution...)
    if n_clicks == 0:
        filename = "assets\BisdeeTier_SkyMap_base.png" 
    elif n_clicks % 2 == 0:
        filename = "assets\BisdeeTier_SkyMap_e.png" 
    else:
        filename = "assets\BisdeeTier_SkyMap_o.png" 
    
    # Save the figure and return the file name
    plt.savefig(filename, format='png')
    return filename
    
# call the function above
if __name__=='__main__':
    when = '2023-05-29 00:00'
    chart_size=12
    max_star_size=200
    skymap = create_star_chart(when, chart_size, max_star_size)
    print(skymap)