"""
functions.py

Author: Thomas Plunkett

Date: 10/06/23

Purpose: 

Defines constants, telescope parameters and various functions needed for the interactive dashboard (app.py)
"""

# Import necessary packages
import numpy as np
import pandas as pd
import datetime as dt
from pytz import timezone
from skyfield import almanac
from skyfield.api import S, E, wgs84, load, Star
from skyfield.framelib import ecliptic_frame
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_constellation
from astropy.time import Time
from astroplan import Observer, FixedTarget
import plotly.express as px

# Define constants and telescope parameters
q_vec = [0.05, 0.1, 0.075, 0.2, 0.1, 0.075, 0.09, 0.15, 0.6, 0.9, 0.92, 0.9, 0.7, 0.4, 0.15]
wav_vec = [200, 220, 230, 240, 260, 280, 300, 340, 410, 500, 580, 660, 800, 900, 1000]
r = 50.8/2 # in cm
rn = 8.7 # e-
R = 34
T = 0.002
D = 0.14 # e/pix/sec
gain = 1.2 # e-/ADU
lati = -42.43
long = 147.3
hght = 646 # m

# Define some filter characteristics (all in nm) from Bessell 2005
i_bp = 123.0
i_eff = 743.9
q_i = np.interp(i_eff, wav_vec, q_vec)
i_zpt = 21.83
i_zpt_er = 0.04

r_bp = 115.0
r_eff = 612.2
q_r = np.interp(r_eff, wav_vec, q_vec)
r_zpt = 22.07
r_zpt_er = 0.07

g_bp = 128.0
g_eff = 463.9
q_g = np.interp(g_eff, wav_vec, q_vec)
g_zpt = 22.22
g_zpt_er = 0.04

V_bp = 84.0
V_eff = 544.8
q_V = np.interp(V_eff, wav_vec, q_vec)
V_zpt = 22.03
V_zpt_er = 0.06

B_bp = 89.0
B_eff = 436.1
q_B = np.interp(B_eff, wav_vec, q_vec)
B_zpt = 20.76
B_zpt_er = 0.06

# Define the observatory
bisdee_tier = EarthLocation(lat=lati*u.deg, lon=long*u.deg, height=hght*u.m)
GHO = Observer(location=bisdee_tier, name="Greenhill Observatory", timezone="Australia/Tasmania")

def convert2Deg(RA, DEC):
    """
    Function to convert sexigesimal RA, DEC to degrees
    
    params:
    
    RA - The right ascension of the target (either sex. or deg.)
    DEC - The declination of the target (either sex. or deg.)
    
    returns:
    con_ra, con_dec - The converted coordinates
    """
    # If : is present, assume a sexigesimal string and convert to degrees
    if ':' in str(RA) or ':' in str(DEC):   
        c = SkyCoord(ra=RA, dec=DEC, unit=(u.hourangle, u.deg))
        con_ra = c.ra.degree
        con_dec = c.dec.degree
        
    # If not present, assume already in degrees and return as floats
    else:
        con_ra = RA
        con_dec = DEC
        
    return float(con_ra), float(con_dec)

def convert2Sex(RA, DEC):
    """
    Function to convert degree RA, DEC to sexigesimal
    
    params:
    
    RA - The right ascension of the target (either sex. or deg.)
    DEC - The declination of the target (either sex. or deg.)
    
    returns:
    con_ra, con_dec - The converted coordinates
    """
    # If : is present, assume a sexigesimal string and keep
    if ':' in str(RA) or ':' in str(DEC):
        con_ra = RA
        con_dec = DEC
        
    # If not present, assume degrees and convert to hour angle and degree in sex. representation
    else:    
        c = SkyCoord(ra=RA, dec=DEC, unit = (u.deg, u.deg))
        con_ra = c.ra.to_string(unit=u.hourangle, sep=':')
        con_dec = c.dec.to_string(unit=u.degree, sep=':')    
        
    return con_ra, con_dec 

def convert2gal(RA, DEC):
    """
    Function to convert degree RA, DEC to sexigesimal
    
    params:
    
    RA - The right ascension of the target (either sex. or deg.)
    DEC - The declination of the target (either sex. or deg.)
    
    returns:
    l - The galatic longitude in degrees
    b - The galactic latitude in degrees
    """
    c = SkyCoord(ra=RA, dec=DEC, unit = (u.deg, u.deg))
    l = c.galactic.l.value
    b = c.galactic.b.value  
        
    return l, b

def get_fltrinfo(fltr):
    """
    Retrieves the zeropoint, effective wavelength, bandpass and quantum efficiency
    for a specific filter
    
    params:
    fltr - The desired filter
    """
    if fltr == 'i':
        zpt = i_zpt
        eff_wl = i_eff
        bp = i_bp
        q = q_i
    elif fltr == 'r':
        zpt = r_zpt
        eff_wl = r_eff
        bp = r_bp
        q = q_r   
    elif fltr == 'g':
        zpt = g_zpt
        eff_wl = g_eff
        bp = g_bp
        q = q_g 
    elif fltr == 'V':
        zpt = V_zpt
        eff_wl = V_eff
        bp = V_bp
        q = q_V   
    elif fltr == 'B':
        zpt = B_zpt
        eff_wl = B_eff
        bp = B_bp
        q = q_B
    else:
        print('Not a valid filter...')
           
    return zpt, eff_wl, bp, q

def calc_flux_density(mag, fltr):
    """
    A function to calculate the flux density of a point source given input mag and filter info.
    
    Params:
    mag - The magnitude in some filter
    Zpt - The magnitude zeropoint for desired filter
    eff_wl - The effective (central) wavelength of the filter
    
    Return:
    F - The flux density of the object
    """
    zpt, eff_wl, bp, q = get_fltrinfo(fltr)
    f = 10**((zpt-mag)/2.5)
    F = f/eff_wl
    return F

def calc_mag_limit(F, fltr):
    """
    A function to calculate a limiting magnitude given the flux density of a point source for a given filter.
    
    Params:
    F - The flux density in desired filter
    Zpt - The magnitude zeropoint for desired filter
    eff_wl - The effective (central) wavelength of the filter
    
    Return:
    mag - The limiting magnitude
    """
    zpt, eff_wl, bp, q = get_fltrinfo(fltr)
    f = F*eff_wl
    mag = -2.5*np.log10(f)+zpt
    return mag

def calc_E(F_obj, r, T, fltr):
    """
    A function to calculate the object signal rate of a point source given flux density 
    and telescope properties.
    
    Params:
    F_obj - The flux density of the object in a certain filter
    r - The radius of the telescope mirror
    T - The transmission efficiency of atmosphere, telescope and camera
    q - The quantum efficiency of the camera
    fltr - The desired filter (str)
    
    Return:
    E - The object signal rate
    """
    zpt, eff_wl, bp, q = get_fltrinfo(fltr)
    E = F_obj*np.pi*(r**2)*T*q*bp
    return E

def calc_Fobj(E, r, T, fltr):
    """
    A function to calculate the object signal rate of a point source given flux density 
    and telescope properties.
    
    Params:
    E - The object signal rate
    r - The radius of the telescope mirror
    T - The transmission efficiency of atmosphere, telescope and camera
    q - The quantum efficiency of the camera
    fltr - The desired filter
    
    Return:
    F_obj - The flux density of the object in a certain filter
    """
    zpt, eff_wl, bp, q = get_fltrinfo(fltr)
    F_obj = E/(np.pi*(r**2)*T*q*bp)
    return F_obj

def calc_SNR(E, t, n_pix, R, D, rn):
    """
    Function to calculate the theoretical signal to noise ratio for the 50cm telescope
    
    params:
    E - The object signal rate
    t - The exposure time
    n_pix - The number of pixels for the aperture
    R - The sky background 
    D - The dark current 
    rn - The read noise 
    """
    SNR = E*t/np.sqrt((E*t)+n_pix*(R*t+D*t+rn**2))
    return SNR

def invert_SNR(snr, t, n_pix, R, D, rn):
    """A function to calculate the object signal rate of a point source given a desired SNR
    and telescope properties.
    
    Params:
    snr - The signal to noise ratio
    t - The exposure time
    n_pix - The number of pixels for the aperture
    R - The sky background 
    D - The dark current 
    rn - The read noise
    
    Return:
    E - The object signal rate
    """
    sigma = (R*t+D*t+rn**2)
    E = ((snr**2)+np.sqrt((snr**4)+(4*n_pix*sigma*(snr**2))))/(2*t)
    return E

def plot_SNRvsTime(mag, fltr, n_pix):
    """
    Plot the signal-to-noise as a function of exposure time, given a target magnitude and filter
    
    params:
    mag - The magnitude of the target
    fltr - The filter to be used
    n_pix - The number of pixels 
    
    return:
    fig - The plotly figure
    """
    # Calculate the snr at times between 0.5 and 600 s
    t_vec = np.arange(0.5, 10*60, 0.5)
    f_density = calc_flux_density(mag, fltr)
    E = calc_E(f_density, r, 0.002, fltr)
    snr = calc_SNR(E, t_vec, n_pix, R, D, rn)    
    
    # Make the figure
    fig = px.line(x=t_vec, y=snr, labels={'x': 'Exposure Time (s)', 'y': 'SNR'}, title='Signal-To-Noise vs Exp. Time for GHO 50 cm', template='plotly_dark')
    fig.add_annotation(x=max(t_vec)/2, y=max(snr)/4,
            text="Mag = "+str(mag), showarrow=False)
    fig.add_annotation(x=max(t_vec)/2, y=max(snr)/2,
            text="Filter = "+fltr, showarrow=False)

    return fig

def plot_MagLimvsTime(snr, fltr, n_pix):
    """
    Plot the magnitude limit as a function of exposure time, given a target SNR and filter
    
    params:
    snr - The desired snr for the target
    fltr - The filter to be used
    
    return:
    fig - The plotly figure
    """
    # Calculate the limiting magnitude for times between 0.5 and 600 s
    t_vec = np.arange(0.5, 10*60, 0.5)
    E = invert_SNR(snr, t_vec, n_pix, R, D, rn)
    F = calc_Fobj(E,r, 0.002, fltr)
    mag_lim = calc_mag_limit(F, fltr)
    
    # Make the plot
    fig = px.line(x=t_vec, y=mag_lim, labels={'x': 'Exposure Time (s)', 'y': 'Magnitude Limit'}, title='Magnitude Limit vs Exp. Time for GHO 50 cm', template='plotly_dark')
    fig.add_annotation(x=max(t_vec)/2, y=max(mag_lim)/1.1,
            text="SNR = "+str(snr), showarrow=False)
    fig.add_annotation(x=max(t_vec)/2, y=max(mag_lim)/1.25,
            text="Filter = "+fltr, showarrow=False)

    return fig
    
def plot_SNRvsMag(exp_time, fltr, n_pix):
    """
    Plot the signal-to-noise of target magnitude, given a fixed exposure time and filter
    
    params:
    exp_time - The exposure time (s)
    fltr - The filter to be used
    
    return:
    fig - The plotly figure
    """
    # Calculate the limiting magnitude as a function of SNR between 0.1 and 1000
    snr_vec = np.arange(0.1, 1000, 0.1)
    E_vec = invert_SNR(snr_vec, exp_time, n_pix, R, D, rn)
    F_vec = calc_Fobj(E_vec,r, 0.002, fltr)
    mag_vec = calc_mag_limit(F_vec, fltr)
    
    # Make the plot
    fig = px.line(x=mag_vec, y=snr_vec, labels={'x': 'Magnitude', 'y': 'SNR'}, title='Signal-To-Noise vs Magnitude for GHO 50 cm', template='plotly_dark')
    fig.add_annotation(x=max(mag_vec), y=max(snr_vec)/4,
            text="Exp. Time = "+str(exp_time), showarrow=False)
    fig.add_annotation(x=max(mag_vec), y=max(snr_vec)/2,
            text="Filter = "+fltr, showarrow=False)

    return fig
    
def moon_phase(dt):
    """
    Find the moon phase and percentage illuminated given a datetime, with timezone given.
    
    params:
    dt - The datetime of the observation
    
    returns:
    phase - The moon phase in degrees
    percent - The percent illumination of the moon
    """
    # Try the easy Astroplan way (might need internet)
    try:
        phase = GHO.moon_phase(dt)
        percent = 100.0 * GHO.moon_illuminated(dt)
    
    # Try the harder Skyfield way, which won't require internet
    except:
        tz = timezone("Australia/Tasmania")
        dt_aware = tz.localize(dt)
        ts = load.timescale()
        t = ts.from_datetime(dt_aware)

        eph = load('de421.bsp')
        sun, moon, earth = eph['sun'], eph['moon'], eph['earth']

        e = earth.at(t)
        s = e.observe(sun).apparent()
        m = e.observe(moon).apparent()

        _, slon, _ = s.frame_latlon(ecliptic_frame)
        _, mlon, _ = m.frame_latlon(ecliptic_frame)
        phase = (mlon.degrees - slon.degrees) % 360.0

        percent = 100.0 * m.fraction_illuminated(sun)

    return phase, percent

def plot_airmass(RA, DEC, dt):
    """
    Plotting airmass for an object at a certain time
    
    params:
    RA, DEC - The coordinates of the target in degrees
    dt - The datetime of observation
    
    returns:
    fig - The plotly figure
    """
    # Define timezone and get UTC time of observations
    tz = timezone("Australia/Tasmania")
    aware1 = tz.localize(dt)
    utcoffset = aware1.utcoffset()
    time_utc = Time(dt) - utcoffset
    time_delta = np.arange(-6, 6, 0.1)*u.hour
    
    # Create sky coordinates for this object
    coord = SkyCoord(ra=RA, dec=DEC, unit=(u.deg, u.deg))
    target = FixedTarget(coord=coord)
    obj_airmass = GHO.altaz(time_utc+time_delta, target).secz
    #obj_azimuth = GHO.altaz(time_utc+time_delta, target).az.value
    
    # Plot the figure 
    fig = px.line(x=time_delta, y=obj_airmass, labels={'x': 'Hours from Midnight on {} (Tas)'.format(str(dt.date())), 'y': 'Airmass [Sec(z)]'}, title='Airmass for (RA,DEC) = ({},{}) at GHO'.format(RA, DEC), template='plotly_dark')
    fig.update_yaxes(range=[3, 1])
    fig.update_xaxes(range=[-6,6])

    return fig

def calc_meridian_transit(RA, DEC, dt):
    """
    Calculate the meridian transit time of an object
    
    params:
    RA, DEC - The coordinates of the target in degrees
    dt - The datetime of observation
    """
    tz = timezone("Australia/Tasmania")
    aware1 = tz.localize(dt)
    utc_offset = aware1.utcoffset()
    coord = SkyCoord(ra=RA*u.deg, dec=DEC*u.deg)
    target = FixedTarget(coord=coord)
    transit_time = GHO.target_meridian_transit_time(dt, target, which="next")
    
    return transit_time + utc_offset

def calc_moon_sep(RA, DEC, dt):
    """
    Calculate the angular seperation between the moon and a target
    
    params:
    RA, DEC - The coordinates of the target in degrees
    dt - The datetime of observation
    
    return: 
    sep.value - The seperation in degrees
    """
    t = Time(str(dt.isoformat()), format = 'isot')
    moon = GHO.moon_altaz(t)
    coord = SkyCoord(ra=RA*u.deg, dec=DEC*u.deg)
    sep = moon.separation(coord)
    
    return sep.value

def find_constellation(RA, DEC):
    """
    Find the constellation that a given coordinate resides in
    
    params:
    RA, DEC - The coordinates of the target in degrees
    
    return:
    
    const - The constellation string
    """
    coord = SkyCoord(ra=RA*u.deg, dec=DEC*u.deg)
    const = get_constellation(coord)
    return const
    
def find_night(dt, horizon = -18):
    """
    Find when it will be dark enough to observe. I have chosen the end of
    astronomical twilight and start of astronomical morning.
    
    params:
    dt - The datetime of the observation
    horizon - The location of the sun below the horizon at ast. twilight/morning
    """
    tz = timezone("Australia/Tasmania")
    aware1 = tz.localize(dt)
    utc_offset = aware1.utcoffset()
    t = Time(str(dt.isoformat()), format = 'isot')
    start_t = GHO.twilight_evening_astronomical(t, which='next', n_grid_points=150)
    end_t = GHO.twilight_morning_astronomical(t, which='next', n_grid_points=150)

    return start_t+utc_offset, end_t+utc_offset

def get_jd(dt):
    """
    Get the julian day of the observation
    
    params:
    dt - The datetime of the observation
    
    return:
    jd - The julian day
    """
    t = Time(str(dt.isoformat()), format = 'isot')
    jd = t.jd
    return jd
    
def get_npix(seeing):
    """
    A function to retrieve the number of pixels given seeing selection
    
    params:
    seeing - The seeing string 
    
    returns:
    n_pix - The number of pixels for SNR calculation (2*IM_FWHM)
    """
    if seeing == 'Average (2")':
        n_pix = 2*(2/0.8)
        
    elif seeing == 'Poor (2.5")':
        n_pix = 2*(2.5/0.8)
        
    elif seeing == 'Good (1.5")':
        n_pix = 2*(1.5/0.8)
        
    return n_pix