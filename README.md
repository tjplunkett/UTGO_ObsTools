# Greenhill Observatory ObsTools

The UTAS Greenhill Observatory (UTGO) ‘ObsTools’ dashboard is a tool to aid in planning and executing
astronomical observations with the UTAS optical telescopes. It is written solely in Python, utilising the
‘Dash’ package to implement an interactive dashboard on your local machine. The main calculations
use the Astropy, Astroplan and Skyfields packages. The dashboard contains a basic sky map, table
with astronomical information and unit conversions, airmass plotter and 3 variations of exposure time
calculation plots. 

### Installation

It is recommended that the user downloads and installs Anaconda or Miniconda to manage their
python distributions. First, download/clone the Github repository to your local machine: ObsTools
Github. You will need to change into the downloaded directory, using:

#### cd [insert path to folder here]

Then, you will need to install the necessary requirements to run the dashboard. To do this, go to
your terminal (or Anaconda Prompt on Windows) and create a new conda environment (for example
‘obs’):

#### conda create -n obs python=3.9.16

You will then need to activate the environment:

#### conda activate obs

Then, assuming that ‘pip’ installed into this environment during the process, type:

#### pip install -r requirements.txt

### How to Use

To run the dashboard, activate
your conda environment, change to the ObsTool directory and type:
python obstools.py
This will launch the dashboard on your local host. Go to a web browser and type the address given in
the command line.

For more detailed help, please see the ObsTools_UserGuide in the Documentation folder.
