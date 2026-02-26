# reads output in csv from swat and makes 3d plots using plotly
# Feb 25 2026
import numpy as np
import math
import plotly.graph_objects as go
import sys
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

EARTH_R_KM = 6371.0
#####
# to do read csv.
