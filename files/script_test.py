import math
import pandas as pd
import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import random
import pickle

# Define the reference start date
REFERENCE_DATE = datetime(1995, 1, 1)
PLOT_CUTTED = False


print("reading a file")
with open(f"files/tle_plots/tle_{1}_with_maneuvers{'_cut' if False else ''}.pkl", 'rb') as f:
    fig = pickle.load(f)
    plt.show()
print("done")