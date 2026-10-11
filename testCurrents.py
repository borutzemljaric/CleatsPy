import numpy as np
from i_magnetic.i_currents import i_currents
from prestatic.create_classes import cable_class



cab=cable_class
t=0
t = np.linspace(0, 1, 10)   # 3 points

Nel=6
cab.tfinal =1 # (s) whole length of simulation
cab.tsim = 0.01  #  (s) delay of simulation from start for better plots
cab.trelay=0.5 # (s) lengt of short circuit - time required to  switch off
cab.interval = .1 # (s) time length of sub-intervals, for easier memory management, must be smaller than tfinal
cab.ste_delitev = int(np.round(cab.tfinal / cab.interval))
cab.note = 0.02  # (s ) interval of time for saved results
cab.Ik=1000  # (A) electric current Ik''
cab.phi=20 # (degrees)surge impedance angle
cab.Nseg=2 # number of sub-segment, minimum 2

for time in t :
    res=i_currents(time, cab, Nel)
    print (res)