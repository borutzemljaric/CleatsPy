from dataclasses import dataclass

# function with declaration of classes used in dynamic part

@dataclass
class cable_class:
    name:str
    material:str
    A:float
    d:float
    m:float
    EI:float
    EA:float
    kn: float
    mi: float
    tfinal: float = None  # Simulation length  # this part is used in dynamic part
    tsim: float = None  # Time delay, real start of simulation
    trelay: float = None # duration of short circuit before relay protection reaction
    interval: int = None  # Splitting time into more intervals, for better efficiency
    ste_delitev: int  = None # Number of divisions, rounded (tfinal/interval)
    note: float  = None # time step to save in file
    Ik: float = None  # electric current
    phi: float = None # surge impedance angle
    Nseg: int = None # number of sub-segments for B calculation
    Ncon: int =None # number of contact points (checks) on element
    a: float=None  # damping coefficient alpha
    b: float=None  # damping coefficient beta


@dataclass
class plot_class:
    k:int # span index to plot (override python indexing, so 1 span is k=1)
    n_files:int  # number of files to be opened and compared
    n_curves:int # number of curves to plot
    plot_type:int #plot type






