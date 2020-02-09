import ctypes as _c
import numpy  as _n
import os     as _os
from sys import platform as _platform


# Find the path to the compiled c-code (only Windows and Linux supported so far.)
if _platform in ['win32']: _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine.dll')
else:                      _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine.so')

# Get the engine.
_engine = _c.CDLL(_path_dll)

class _domain(_c.Structure):
    """
    Structure for sending all the simulation parameters to the c library.
    The user should not interact with this for the most part.
    """
    
    # NOTE: This thing REALLY has to match the struct in the C-code.
    _fields_ = [
        
        # Magnitude of the gyromagnetic ratio [radians / (sec T)]
        ('gamma', _c.c_double), ("gammas", _c.POINTER(_c.c_double)),
        
        #Gilbert damping parameter [unitless]
        ('M', _c.c_double), ("Ms", _c.POINTER(_c.c_double)),
        
        # Gilbert damping
        ('alpha', _c.c_double), ('alphas', _c.POINTER(_c.c_double)),
        
        # Exchange-like field strength [T], applied in the direction of the other domain's unit vector
        ('X', _c.c_double), ("Xs", _c.POINTER(_c.c_double)),
        
        # Spin transfer torque (rate) parallel to other domain [rad / s]
        ('s', _c.c_double), ("ss", _c.POINTER(_c.c_double)),
        
        # Other torque (rate) unrelated to either domain [rad / s]
        ('tx', _c.c_double), ("txs", _c.POINTER(_c.c_double)),
        ('ty', _c.c_double), ("tys", _c.POINTER(_c.c_double)),
        ('tz', _c.c_double), ("tzs", _c.POINTER(_c.c_double)),
        
        # Externally applied field
        ('Bx', _c.c_double), ("Bxs", _c.POINTER(_c.c_double)),
        ('By', _c.c_double), ("Bys", _c.POINTER(_c.c_double)),
        ('Bz', _c.c_double), ("Bzs", _c.POINTER(_c.c_double)),
        
        # Anisotropy tensor elements [unitless], defined such that Nxx+Nyy+Nzz=1 for an aligned ellipsoid
        ('Nxx', _c.c_double), ("Nxxs", _c.POINTER(_c.c_double)),
        ('Nxy', _c.c_double), ("Nxys", _c.POINTER(_c.c_double)),
        ('Nxz', _c.c_double), ("Nxzs", _c.POINTER(_c.c_double)),
        
        ('Nyx', _c.c_double), ("Nyxs", _c.POINTER(_c.c_double)),
        ('Nyy', _c.c_double), ("Nyys", _c.POINTER(_c.c_double)),
        ('Nyz', _c.c_double), ("Nyzs", _c.POINTER(_c.c_double)),
        
        ('Nzx', _c.c_double), ("Nzxs", _c.POINTER(_c.c_double)),
        ('Nzy', _c.c_double), ("Nzys", _c.POINTER(_c.c_double)),
        ('Nzz', _c.c_double), ("Nzzs", _c.POINTER(_c.c_double)),
        
        # Dipole tensor [unitless], representing the fraction of the other layer's saturation magnetization
        ('Dxx', _c.c_double), ("Dxxs", _c.POINTER(_c.c_double)),
        ('Dxy', _c.c_double), ("Dxys", _c.POINTER(_c.c_double)),
        ('Dxz', _c.c_double), ("Dxzs", _c.POINTER(_c.c_double)),
        
        ('Dyx', _c.c_double), ("Dyxs", _c.POINTER(_c.c_double)),
        ('Dyy', _c.c_double), ("Dyys", _c.POINTER(_c.c_double)),
        ('Dyz', _c.c_double), ("Dyzs", _c.POINTER(_c.c_double)),
        
        ('Dzx', _c.c_double), ("Dzxs", _c.POINTER(_c.c_double)),
        ('Dzy', _c.c_double), ("Dzys", _c.POINTER(_c.c_double)),
        ('Dzz', _c.c_double), ("Dzzs", _c.POINTER(_c.c_double)),
        
        # Solver mode: 0=disabled, 1=LLG
        ('mode', _c.c_int),
        
        # Solution arrays
        ("x", _c.POINTER(_c.c_double)),
        ("y", _c.POINTER(_c.c_double)),
        ("z", _c.POINTER(_c.c_double)),
        
    ] # End of _data structure
    

class solver():
    """
    Class for interfacing with the macrospin simulation c-library.
    """
    
    def __init__(self, dt=0.005, steps=1e5):

        # Store the run parameters
        self.dt    = dt
        self.steps = steps
        
        # Initial conditions
        self.ax0   = 1.0
        self.ay0   = 0.0
        self.az0   = 0.0
        self.bx0   = 1.0
        self.by0   = 0.0
        self.bz0   = 0.0

        # Create the settings structure, and set the default values.
        self.a = _domain()
        self.b = _domain()
        
        # Null all the array pointers, just to be safe 
        # (different platforms, Python versions, etc...)
        self.a.gammas = self.b.gammas = None
        self.a.Ms     = self.b.Ms     = None
        self.a.alphas = self.b.alphas = None
        self.a.Xs     = self.b.Xs     = None
        self.a.ss     = self.b.ss     = None
        
        self.a.Bxs    = self.b.Bxs    = None
        self.a.Bys    = self.b.Bys    = None
        self.a.Bzs    = self.b.Bzs    = None
        
        self.a.txs    = self.b.txs    = None
        self.a.tys    = self.b.tys    = None
        self.a.tzs    = self.b.tzs    = None
        
        self.a.Nxxs   = self.b.Nxxs   = None
        self.a.Nxys   = self.b.Nxys   = None
        self.a.Nxzs   = self.b.Nxzs   = None
        self.a.Nyxs   = self.b.Nyxs   = None
        self.a.Nyys   = self.b.Nyys   = None
        self.a.Nyzs   = self.b.Nyzs   = None
        self.a.Nzxs   = self.b.Nzxs   = None
        self.a.Nzys   = self.b.Nzys   = None
        self.a.Nzzs   = self.b.Nzzs   = None
        
        self.a.Dxxs   = self.b.Dxxs   = None
        self.a.Dxys   = self.b.Dxys   = None
        self.a.Dxzs   = self.b.Dxzs   = None
        self.a.Dyxs   = self.b.Dyxs   = None
        self.a.Dyys   = self.b.Dyys   = None
        self.a.Dyzs   = self.b.Dyzs   = None
        self.a.Dzxs   = self.b.Dzxs   = None
        self.a.Dzys   = self.b.Dzys   = None
        self.a.Dzzs   = self.b.Dzzs   = None
        
        
        # By default, enable the two domains.
        self.a.mode=1
        self.b.mode=1
        
    def __repr__(self):
        """
        Returns the string that appears when you inspect the object.
        """
        s = "solver instance"
        
        return s
    
    def run(self):
        """
        Creates the solution arrays and runs the solver.
        """
        self.steps = int(self.steps)
        
        # Create the solution arrays
        self.ax = _n.zeros(self.steps); self.ax[0] = self.ax0
        self.ay = _n.zeros(self.steps); self.ay[0] = self.ay0
        self.az = _n.zeros(self.steps); self.az[0] = self.az0
        self.bx = _n.zeros(self.steps); self.bx[0] = self.bx0
        self.by = _n.zeros(self.steps); self.by[0] = self.by0
        self.bz = _n.zeros(self.steps); self.bz[0] = self.bz0
        
        # Provide the c-code with access to the numpy array data
        self.a.x = _n.ctypeslib.as_ctypes(self.ax)
        self.a.y = _n.ctypeslib.as_ctypes(self.ay)
        self.a.z = _n.ctypeslib.as_ctypes(self.az)
        self.b.x = _n.ctypeslib.as_ctypes(self.bx)
        self.b.y = _n.ctypeslib.as_ctypes(self.by)
        self.b.z = _n.ctypeslib.as_ctypes(self.bz)
        
        # Solve it.
        _engine.solve_heun.restype = None
        _engine.solve_heun(_c.byref(self.a), 
                           _c.byref(self.b), 
                           _c.c_double(self.dt), 
                           _c.c_int(self.steps))
        
        return self




m = solver()
m.steps=1000

m.a.gamma=2*_n.pi; m.a.By=1.0
m.b.gamma=2*_n.pi; m.b.By=1.0

m.run()

import spinmob as sm
sm.plot.xy.data(m.dt, [m.ax,m.ay,m.az,m.bx,m.by,m.bz],
                label=['ax','ay','az','bx','by','bz'])
    
    
    
