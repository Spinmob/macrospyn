import ctypes as _c
import numpy  as _n
import os     as _os
from sys import platform as _platform


# Find the path to the compiled c-code (only Windows and Linux supported so far.)
if _platform in ['win32']: _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine.dll')
else:                      _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine.so')

# Get the engine.
_engine = _c.cdll.LoadLibrary(_path_dll)

def to_pointer_double(numpy_array): 
    """
    Converts the supplied numpy_array (assumed to be the usual 64-bit float)
    to a pointer, allowing it ot be "connected" to the C-code engine.
    
    SUPER IMPORTANT
    ---------------
    Make sure you assign your array to a variable name prior to calling this.
    If you do not keep a handle on the numpy array, garbage collection can
    delete it while the simulation is running!
    
    Parameters
    ----------
    numpy_array
        1D ndarray of 64-bit floats.
        
    Returns
    -------
    C-pointer to the first element.
    
    Examples
    --------
    my_solver.Bxs = to_pointer_double(my_Bxs).
    
    """
    return numpy_array.ctypes.data_as(_c.POINTER(_c.c_double))



class _domain(_c.Structure):
    """
    Structure for sending all the simulation parameters to the c library.
    The user should not interact with this for the most part.
    """
    
    # NOTE: The order here and structure here REALLY has to match the struct 
    # in the C-code! This class will be sent by reference, and the c-code 
    # will expect everything to be in its place.
    _fields_ = [
        
        # Magnitude of the gyromagnetic ratio [radians / (sec T)]
        ('gamma', _c.c_double), ("_gammas", _c.POINTER(_c.c_double)),
        
        #Gilbert damping parameter [unitless]
        ('M', _c.c_double), ("_Ms", _c.POINTER(_c.c_double)),
        
        # Gilbert damping
        ('alpha', _c.c_double), ('_alphas', _c.POINTER(_c.c_double)),
        
        # Exchange-like field strength [T], applied in the direction of the other domain's unit vector
        ('X', _c.c_double), ('_Xs', _c.POINTER(_c.c_double)),
        
        # Spin transfer torque (rate) parallel to other domain [rad / s]
        ('s', _c.c_double), ('_ss', _c.POINTER(_c.c_double)),
        
        # Other torque (rate) unrelated to either domain [rad / s]
        ('tx', _c.c_double), ('_txs', _c.POINTER(_c.c_double)),
        ('ty', _c.c_double), ('_tys', _c.POINTER(_c.c_double)),
        ('tz', _c.c_double), ('_tzs', _c.POINTER(_c.c_double)),
        
        # Externally applied field
        ('Bx', _c.c_double), ('_Bxs', _c.POINTER(_c.c_double)),
        ('By', _c.c_double), ('_Bys', _c.POINTER(_c.c_double)),
        ('Bz', _c.c_double), ('_Bzs', _c.POINTER(_c.c_double)),
        
        # Anisotropy tensor elements [unitless], defined such that Nxx+Nyy+Nzz=1 for an aligned ellipsoid
        ('Nxx', _c.c_double), ('_Nxxs', _c.POINTER(_c.c_double)),
        ('Nxy', _c.c_double), ('_Nxys', _c.POINTER(_c.c_double)),
        ('Nxz', _c.c_double), ('_Nxzs', _c.POINTER(_c.c_double)),
    
        ('Nyx', _c.c_double), ('_Nyxs', _c.POINTER(_c.c_double)),
        ('Nyy', _c.c_double), ('_Nyys', _c.POINTER(_c.c_double)),
        ('Nyz', _c.c_double), ('_Nyzs', _c.POINTER(_c.c_double)),
        
        ('Nzx', _c.c_double), ('_Nzxs', _c.POINTER(_c.c_double)),
        ('Nzy', _c.c_double), ('_Nzys', _c.POINTER(_c.c_double)),
        ('Nzz', _c.c_double), ('_Nzzs', _c.POINTER(_c.c_double)),
    
        # Dipole tensor [unitless], representing the fraction of the other layer's saturation magnetization
        ('Dxx', _c.c_double), ('_Dxxs', _c.POINTER(_c.c_double)),
        ('Dxy', _c.c_double), ('_Dxys', _c.POINTER(_c.c_double)),
        ('Dxz', _c.c_double), ('_Dxzs', _c.POINTER(_c.c_double)),
        
        ('Dyx', _c.c_double), ('_Dyxs', _c.POINTER(_c.c_double)),
        ('Dyy', _c.c_double), ('_Dyys', _c.POINTER(_c.c_double)),
        ('Dyz', _c.c_double), ('_Dyzs', _c.POINTER(_c.c_double)),
        
        ('Dzx', _c.c_double), ('_Dzxs', _c.POINTER(_c.c_double)),
        ('Dzy', _c.c_double), ('_Dzys', _c.POINTER(_c.c_double)),
        ('Dzz', _c.c_double), ('_Dzzs', _c.POINTER(_c.c_double)),
        
        # Solver mode: 0=disabled, 1=LLG
        ('mode', _c.c_int),
        
        # Initial conditions
        ('x0',   _c.c_double),
        ('y0',   _c.c_double),
        ('z0',   _c.c_double),
        
        # Solution arrays
        ('_x', _c.POINTER(_c.c_double)),
        ('_y', _c.POINTER(_c.c_double)),
        ('_z', _c.POINTER(_c.c_double)),
        
    ] # End of _data structure
    
    def set(self, **kwargs): 
        """
        Sets any number of parameters for this domain.

        Parameters
        ----------
        **kwargs : keyword will be set by evaluating self.keyword = value. In 
        the case of an array, it will connect the appropriate (64-bit float) 
        pointer. This will also save a "local" copy of the supplied value, such
        that garbage collection doesn't automatically delete it.
        
        Example
        -------
        my_solver.a.set(Bys=numpy.linspace(0,5,my_solver.steps), gamma=27)

        Returns
        -------
        self

        """
        
        for k in kwargs:

            # Store the "local" copy, to make sure we keep it from getting
            # deleted by the garbage collection.
            exec('self.'+k+"=v", dict(self=self,v=kwargs[k]))
            
            # If it's an array, convert it before setting
            if type(kwargs[k])==_n.ndarray:
                exec('self._'+k+"=v", dict(self=self,v=to_pointer_double(kwargs[k])))
            
            # Otherwise, if it's None, update the pointer
            elif kwargs[k] is None:
                exec('self._'+k+"=v", dict(self=self,v=None))
            
        return self
    
    __call__ = set



class solver():
    """
    Class for interfacing with the compiled C-code engine.
    """
    
    def __init__(self, dt=0.005, steps=1e5):

        # Store the run parameters
        self.dt    = dt
        self.steps = steps
        
        # Create the settings structure, and set the default values.
        self.a = _domain()
        self.b = _domain()
        
        # Initial conditions
        self.a.x0   = 1.0
        self.a.y0   = 0.0
        self.a.z0   = 0.0
        self.b.x0   = 1.0
        self.b.y0   = 0.0
        self.b.z0   = 0.0

        # Null all the array pointers, just to be safe 
        # (different platforms, Python versions, etc...)
        self.a._gammas = self.b._gammas = None
        self.a._Ms     = self.b._Ms     = None
        self.a._alphas = self.b._alphas = None
        self.a._Xs     = self.b._Xs     = None
        self.a._ss     = self.b._ss     = None
        
        self.a._Bxs    = self.b._Bxs    = None
        self.a._Bys    = self.b._Bys    = None
        self.a._Bzs    = self.b._Bzs    = None
        
        self.a._txs    = self.b._txs    = None
        self.a._tys    = self.b._tys    = None
        self.a._tzs    = self.b._tzs    = None
        
        self.a._Nxxs   = self.b._Nxxs   = None
        self.a._Nxys   = self.b._Nxys   = None
        self.a._Nxzs   = self.b._Nxzs   = None
        self.a._Nyxs   = self.b._Nyxs   = None
        self.a._Nyys   = self.b._Nyys   = None
        self.a._Nyzs   = self.b._Nyzs   = None
        self.a._Nzxs   = self.b._Nzxs   = None
        self.a._Nzys   = self.b._Nzys   = None
        self.a._Nzzs   = self.b._Nzzs   = None
        
        self.a._Dxxs   = self.b._Dxxs   = None
        self.a._Dxys   = self.b._Dxys   = None
        self.a._Dxzs   = self.b._Dxzs   = None
        self.a._Dyxs   = self.b._Dyxs   = None
        self.a._Dyys   = self.b._Dyys   = None
        self.a._Dyzs   = self.b._Dyzs   = None
        self.a._Dzxs   = self.b._Dzxs   = None
        self.a._Dzys   = self.b._Dzys   = None
        self.a._Dzzs   = self.b._Dzzs   = None
        
        # By default, enable the two domains.
        self.a.mode=1
        self.b.mode=1
        
    def set(self, **kwargs): 
        """
        Sets any number of parameters for the solver. Magnetic parameters will
        be applied to both domains.

        Parameters
        ----------
        **kwargs : keyword will be set by evaluating self.keyword = value. In 
        the case of an array, it will connect the appropriate (64-bit float) 
        pointer. This will also save a "local" copy of the supplied value, such
        that garbage collection doesn't automatically delete it.
        
        Example
        -------
        my_solver.set(dt=0.005)

        Returns
        -------
        self

        """
        
        for k in kwargs:

            # If it's a property of the solver, store it in the solver
            if k in ['dt', 'steps', 'ax0', 'ay0', 'az0', 'bx0', 'by0', 'bz0']:
                exec('self.'+k+"=v", dict(self=self,v=kwargs[k]))
        
            # Otherwise, send it to a and b.
            else:
                self.a.set(**{k:kwargs[k]})
                self.b.set(**{k:kwargs[k]})
            
        return self
    
    __call__ = set
    
    def run(self, update_initial_condition=True):
        """
        Creates the solution arrays and runs the solver to fill them up.
        Afterward, the initial conditions are (by default) set to the 
        last value of the solution arrays.
        
        Parameters
        ----------
        update_initial_condition=True
            If True, the initial conditions (self.ax0, self.bx0, ...) will be
            set to the last value of the solution arrays.
        """
        self.steps = int(self.steps)
        
        # Create the solution arrays
        self.ax = _n.zeros(self.steps); self.ax[0] = self.a.x0
        self.ay = _n.zeros(self.steps); self.ay[0] = self.a.y0
        self.az = _n.zeros(self.steps); self.az[0] = self.a.z0
        self.bx = _n.zeros(self.steps); self.bx[0] = self.b.x0
        self.by = _n.zeros(self.steps); self.by[0] = self.b.y0
        self.bz = _n.zeros(self.steps); self.bz[0] = self.b.z0
        
        # Provide the c-code with access to the array data
        self.a._x = to_pointer_double(self.ax)
        self.a._y = to_pointer_double(self.ay)
        self.a._z = to_pointer_double(self.az)
        self.b._x = to_pointer_double(self.bx)
        self.b._y = to_pointer_double(self.by)
        self.b._z = to_pointer_double(self.bz)
        
        # Solve it.
        _engine.solve_heun.restype = None
        _engine.solve_heun(_c.byref(self.a), _c.byref(self.b), 
                           _c.c_double(self.dt), _c.c_int(self.steps))
        
        # Set the initial conditions for the next run
        if update_initial_condition:
            self.a.x0 = self.ax[-1]
            self.a.y0 = self.ay[-1]
            self.a.z0 = self.az[-1]
            self.b.x0 = self.bx[-1]
            self.b.y0 = self.by[-1]
            self.b.z0 = self.bz[-1]
        
        return self

    def plot(self, n1=0, n2=None, t0=0, **kwargs):
        """
        Creates a plot for inspecting trajectories.
        

        Parameters
        ----------
        n1=0, n2=None
            Start and stop indices respectively.
        
        **kwargs
            These are sent to spinmob.plot.xy.data().
        
        Returns
        -------
        self

        """
        # Find the upper index if needed
        if n2 is None: n2 = self.steps

        # Make the plot
        import spinmob as sm
        ts = m.dt*_n.array(range(len(m.ax))) + t0
        sm.plot.xy.data(ts, [m.ax[n1:n2],m.ay[n1:n2],m.az[n1:n2],
                               m.bx[n1:n2],m.by[n1:n2],m.bz[n1:n2]],
                        label=['ax','ay','az','bx','by','bz'], 
                        xlabel='Time', ylabel='Projection', **kwargs)
    
        return self

if __name__ == '__main__':
    
    # Create a solver instance
    m = solver(dt=0.005, steps=870)
    
    # Set up the physical parameters
    m.set(By=10.0, gamma=_n.pi*2, dt=0.005, zzz=300, steps=270, alpha=0.1/_n.pi/2)
    m.a.set(Bys=_n.linspace(0,2,m.steps))
    m.set(bz0=2)
    
    # Run it & plot.
    m.run().plot()
    m.run().plot(clear=0, t0=(m.steps-1)*m.dt)
    m.run().plot(clear=0, t0=(m.steps-1)*m.dt*2)
    
    
    
    
    
