import ctypes  as _c
import numpy   as _n
import os      as _os
import spinmob as _s
import spinmob.egg as _egg; _g = _egg.gui
import pyqtgraph.opengl as _gl
import pyqtgraph        as _pg
from sys import platform as _platform
import traceback as _t
_p = _t.print_last

# Find the path to the compiled c-code (only Windows and Linux supported so far.)
if _platform in ['win32']: _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine.dll')
else:                      _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine.so')

# Get the engine.
_engine = _c.cdll.LoadLibrary(_path_dll)


## TO DO: Get check boxes working!

def _to_pointer(numpy_array): 
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
    my_solver.Bxs = _to_pointer(my_Bxs).
    
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
        ('STT', _c.c_double), ('_STTs', _c.POINTER(_c.c_double)),
        
        # Other torque (rate) unrelated to either domain [rad / s]
        ('Tx', _c.c_double), ('_Txs', _c.POINTER(_c.c_double)),
        ('Ty', _c.c_double), ('_Tys', _c.POINTER(_c.c_double)),
        ('Tz', _c.c_double), ('_Tzs', _c.POINTER(_c.c_double)),
        
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
        the case of an array, it will convert it to a pointer and (64-bit float)
        use an underscore, as needed, saving a "local" copy of the supplied 
        value, such that garbage collection doesn't automatically delete it.
        
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
            
            # If it's an array, convert it before setting and use _
            if type(kwargs[k])==_n.ndarray:
                exec('self._'+k+"=v", dict(self=self,v=_to_pointer(kwargs[k])))
            
            # Otherwise, if it's None, update the pointer
            elif kwargs[k] is None:
                exec('self._'+k+"=v", dict(self=self,v=None))
            
        return self
    
    __call__ = set



class solver_api():
    """
    Scripted interface for the solver engine. Keyword arguments are 
    sent to self.set().
    """
    
    def __init__(self, **kwargs):

        # Store the run parameters
        self.dt    = 1e-12
        self.steps = 1e3
        
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
        self.a._STTs   = self.b._STTs   = None
        
        self.a._Bxs    = self.b._Bxs    = None
        self.a._Bys    = self.b._Bys    = None
        self.a._Bzs    = self.b._Bzs    = None
        
        self.a._Txs    = self.b._Txs    = None
        self.a._Tys    = self.b._Tys    = None
        self.a._Tzs    = self.b._Tzs    = None
        
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
        
        # No arrays initially
        self.ax = self.ay = self.az = None
        self.bx = self.by = self.bz = None
        
        self.set(**kwargs)
        
    def set(self, **kwargs): 
        """
        Sets any number of parameters for the solver. Magnetic parameters will
        be applied to both domains.

        Parameters
        ----------
        **kwargs : keyword will be set by evaluating self.keyword = value or
        self.a.keyword = value and self.b.keyword=value, as appropriate. 
        
        In the case of an array value, it will convert to pointer (64-bit float) 
        pointer and evaluate self.a._keyword=pointer. This will also save a 
        "local" copy of the supplied value, such that garbage collection 
        doesn't automatically delete it.
        
        Example
        -------
        my_solver.set(dt=0.002)

        Returns
        -------
        self

        """
        
        for k in kwargs:

            # If it's a property of the solver, store it in the solver
            if k in ['dt', 'steps']:
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
            If True, the initial conditions (self.a.x0, self.b.x0, ...) will be
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
        self.a._x = _to_pointer(self.ax)
        self.a._y = _to_pointer(self.ay)
        self.a._z = _to_pointer(self.az)
        self.b._x = _to_pointer(self.bx)
        self.b._y = _to_pointer(self.by)
        self.b._z = _to_pointer(self.bz)
        
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


class solver():
    """
    Graphical and scripted interface for the solver engine.
    
    Keyword arguments are sent to self.set()
    """    
    
    def __init__(self):
        
        # Solver application programming interface.
        self.api = solver_api()
        
        # Graphical interface
        self.window = _g.Window(title='Macrospin(mob)', autosettings_path='solver.window.txt', size=[1000,550])
        
        # Top row controls for the "go" button, etc
        self.grid_top          = self.window  .place_object(_g.GridLayout(False), alignment=1) 
        self.button_go         = self.grid_top.place_object(_g.Button('Go!'))
        self.number_iterations = self.grid_top.place_object(_g.NumberBox(1, bounds=(0,None), int=True))
        self.label_iteration   = self.grid_top.place_object(_g.Label(''))
        
        # Bottom row controls for settings and plots.
        self.window.new_autorow()
        self.grid_bottom  = self.window     .place_object(_g.GridLayout(False), alignment=0)
        
        # Settings
        self.settings     = self.grid_bottom.place_object(_g.TreeDictionary(autosettings_path='solver.settings.txt'))
        
        self.settings.add_parameter('solver/dt',    1e-12, dec=True,                  siPrefix=True, suffix='s')
        self.settings.add_parameter('solver/steps', 5000,  dec=True, limits=(2,None), siPrefix=True, suffix='steps')
        self.settings.add_parameter('solver/reset', True)
        
        self.settings.add_parameter('a/mode', 1, limits=(0,1), tip='0=disabled, 1=LLG')
        
        self.settings.add_parameter('a/initial_condition/x0', 1.0, tip='Initial magnetization direction (will be normalized to unit length)')
        self.settings.add_parameter('a/initial_condition/y0', 1.0, tip='Initial magnetization direction (will be normalized to unit length)')
        self.settings.add_parameter('a/initial_condition/z0', 0.0, tip='Initial magnetization direction (will be normalized to unit length)')
        
        self.settings.add_parameter('a/material/gamma', 1.760859644e11, siPrefix=True, suffix='rad/(s*T)', tip='Magnitude of gyromagnetic ratio')
        self.settings.add_parameter('a/material/M',     1.0, siPrefix=True, suffix='T', tip='Saturation magnetization (u0*Ms)')
        self.settings.add_parameter('a/material/alpha', 0.01, tip='Gilbert damping parameter')        
        
        self.settings.add_parameter('a/applied_field', True)
        self.settings.add_parameter('a/applied_field/Bx', 0.0, siPrefix=True, suffix='T', tip='Externally applied magnetic field')
        self.settings.add_parameter('a/applied_field/By', 0.0, siPrefix=True, suffix='T', tip='Externally applied magnetic field')
        self.settings.add_parameter('a/applied_field/Bz', 0.0, siPrefix=True, suffix='T', tip='Externally applied magnetic field')
        
        self.settings.add_parameter('a/other_torques', True)
        self.settings.add_parameter('a/other_torques/X', 0.0, siPrefix=True, suffix='T', tip='Exchange field parallel to domain b\'s magnetization')
        self.settings.add_parameter('a/other_torques/STT', 0.0, siPrefix=True, suffix='rad/s', tip='Spin-transfer-like torque, parallel to domain b\'s magnetization')
        
        self.settings.add_parameter('a/other_torques/Tx', 0.0, siPrefix=True, suffix='rad/s', tip='Other externally applied torque')
        self.settings.add_parameter('a/other_torques/Ty', 0.0, siPrefix=True, suffix='rad/s', tip='Other externally applied torque')
        self.settings.add_parameter('a/other_torques/Tz', 0.0, siPrefix=True, suffix='rad/s', tip='Other externally applied torque')
        
        self.settings.add_parameter('a/anisotropy', True)
        self.settings.add_parameter('a/anisotropy/Nxx', 0.01, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('a/anisotropy/Nyy', 0.10, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('a/anisotropy/Nzz', 0.89, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('a/anisotropy/Nxy', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('a/anisotropy/Nxz', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('a/anisotropy/Nyx', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('a/anisotropy/Nyz', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('a/anisotropy/Nzx', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('a/anisotropy/Nzy', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        
        self.settings.add_parameter('a/dipole', True)
        self.settings.add_parameter('a/dipole/Dxx', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('a/dipole/Dyy', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('a/dipole/Dzz', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('a/dipole/Dxy', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('a/dipole/Dxz', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('a/dipole/Dyx', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('a/dipole/Dyz', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('a/dipole/Dzx', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('a/dipole/Dzy', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        
        self.settings.add_parameter('b/mode', 1, limits=(0,1), tip='0=disabled, 1=LLG')
        
        self.settings.add_parameter('b/initial_condition/x0',-1.0, tip='Initial magnetization direction (will be normalized to unit length)')
        self.settings.add_parameter('b/initial_condition/y0', 1.0, tip='Initial magnetization direction (will be normalized to unit length)')
        self.settings.add_parameter('b/initial_condition/z0', 0.0, tip='Initial magnetization direction (will be normalized to unit length)')
        
        self.settings.add_parameter('b/material/gamma', 1.760859644e11, siPrefix=True, suffix='rad/(s*T)', tip='Magnitude of gyromagnetic ratio')
        self.settings.add_parameter('b/material/M',     1.0, siPrefix=True, suffix='T', tip='Saturation magnetization (u0*Ms)')
        self.settings.add_parameter('b/material/alpha', 0.01, tip='Gilbert damping parameter')        
        
        self.settings.add_parameter('b/applied_field', True)
        self.settings.add_parameter('b/applied_field/Bx', 0.0, siPrefix=True, suffix='T', tip='Externally applied magnetic field')
        self.settings.add_parameter('b/applied_field/By', 0.0, siPrefix=True, suffix='T', tip='Externally applied magnetic field')
        self.settings.add_parameter('b/applied_field/Bz', 0.0, siPrefix=True, suffix='T', tip='Externally applied magnetic field')
        
        self.settings.add_parameter('b/other_torques', True)
        self.settings.add_parameter('b/other_torques/X', 0.0, siPrefix=True, suffix='T', tip='Exchange field parallel to domain b\'s magnetization')
        self.settings.add_parameter('b/other_torques/STT', 0.0, siPrefix=True, suffix='rad/s', tip='Spin-transfer-like torque, parallel to domain b\'s magnetization')
        
        self.settings.add_parameter('b/other_torques/Tx', 0.0, siPrefix=True, suffix='rad/s', tip='Other externally applied torque')
        self.settings.add_parameter('b/other_torques/Ty', 0.0, siPrefix=True, suffix='rad/s', tip='Other externally applied torque')
        self.settings.add_parameter('b/other_torques/Tz', 0.0, siPrefix=True, suffix='rad/s', tip='Other externally applied torque')
        
        self.settings.add_parameter('b/anisotropy', True)
        self.settings.add_parameter('b/anisotropy/Nxx', 0.01, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('b/anisotropy/Nyy', 0.10, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('b/anisotropy/Nzz', 0.89, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('b/anisotropy/Nxy', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('b/anisotropy/Nxz', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('b/anisotropy/Nyx', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('b/anisotropy/Nyz', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('b/anisotropy/Nzx', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        self.settings.add_parameter('b/anisotropy/Nzy', 0.0, tip='Anisotropy matrix (diagonal matrix has values adding to 1)')
        
        self.settings.add_parameter('b/dipole', True)
        self.settings.add_parameter('b/dipole/Dxx', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('b/dipole/Dyy', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('b/dipole/Dzz', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('b/dipole/Dxy', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('b/dipole/Dxz', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('b/dipole/Dyx', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('b/dipole/Dyz', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('b/dipole/Dzx', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        self.settings.add_parameter('b/dipole/Dzy', 0.0, tip='Dipolar field matrix exerted by domain b, expressed\nas a fraction of b\'s saturation magnetization.')
        
        # Plot tabs
        self.tabs         = self.grid_bottom.place_object(_g.TabArea(autosettings_path='solver.tabs.txt'), alignment=0)
        
        # Inspection plot for all arrays
        self.tab_inspect  = self.tabs.add_tab('Inspect')
        self.plot_inspect = self.tab_inspect.place_object(_g.DataboxPlot(autoscript=5, autosettings_path='solver.plot_inspect.txt'), alignment=0)
        self.plot_inspect['t']  = []
        self.plot_inspect['ax'] = []
        self.plot_inspect['ay'] = []
        self.plot_inspect['az'] = []
        self.plot_inspect['bx'] = []
        self.plot_inspect['by'] = []
        self.plot_inspect['bz'] = []
        
        # 3D plot
        self.tab_3d = self.tabs.add_tab('3D')
        self.button_3d_a   = self.tab_3d.place_object(_g.Button('a',   checkable=True, checked=True)) 
        self.button_3d_b   = self.tab_3d.place_object(_g.Button('b',   checkable=True, checked=True)) 
        self.button_3d_sum = self.tab_3d.place_object(_g.Button('Sum', checkable=True, checked=True)) 
        self.button_plot_3d = self.tab_3d.place_object(_g.Button('Update Plot'))
        self.tab_3d.new_autorow()
        
        # Make the 3D plot window
        self._widget_3d = _gl.GLViewWidget()
        self._widget_3d.opts['distance'] = 50
        
        # Make the grids
        self._gridx_3d = _gl.GLGridItem()
        self._gridx_3d.rotate(90,0,1,0)
        self._gridx_3d.translate(-10,0,0)
        self._widget_3d.addItem(self._gridx_3d)
        self._gridy_3d = _gl.GLGridItem()
        self._gridy_3d.rotate(90,1,0,0)
        self._gridy_3d.translate(0,-10,0)
        self._widget_3d.addItem(self._gridy_3d)
        self._gridz_3d = _gl.GLGridItem()
        self._gridz_3d.translate(0,0,-10)
        self._widget_3d.addItem(self._gridz_3d)
        
        # Trajectories
        color_a = _pg.glColor(100,100,255)
        color_b = _pg.glColor(255,100,100)
        color_n = _pg.glColor(50,255,255)
        self._trajectory_a_3d   = _gl.GLLinePlotItem(color=color_a, width=2.5, antialias=True)
        self._trajectory_b_3d   = _gl.GLLinePlotItem(color=color_b, width=2.5, antialias=True)
        self._trajectory_sum_3d = _gl.GLLinePlotItem(color=color_n, width=2.5, antialias=True)
        self._widget_3d.addItem(self._trajectory_a_3d)
        self._widget_3d.addItem(self._trajectory_b_3d)
        self._widget_3d.addItem(self._trajectory_sum_3d)
        
        # Other items
        self._start_dot_a_3d   = _gl.GLScatterPlotItem(color=color_a, size=7.0, pos=_n.array([[10,0,0]]))
        self._start_dot_b_3d   = _gl.GLScatterPlotItem(color=color_b, size=7.0, pos=_n.array([[-10,0,0]]))
        self._start_dot_sum_3d = _gl.GLScatterPlotItem(color=color_n, size=7.0, pos=_n.array([[-10,0,0]]))
        self._widget_3d.addItem(self._start_dot_a_3d)
        self._widget_3d.addItem(self._start_dot_b_3d)
        self._widget_3d.addItem(self._start_dot_sum_3d)
        self._update_start_dots()
        
        # Add the 3D plot window to the tab
        self.tab_3d.place_object(self._widget_3d, column_span=4, alignment=0)
        self.tab_3d.set_column_stretch(3)
        
        # Dump all the (autoloaded already) settings to the API
        for k in self.settings.keys():
            s = k.split('/')
            self._set_domain_parameter(s[0], s[-1], self.settings[k])
        
        # When one of the settings changes, make sure to update the solver_api
        self.settings.connect_any_signal_changed(self._setting_changed)
        
        # Connect the other controls
        self.button_go     .signal_clicked.connect(self.button_go_clicked)
        self.button_plot_3d.signal_clicked.connect(self.button_plot_3d_clicked)
        
        # Transfer the solver data based on the check boxes.
        self.settings.emit_signal_changed('a/applied_field')        
        self.settings.emit_signal_changed('a/other_torques')        
        self.settings.emit_signal_changed('a/anisotropy')        
        self.settings.emit_signal_changed('a/dipole')        
        self.settings.emit_signal_changed('b/applied_field')        
        self.settings.emit_signal_changed('b/other_torques')        
        self.settings.emit_signal_changed('b/anisotropy')        
        self.settings.emit_signal_changed('b/dipole')        
        
        # Let's have a look!
        self.window.show()

    def _update_start_dots(self):
        """
        Gets the initial condition from the settings and updates the start
        dot positions.
        """
        ax = self.settings['a/initial_condition/x0']
        ay = self.settings['a/initial_condition/y0']
        az = self.settings['a/initial_condition/z0']
        an = 1.0/_n.sqrt(ax*ax+ay*ay+az*az)
        ax = ax*an
        ay = ay*an
        az = az*an
        
        
        bx = self.settings['b/initial_condition/x0']
        by = self.settings['b/initial_condition/y0']
        bz = self.settings['b/initial_condition/z0']
        bn = 1.0/_n.sqrt(bx*bx+by*by+bz*bz)
        bx = bx*bn
        by = by*bn
        bz = bz*bn
        
        self._start_dot_a_3d  .setData(pos=10*_n.array([[ax, ay, az]]))
        self._start_dot_b_3d  .setData(pos=10*_n.array([[bx, by, bz]]))
        self._start_dot_sum_3d.setData(pos=10*_n.array([[ax+bx, ay+by, az+bz]]))

    def button_go_clicked(self, *a):
        """
        Go button pressed: Run the simulation!
        """

        for n in range(self.number_iterations.get_value()):
            self.label_iteration.set_text(str(n+1))
            self.go(self.settings['solver/reset'])
            self.button_plot_3d_clicked()
            
    def button_plot_3d_clicked(self, *a):
        """
        Plot 3d button pressed: Update the plot!
        """
        d = self.plot_inspect
        if self.button_3d_a.is_checked():   self._trajectory_a_3d  .setData(pos=10*_n.vstack([d['ax'],d['ay'],d['az']]).transpose())
        if self.button_3d_b.is_checked():   self._trajectory_b_3d  .setData(pos=10*_n.vstack([d['bx'],d['by'],d['bz']]).transpose())
        if self.button_3d_sum.is_checked(): self._trajectory_sum_3d.setData(pos=10*_n.vstack([d['ax']+d['bx'],d['ay']+d['by'],d['az']+d['bz']]).transpose())
        
        self.window.process_events()
        
        
    def go(self, reset=False):
        """
        Run the specified simulation.
        
        Parameters
        ----------
        reset=False
            After the run, update the initial conditions to match the 
            last point calculated.
        """
        
        # Run it.
        self.api.run(not self.settings['solver/reset'])
        
        # Transfer to the initial condition
        self.settings['a/initial_condition/x0'] = self.api.a.x0
        self.settings['a/initial_condition/y0'] = self.api.a.y0
        self.settings['a/initial_condition/z0'] = self.api.a.z0
        self.settings['b/initial_condition/x0'] = self.api.b.x0
        self.settings['b/initial_condition/y0'] = self.api.b.y0
        self.settings['b/initial_condition/z0'] = self.api.b.z0
        
        # Transfer the results to the inspector
        self.plot_inspect['t']  = self.api.dt*_n.array(range(self.api.steps))
        self.plot_inspect['ax'] = self.api.ax
        self.plot_inspect['ay'] = self.api.ay
        self.plot_inspect['az'] = self.api.az
        self.plot_inspect['bx'] = self.api.bx
        self.plot_inspect['by'] = self.api.by
        self.plot_inspect['bz'] = self.api.bz
        self.plot_inspect.plot()
        
        self.window.process_events()

    def _before_settings_set_value(self, key, value, ignore_error, block_user_signals):
        """
        Function called right before someone sets a settings value. If the
        value is an array, this will send it to the Inspect databox as a column
        of data with the specified key, and abort the "usual" set. Otherwise, it 
        will try to remove the column from the databox (if it exists) and
        set the value normally.
        """
        return True
        # s = key.split('/')
        
        # if type(value) == _n.ndarray:
            
        #     # Set it in the databox Inspector
        #     self.plot_inspect[key] = value
            
        #     # Send the array to the solver
        #     if   s[0] == 'a': self.api.a.set(**{s[-1]:value})
        #     elif s[0] == 'b': self.api.b.set(**{s[-1]:value})
                        
        #     return False # Don't set it as usual!
        
        # # Normal single value. Better set the arrays to None!
        # else:
        #     # Pop the column if it exists
        #     if s[-1] in self.plot_inspect.ckeys: self.plot_inspect.pop(s[-1])
            
        #     # Make sure the solver array value is None
        #     if   s[0] == 'a': 
            
        #     # Now do the thing.
        #     return True

    def _setting_changed(self, *a):
        """
        Called whenever settings are changed.
        """
        
        self.x = a
        
        # Get the associated item name and value
        domain    = a[0].name()
        parameter = a[1][0][0].name()
        value     = a[1][0][2]
        self._transfer_to_solver(domain, parameter, value)
    
    def _transfer_to_solver(self, domain, parameter, value):
        """
        Transfers the domain's parameter to the solver data.
        """
        # If we're enabling / disabling one of the settings, set all the 
        # sub-parameters
        if parameter in ['applied_field', 'other_torques', 'anisotropy', 'dipole']:
            
            # Assemble the root of the keys we wish to update
            root = domain+'/'+parameter
            
            # Loop over the keys
            for k in self.settings.keys():
                
                # Set the values if the root matches
                if k.find(root) == 0:
                    
                    # Get the associated item name and value
                    s = k.split('/')
                    domain    = s[0]
                    parameter = s[-1]
                    if value: value = self.settings[k]
                    else:     value = 0
                    
                    # Set it.
                    self._set_domain_parameter(domain, parameter, value)
                                
        # Otherwise it's a normal parameter.
        else: self._set_domain_parameter(domain, parameter, value)
        
        # Update the start dots
        self._update_start_dots()
        
    def _set_domain_parameter(self, domain, parameter, value):
        """
        Sets the value given domain and parameter.
        """
        # Come up with the command to execute
        if domain=='solver': command =            'self.api.set('+parameter+'='+str(value)+')'
        else:                command = 'self.api.'+domain+'.set('+parameter+'='+str(value)+')'   
        
        # Try it!
        try:    exec(command, dict(self=self))
        except: print('FAIL: "'+command+'"')
        

    def set(self, **kwargs):
        return

if __name__ == '__main__':
    
    # # Create a solver instance
    # m = solver_api(dt=0.005, steps=870)
    
    # # Set up the physical parameters
    # m.set(By=10.0, gamma=_n.pi*2, dt=0.005, zzz=300, steps=270, alpha=0.1/_n.pi/2)
    # m.a.set(Bys=_n.linspace(0,2,m.steps))
    # m.set(bz0=2)
    
    # # Run it & plot.
    # m.run().plot()
    # m.run().plot(clear=0, t0=(m.steps-1)*m.dt)
    # m.run().plot(clear=0, t0=(m.steps-1)*m.dt*2)
    
    self = solver()
    
    
    
    
    
    
