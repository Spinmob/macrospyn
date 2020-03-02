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
if   _platform in ['win32']:  _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine-windows.dll')
elif _platform in ['darwin']: _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine-osx.so')
else:                         _path_dll = _os.path.join(_os.path.split(__file__)[0],'engine-linux.so')

if _platform == 'win32': _3d_enabled = False
else:                    _3d_enabled = True

# Get the engine.
_engine = _c.cdll.LoadLibrary(_path_dll)

# Constants
pi   = _n.pi
u0   = 1.25663706212e-6 # Vacuum permeability [H/m | N/A^2]
ec   = 1.60217662e-19   # Elementary charge [Coulomb]
me   = 9.1093837015e-31 # Electron mass [kg]
hbar = 1.0545718e-34    # [m^2 kg/s | J s]
uB   = 9.274009994e-24  # Bohr magneton [J/T]
c    = 299792458.0      # Speed of light [m/s]
kB   = 1.38064852e-23   # Boltzmann constant [J/K]

## TO DO: Get check boxes working!

def _to_pointer(numpy_array): 
    """
    Converts the supplied numpy_array (assumed to be the usual 64-bit float)
    to a pointer, allowing it ot be "connected" to the C-code engine. If None
    is supplied, this returns None.
    
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
    C-pointer to the first element, or None if numpy_array=None
    
    Examples
    --------
    my_solver.Bxs = _to_pointer(my_Bxs).
    
    """
    if numpy_array is None: return None
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
    
    def keys(self):
        """
        Returns a list of keys that can be used in set() or get().
        """
        return ['gamma', 'M', 'alpha', 'X', 'STT',
                'Tx', 'Ty', 'Tz', 'Bx', 'By', 'Bz',
                'Nxx', 'Nxy', 'Nxz', 'Nyx', 'Nyy', 'Nyz', 'Nzx', 'Nzy', 'Nzz']
    
    def set(self, key, value): 
        """
        Sets the specified parameter (key) to the specified value. Specifically,
        will be set by evaluating self.keyword = value. In the case of an array, 
        it will convert it to a pointer and (64-bit float) use an underscore, as needed, 
        saving a "local" copy of the supplied value, such that garbage 
        collection doesn't automatically delete it.

        Parameters
        ----------
        key:
            Parameter name to set, e.g. 'Bx'. 
        value:
            Value to set it to. Can be a number or ndarray.
        
        Example
        -------
        my_solver.a.set(By=numpy.linspace(0,5,my_solver.steps), gamma=27)

        Returns
        -------
        self

        """
    
        # Store the "local" copy, to make sure we keep it from getting
        # deleted by the garbage collection.
        exec('self.'+key+"s=v", dict(self=self,v=value))
        
        # If it's an array, convert it before setting and use _*s
        if type(value)==list: value = _n.ndarray(value)
        if type(value)==_n.ndarray:
            exec('self._'+key+"s=v", dict(self=self,v=_to_pointer(value)))
        
        # Otherwise it's a number, so we need to kill the old pointer and
        # array
        else: 
            exec('self.'+key+'=v', dict(self=self,v=value))
            exec('self.' +key+'s=None', dict(self=self)) # Local copy
            exec('self._'+key+"s=None", dict(self=self)) # Converted copy
            
        return self
    
    __setitem__ = set
    
    def set_multiple(self, **kwargs):
        """
        Sends all keyword arguments to self.set().
        """
        for k in kwargs: self[k] = kwargs[k]
    
    __call__ = set_multiple
    
    def get(self, key='Bx'):
        """
        Returns the specified parameter. Will return the array (e.g., Bxs) 
        if there is one, and the value if there is not.
        """
        # If it's an array, return that
        if hasattr(self, key+'s'):
            x = eval('self.'+key+'s', dict(self=self))
            if x is not None: return x
        
        # Otherwise just return the value.
        return eval('self.'+key, dict(self=self))
    
    __getitem__ = get
    
    def clear_arrays(self):
        """
        This sets all array pointers to NULL (None).
        """
        self['gamma'] = self.gamma
        self['M']     = self.M
        self['alpha'] = self.alpha
        self['X']     = self.X
        self['STT']   = self.STT
        
        self['Bx']    = self.Bx
        self['By']    = self.By
        self['Bz']    = self.Bz
        
        self['Nxx']   = self.Nxx
        self['Nxy']   = self.Nxy
        self['Nxz']   = self.Nxz
        self['Nyx']   = self.Nyx
        self['Nyy']   = self.Nyy
        self['Nyz']   = self.Nyz
        self['Nzx']   = self.Nzx
        self['Nzy']   = self.Nzy
        self['Nzz']   = self.Nzz
        
        self['Dxx']   = self.Dxx
        self['Dxy']   = self.Dxy
        self['Dxz']   = self.Dxz
        self['Dyx']   = self.Dyx
        self['Dyy']   = self.Dyy
        self['Dyz']   = self.Dyz
        self['Dzx']   = self.Dzx
        self['Dzy']   = self.Dzy
        self['Dzz']   = self.Dzz
    



class solver_api():
    """
    Scripted interface for the solver engine. 
    """
    
    _solver_keys  = ['dt', 'steps']
    
    def __init__(self):

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
        
    def set(self, key, value): 
        """
        Sets a parameter for the solver. Magnetic parameters, e.g., 'gamma' will
        be applied to both domains.

        Parameters
        ----------
        key:
            Parameter name to set, e.g. 'Bx'. 
        value:
            Value to set it to. Can be a number or ndarray.
        
        Example
        -------
        my_solver.a.set(By=numpy.linspace(0,5,my_solver.steps), gamma=27)

        Returns
        -------
        self

        """
        s = key.split('/')
        
        # If it's a property of the solver, store it in the solver
        if key in self._solver_keys:
            exec('self.'+key+"=v", dict(self=self,v=value))
    
        # If we've specified a domain
        elif s[0] == 'a': self.a.set(s[-1], value)
        elif s[0] == 'b': self.b.set(s[-1], value)
        
        # Otherwise, if it's something that we send to both a and b.
        elif key in self.a.keys():
            self.a.set(key, value)
            self.b.set(key, value)
        
        else: 
            print('OOPS api.set(): Cannot find key "'+key+'"')
        return self
    
    __setitem__ = set
    
    def set_multiple(self, **kwargs):
        """
        Sends all keyword arguments to self.set().
        """
        for k in kwargs: self[k] = kwargs[k]
    
    __call__ = set_multiple
    
    def get(self, key='a/Bx'):
        """
        Returns the specified parameter. Will return the array (e.g., Bxs) 
        if there is one, and the value if there is not.
        
        Parameters
        ----------
        key='a/Bx'
            Key of item to retrieve. Specify domain 'a' or 'b' with a '/' as in
            the example for domain-specific items, and just the parameter for
            the solver itself, e.g., 'steps'.
        
        Returns
        -------
        The value (or array if present).
        """
        
        # First split by '/'
        s = key.split('/')
        
        # If length is 1, this is a solver parameter
        if len(s) == 1: return eval('self.'+key)
        
        # Otherwise, it's a domain parameter
        domain = eval('self.'+s[0])
        key    = s[1]
        
        return domain.get(key)
    
    __getitem__ = get
    
    def run(self, update_initial_condition=False):
        """
        Creates the solution arrays and runs the solver to fill them up.
        Afterward, the initial conditions are (by default) set to the 
        last value of the solution arrays.
        
        NOTE: BE CAREFUL about update_initial_condition=True. This will ONLY
        transfer the 6 magnetization components, not the other arrays you have
        specified. For example, if you're using a Langevin (randomly fluctuating)
        field, you will need to transfer the last element of the field arrays to 
        the first element of the subsequent run's arrays manually, to maintain
        a consistent / continuous simulation. The last value of any user-
        specified array is used to calculate the last step.
        
        Parameters
        ----------
        update_initial_condition=True
            If True, the initial conditions (self.a.x0, self.b.x0, ...) will be
            set to the last value of the solution arrays. SEE ABOVE NOTE.
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
        
        # Previous values for the thermal field
        self._a_BTx = None
        self._a_BTy = None
        self._a_BTz = None
        
        self._b_BTx = None
        self._b_BTy = None
        self._b_BTz = None
        
        # Previous parameters relevant to thermal field
        self._dt       = None
        self._a_T      = None
        self._a_M      = None
        self._a_volume = None
        self._a_gamma  = None
        self._a_alpha  = None
        
        self._b_T      = None
        self._b_M      = None
        self._b_volume = None
        self._b_gamma  = None
        self._b_alpha  = None
        
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
        self.settings.add_parameter('solver/T',      0.0,  limits=(0,None),           siPrefix=True, suffix='K')
        self.settings.add_parameter('solver/reset', True)
        
        self.settings.add_parameter('a/mode', 1, limits=(0,1), tip='0=disabled, 1=LLG')
        
        self.settings.add_parameter('a/initial_condition/x0', 1.0, tip='Initial magnetization direction (will be normalized to unit length)')
        self.settings.add_parameter('a/initial_condition/y0', 1.0, tip='Initial magnetization direction (will be normalized to unit length)')
        self.settings.add_parameter('a/initial_condition/z0', 0.0, tip='Initial magnetization direction (will be normalized to unit length)')
        
        self.settings.add_parameter('a/material/gamma',      1.760859644e11, siPrefix=True, suffix='rad/(s*T)', tip='Magnitude of gyromagnetic ratio')
        self.settings.add_parameter('a/material/M',          1.0,  siPrefix=True, suffix='T', tip='Saturation magnetization (u0*Ms)')
        self.settings.add_parameter('a/material/volume', 100*50*3, siPrefix=False, suffix=' nm^3', tip='Volume of domain (nm^3). Relevant only for thermal and STT.')
        self.settings.add_parameter('a/material/alpha',      0.01, tip='Gilbert damping parameter')        
        
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
        self.settings.add_parameter('b/material/volume', 100*50*3, siPrefix=False, suffix=' nm^3', tip='Volume of domain (nm^3). Relevant only for thermal and STT.')
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
        self.initialize_plot_inspect()
        self.plot_inspect.after_clear = self.initialize_plot_inspect
        
        # 3D plot
        if _3d_enabled:
            self.tab_3d = self.tabs.add_tab('3D')
            self.button_3d_a   = self.tab_3d.place_object(_g.Button('a',   checkable=True, checked=True)) 
            self.button_3d_b   = self.tab_3d.place_object(_g.Button('b',   checkable=True, checked=False)) 
            self.button_3d_sum = self.tab_3d.place_object(_g.Button('Sum', checkable=True, checked=False)) 
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
            
            self.button_3d_a   .signal_clicked.connect(self._button_plot_3d_clicked)
            self.button_3d_b   .signal_clicked.connect(self._button_plot_3d_clicked)
            self.button_3d_sum .signal_clicked.connect(self._button_plot_3d_clicked)
            self.button_plot_3d.signal_clicked.connect(self._button_plot_3d_clicked)
            
        # Connect the other controls
        self.button_go     .signal_clicked.connect(self._button_go_clicked)
        
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

    def initialize_plot_inspect(self):
        """
        Clears the Inspect plot and populates it with empty arrays.
        """
        self.plot_inspect.clear()
        self.plot_inspect['t']  = []
        self.plot_inspect['ax'] = []
        self.plot_inspect['ay'] = []
        self.plot_inspect['az'] = []
        self.plot_inspect['bx'] = []
        self.plot_inspect['by'] = []
        self.plot_inspect['bz'] = []
        

    def _update_start_dots(self):
        """
        Gets the initial condition from the settings and updates the start
        dot positions.
        """
        if _3d_enabled:
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

    def _button_plot_3d_clicked(self, *a):
        """
        Plot 3d button pressed: Update the plot!
        """
        if _3d_enabled:
            d = self.plot_inspect
            self._trajectory_a_3d  .setData(pos=10*_n.vstack([d['ax'],d['ay'],d['az']]).transpose())
            self._trajectory_b_3d  .setData(pos=10*_n.vstack([d['bx'],d['by'],d['bz']]).transpose())
            self._trajectory_sum_3d.setData(pos=10*_n.vstack([d['ax']+d['bx'],d['ay']+d['by'],d['az']+d['bz']]).transpose())
            
            self._trajectory_a_3d  .setVisible(self.button_3d_a  .is_checked())
            self._trajectory_b_3d  .setVisible(self.button_3d_b  .is_checked())
            self._trajectory_sum_3d.setVisible(self.button_3d_sum.is_checked())
            
            self._start_dot_a_3d  .setVisible(self.button_3d_a  .is_checked())
            self._start_dot_b_3d  .setVisible(self.button_3d_b  .is_checked())
            self._start_dot_sum_3d.setVisible(self.button_3d_sum.is_checked())
            
            self.window.process_events()
    
    def _same_thermal_settings(self):
        """
        Checks the current thermal field settings against previous ones (if any)
        and returns true if they are the same.
        """
        return  self._dt       == self['dt']       and \
                self._T        == self['T']        and \
                self._a_M      == self['a/M']      and \
                self._a_volume == self['a/volume'] and \
                self._a_gamma  == self['a/gamma']  and \
                self._a_alpha  == self['a/alpha']  and \
                self._b_M      == self['b/M']      and \
                self._b_volume == self['b/volume'] and \
                self._b_gamma  == self['b/gamma']  and \
                self._b_alpha  == self['b/alpha']
       

    def _button_go_clicked(self, *a):
        """
        Go button pressed: Run the simulation!
        """

        for n in range(self.number_iterations.get_value()):
            self.label_iteration.set_text(str(n+1))
            self.run(self.settings['solver/reset'])
            
            
    def run(self, reset=False):
        """
        Run the specified simulation.
        
        IMPORTANT: See self.api.run()'s documentation for a concern about 
        setting reset=True.
        
        Parameters
        ----------
        reset=False
            After the run, update the initial conditions to match the 
            last point calculated. 
            
        Returns
        -------
        self
        """
        # Clear the api arrays, and transfer all the values from the 
        # TreeDictionary to the engine (API). This skips non-api entries 
        # like 'solver', 'a', 'b', and 'solver/T', and will transfer arrays
        # like 'a/Bx' if they exist in the plot_inspect.
        self.api.a.clear_arrays()
        self.api.b.clear_arrays()
        self._transfer_all_to_api()
        
        # If we have a Langevin field, add it in.
        if self['T'] > 0:
            
            # If domain a is enabled
            if self['a/mode']:
                
                # Generate new langevin field arrays
                BTx = self.thermal_field('a', self['T'], self['a/volume']*1e-27)
                BTy = self.thermal_field('a', self['T'], self['a/volume']*1e-27)
                BTz = self.thermal_field('a', self['T'], self['a/volume']*1e-27)
                
                # If thermal field-relevant parameters are the same
                # and we're doing a continuous simulation, and we have
                # previous thermal number, set it!
                if  self._same_thermal_settings() and \
                not self['reset'] and self._a_BTx is not None:
                    BTx[0] = self._a_BTx
                    BTy[0] = self._a_BTy
                    BTz[0] = self._a_BTz
                
                # Now add the thermal field to the existing field.
                self.api['a/Bx'] += BTx
                self.api['a/By'] += BTy
                self.api['a/Bz'] += BTz
                    
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
        if _3d_enabled: self._button_plot_3d_clicked()
            
        self.window.process_events()
        
        return self

    def _transfer_all_to_api(self):
        """
        Loops over the settings keys and transfers them to the api.
        """
        for key in self.settings.keys(): self._transfer(key)

    def _api_key_exists(self, key):
        """
        See if the key (from the TreeDictionary only!) exists in the API.
        """
        s = key.split('/')
        
        # Categories (solver, a, b)
        if len(s) == 1: return False
        
        # Solver parameters
        if s[-1] in self.api._solver_keys: return True
        
        # Domain parameters
        return s[-1] in self.api.a.keys()
        

    def _transfer(self, key):
        """
        Transfers the domain's parameter to the solver data.
        """
        # Ignore some keys that aren't in the api (defined above)
        if not self._api_key_exists(key): return        
        
        # s[0] is the domain, s[-1] is the parameter name
        s = key.split('/')

        # Array value
        short_key = s[0]+'/'+s[-1]
        if short_key in self.plot_inspect.ckeys: 
            
            # Use the array from the databox plotter
            value = self.plot_inspect[short_key]
            
            # Pluralize the variable name.
            s[-1] = s[-1]+'s'
                    
        # Normal value
        else: value = self.settings[key]
        
        # If it's under an unchecked category, zero it.
        if s[1] in ['applied_field', 'other_torques', 'anisotropy', 'dipole']:
            
            # If this is the category itself, break.
            if len(s) == 2: return
            
            # Otherwise it's a parameter. Zero the value if it's in an 
            # Unchecked category.
            elif not self.settings[s[0]+'/'+s[1]]: value = 0
                
        # Come up with the command for sending the parameter to the api
        if s[0]=='solver': command =            'self.api.set_multiple('+s[-1]+'= value)'
        else:              command = 'self.api.'+s[0]+  '.set_multiple('+s[-1]+'= value)'   
        
        # Try it!
        try:    exec(command, dict(self=self, value=value))
        except: print('FAIL: "'+command+'"')
        
        # Update the start dots in the 3D plot
        self._update_start_dots()
    
    def _elongate_domain_key(self, key):
        """
        Returns a key that is the long form (inserting any sub-headings) to 
        make it work on settings. Assumes it's of the form 'a/something' or
        'b/something'
        """
        split = key.split('/')
        
        for k in self.settings.keys():
            s = k.split('/')
            if s[0] == split[0] and s[-1] == split[-1]: return k
        
        print('UH OH, could not elongate "'+key+'"')
        return key
        
    def _get_root(self, key):
        """
        Returns the first root found for the specified short key, or None
        """
        for k in self.settings.keys():
            s = k.split('/')
            if s[-1] == key: return s[0]
            
        return None
        
    def set(self, key, value):
        """
        Sets self.settings[key] = value, unless value is an array. In that case
        it sets self.plot_inspect[key] = value. 
        
        You can also skip the sub-heading, so self.set('a/x0',0.5) will work 
        the same as self.set('a/initial_condition/x0', 0.5)

        Also, if you skip the root, it will assume either 'solver' or 'a' by
        default, so 'Tx' is the same as 'a/Tx'.

        Parameters
        ----------
        key : string
            Parameter key to set.
        value: string
            Parameter value to set.
            
        Returns
        -------
        self

        """
        s = key.split('/')
        
        # If we're using a shortcut key, like 'dt'
        if len(s) == 1 and s[0] not in ['solver', 'a', 'b']:
            
            # Insert the first root we find, if it's simple.
            s.insert(0,self._get_root(s[0]))
            key = '/'.join(s)
        
        # By this stage it better be length 2 or we're hosed.
        if len(s) < 2: 
            print('ERROR: Cannot set', key)
            return
        
        # Arrays go to the inspect plotter, values go to settings
        if type(value) == _n.ndarray: 
            s = key.split('/')
            self.plot_inspect[s[0]+'/'+s[-1]] = value
        
        # Otherwise it's just a value. Update the tree.
        else:                 
            self.settings[self._elongate_domain_key(key)] = value
        
        # Update plot
        self.plot_inspect.plot()
        self.window.process_events()
        
        return self
    
    __setitem__ = set    
    
    def set_multiple(self, **kwargs):
        """
        Sends all keyword arguments to self.set().
        """
        for k in kwargs: self[k] = kwargs[k]
    
    __call__ = set_multiple
    
    def get(self, key):
        """
        Returns the value (or array if available) for the specified key. Key
        can be "short", not including the sub-heading, e.g. 'a/x0' is the same
        as 'a/initial_condition/x0'.

        Parameters
        ----------
        key : string
            Key of parameter to retrieve.

        Returns
        -------
        value of parameter or array if available.

        """
        # For solver settings, return them with simple keys
        if key in ['dt', 'steps', 'reset', 'T']: return self.settings['solver/'+key]
        
        # s[0] is the domain, s[-1] is the parameter name
        s = key.split('/')

        # Don't do anything for the roots
        if len(s) < 2: 
            print('WHOOPS. WHOOPS. "'+key+'" is invalid or not specific enough (forgot the domain?)')
            return
        
        # Array value
        short_key = s[0]+'/'+s[-1]
        if short_key in self.plot_inspect.ckeys: return self.plot_inspect[short_key]
        else: return self.settings[self._elongate_domain_key(key)]
    
    __getitem__ = get
    
    def ns(self, n=0):
        """
        Returns an array of integer indices for the n'th simulation iteration.
        
        Specifically, it returns 
        
        _n.array(range(self['steps'])) + n*(self.settings['steps']-1)
          
        The '-1' in the above is because the last step of the previous iteration
        is used to initialize the first step of the following iteration.
          
        Parameters
        ----------
        n=0 [integer]
            Which simulation iteration this corresponds to.
        """
        return _n.array(range(self.settings['solver/steps'])) \
               + n*(self.settings['solver/steps']-1)

    def ts(self, n=0):
        """
        Returns a time array for the n'th iteration as per the solver settings.
        
        Specifically, it returns self.ns()*self.settings['solver/dt']
        
        Parameters
        ----------
        n=0 [integer]
            Which iteration of the simulation
        """
        return self.ns(n)*self.settings['solver/dt']

    def zeros(self):
        """
        Returns a self['steps']-length array of zeros.
        """
        return _n.zeros(self['steps'])

    def ones(self):
        """
        Returns a self['steps']-length array of ones.
        """
        return _n.ones(self['steps'])

    def pulse(self, n0, n1):
        """
        Returns an array of length self['steps'] that is zero everywhere except
        from n0 to n1-1.
        """
        z = self.zeros()
        z[n0:n1] = 1.0
        return z
    
    def sin(self, frequency_GHz=1.0, phase_rad=0.0, n=0):
        """
        Returns an appropriately sized array of sinusoidal oscillations at
        frequency f_GHz [GHz] with phase p_rad [radians]. The integer n adds
        the appropriate phase for the n'th iteration of the simulation, allowing
        you to generate a continuous sinusoid over many iterations.
        
        Specifically, the returned array is 
        sin(2*pi*frequency_GHz*self.ts(n) + phase_rad)
        
        Parameters
        ----------
        frequency_GHz [GHz]
            Frequency of oscillation in GHz.
        phase_rad [radians]
            Phase in radians.
        n=0 [integer]
            Which iteration this is. You can use this integer to make a 
            continuous sinusoid over multiple iterations of the solver. 
        """
        return _n.sin(2*_n.pi*frequency_GHz*1e9*self.ts(n)+phase_rad)
    
    def cos(self, frequency_GHz=1.0, phase_rad=0.0, n=0):
        """
        Returns an appropriately sized array of sinusoidal oscillations at
        frequency f_GHz [GHz] with phase p_rad [radians]. The integer n adds
        the appropriate phase for the n'th iteration of the simulation, allowing
        you to generate a continuous sinusoid over many iterations.
        
        Specifically, the returned array is 
        cos(2*pi*frequency_GHz*self.ts(n) + phase_rad)
        
        Parameters
        ----------
        f_GHz=1.0 [GHz]
            Frequency of oscillation in GHz.
        p_rad=0.0 [radians]
            Phase in radians.
        n=0 [integer]
            Which iteration this is. You can use this integer to make a 
            continuous sinusoid over multiple iterations of the solver. 
        """
        return _n.cos(2*_n.pi*frequency_GHz*1e9*self.ts(n)+phase_rad)
    
    def gaussian_noise(self, mean=0, standard_deviation=1.0):
        """
        Returns an appropriately sized (self['steps']) array of random values
        drawn from a Gaussian distribution of the specified mean and standard_deviation.
        
        Parameters
        ----------
        mean=0
            Center of the distribution.
        standard_deviation=1.0
            Standard deviation of the distribution.
        """
        return _n.random.normal(mean, standard_deviation, self.settings['solver/steps'])
        
    def thermal_field_rms(self, domain='a', T_K=295.0, volume_nm3=100*100*10):
        """
        Returns an array of length self['steps'] containing random, uncorrelated 
        fields [T], drawn from a Gaussian distribution appropriate for modeling
        thermal fluctuations at temperature T_K [K] for the specified domain,
        assuming a magnetic volume volume_nm3 [nm^3]. Note, that each component
        of the magnetic field should have an independent langevin array added
        to it to properly model thermal fluctuations (i.e., there should be
        three calls to this function per domain, per simulation).
        
        IMPORTANT: If you're conducting a "continuous" simulation that uses
        the last value of the previous run to initialize the first value of the 
        next run, you will need to also overwrite the first value of the 
        next langevin array with the last value of the previous run's array, 
        or else the simulation is not consistent (the last value of the field
        is used to calculate the last step).
        
        Parameters
        ----------
        domain='a'
            Which domain receives the Langevin field The only domain-
            specific parameters used are gamma and M.
            
        efficiency=1.0 [unitless]
            How electron spins are deposited on average per passing electron.
            
        volume_nm3=100*100*10 [nm^3]
            Total volume of the magnetic domain in cubic nanometers.
        
        Returns
        -------
        """
        return _n.sqrt(4*self.settings[domain+'/material/alpha']*kB*T_K     \
                     /  (self.settings[domain+'/material/gamma']*           \
                         self.settings[domain+'/material/M']/u0*volume_nm3* \
                         self.settings['solver/dt']))
    
    def thermal_field(self, domain='a', T_K=295.0, volume_nm3=100*100*10):
        """
        Returns an array of size self['steps'] containing random values drawn
        from a Gaussian distribution having standard deviation
        self.thermal_field_rms(domain, T_K, volume_nm3).
        
        IMPORTANT: To keep the solver consistent when transferring the previous run's 
        value to the initial condition of the current run, you must also transfer
        the previous run's last thermal field values to the first value of the
        current run.
        """
        return self.gaussian_noise()*self.thermal_field_rms(domain, T_K, volume_nm3)
    
    def spin_torque_per_mA(self, domain='a', efficiency=1.0, volume_nm3=100*100*10):
        """
        Returns the torque per milliamp [rad/(s*mA)] that would be applied when
        spin polarization is perpendicular to a magnetization 
        self[domain+'/M'] [T] filling the specified volume_nm3 [nm^3]. This can 
        basically be used to calculate the prefactor on the a x (a x s) term in
        the LLG equation, where a is the domain unit vector, and s is the 
        spin polarization unit vector, which, for this simulation has units of
        [rad/s].
        
        Parameters
        ----------
        domain='a'
            Which domain receives the spin transfer torque. The only domain-
            specific parameters used are gamma and M.
            
        efficiency=1.0 [unitless]
            How electron spins are deposited on average per passing electron.
            
        volume_nm3=100*100*10 [nm^3]
            Total volume of the magnetic domain in cubic nanometers.
        
        Returns
        -------
        The torque per mA [rad/(s*mA)] applied to a unit vector.
        """
        return efficiency*self[domain+'/gamma']*hbar*1e-3 / \
               (2*ec*(self[domain+'/M']/u0)*volume_nm3*1e-27)
    

if __name__ == '__main__':
    
    #######################
    # API Playground
    #######################
    
    # # Create a solver instance
    # m = solver_api()
    
    # # Set up the physical parameters
    # m.set_multiple(By=10.0, gamma=_n.pi*2, dt=0.005, zzz=300, steps=777, alpha=0.1/_n.pi/2)
    # m.a['By']  = _n.linspace(0,2,m.steps)
    # m.a['bz0'] = 2
    
    # # Run it & plot.
    # m.run(True); _s.plot.xy.data(None, [m.ax, m.by], clear=0, label=['ax','by'])
    # m.run(True); _s.plot.xy.data(None, [m.ax, m.by], clear=0, label=['ax','by'], xshift=m.steps-1, xshift_every=0)
    # m.run(True); _s.plot.xy.data(None, [m.ax, m.by], clear=0, label=['ax','by'], xshift=2*(m.steps-1), xshift_every=0)
    
    
    
    
    
    #############################
    # GUI Playground
    #############################

    self = solver()

    # # Pulse sequence
    # n1=1000; 
    # self['Tx'] = 4e9*self.pulse(0,n1); 
    # self['Tz'] = 1e9*self.pulse(n1,self['steps'])
    # self.run()

    # Thermal field
    self['T'] = 295.0
    self.run()
        
    
    
    
