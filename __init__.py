import ctypes as _c
import numpy  as _n
import os     as _os
from sys import platform as _platform


# Find the path to the compiled c-code (only Windows and Linux supported so far.)
if _platform in ['win32']: _path_dll = _os.path.join(_os.path.split(__file__)[0],'macrospin.dll')
else:                      _path_dll = _os.path.join(_os.path.split(__file__)[0],'macrospin.so')

class _data(_c.Structure):
    """
    Structure for sending all the simulation parameters to the c library.
    """
    _fields_ = [
        # Physical system
        ('time_scale',     _c.c_double), # 1e-9 for nanoseconds
        ('T',              _c.c_double), # Temperature (K)
        ('thickness_fm',   _c.c_double), # Thickness of ferromagnetic layer (nm)
        ('thickness_nm',   _c.c_double), # Thickness of normal metal (nm)
        ('resistivity_fm', _c.c_double), # Resistivity of ferromagnet (Ohm/m^2)
        ('resistivity_nm', _c.c_double), # Resistivity of normal metal (Ohm/m^2)
        ('length',         _c.c_double), # Length along x (nm)
        ('width',          _c.c_double), # Width along y (nm)
        ('ms',             _c.c_double), # Saturation magnetization (T) used only for spin transfer efficiency and Langevin
        ('Byx',            _c.c_double), # In-plane anisotropy (T)
        ('Bzx',            _c.c_double), # Out-of-plane anisotropy (T)
        ('damping',        _c.c_double), # Gilber damping (unitless)
    
        # Torque stuff
        ('torque_type',    _c.c_int),    # 1 for spin Hall, 2 for sinusoid
        ('hall_angle',     _c.c_double), # Conversion from charge to spin current density (which also has units of "Amps/m^2")
        ('g_factor',       _c.c_double), # 1e-9 for nanoseconds
        ('efficiency',     _c.c_double), # 1e-9 for nanoseconds
        ('gyromagnetic_magnitude', _c.c_double), # Magnitude of the gyromagnetic ratio, nominally 8.794033700300592e10*d->g_factor (rad / s T)
        
        # Current drive
        ('I0',             _c.c_double), # DC current (mA)
        ('I1',             _c.c_double), # AC amplitude (mA)
        ('f1',             _c.c_double), # AC frequency (GHz)
        ('Bx_per_mA',      _c.c_double), # Conversion from current to field (T/mA)
        ('By_per_mA',      _c.c_double),
        ('Bz_per_mA',      _c.c_double),
    
        # Applied field
        ('B0',             _c.c_double), # DC field (T)
        ('Bx_hat',         _c.c_double), # Unit vector components
        ('By_hat',         _c.c_double),
        ('Bz_hat',         _c.c_double),
        
        # Solver stuff
        ('steps',  _c.c_int),    # Number of steps
        ('t0',     _c.c_double), # Solver starting time (nanoseconds if time_scale=1e-9)
        ('dt',     _c.c_double), # Solver time step
        
        # Instantaneous values (don't mess with)
        ('n',   _c.c_int),
        ('mx',  _c.c_double),
        ('my',  _c.c_double),
        ('mz',  _c.c_double),
        ('Mx',  _c.c_double),
        ('My',  _c.c_double),
        ('Mz',  _c.c_double),
        ('Bx',  _c.c_double),
        ('By',  _c.c_double),
        ('Bz',  _c.c_double),
        ('I',   _c.c_double),
        ('Js',  _c.c_double),
        
        # Solution arrays
        ("solution_mx", _c.POINTER(_c.c_double)),
        ("solution_my", _c.POINTER(_c.c_double)),
        ("solution_mz", _c.POINTER(_c.c_double)),
        
        # Bureaucracy
        ("log_level",   _c.c_int)
        
    ] # End of _data structure


class solver():
    """
    Class for interfacing with the macrospin simulation c-library.
    
    **kwargs are sent to self.set().
    
    Parameters
    ----------
    
    T = 0 (Kelvin)
        Temperature used for Langevin field.

    length = 6000 (nanometers)
        Length (along x) of all device layers. 
        
    width = 500 (nanometers)
        Width (along y) of all device layers. 

    thickness_nm = 10 (nanometers)
        Thickness (along z) of a normal metal film on top of the 
        ferromagnetic film. This, combined with the layer resistivities, 
        is used only to convert mA of current into transverse spin current 
        in the spin Hall torque mode (torque_type=1).    
    
    thickness_fm = 10 (nanometers)
        Thickness (along z) of the "free layer" ferromagnet.
            
    ms = 0.8 (Tesla)
        Saturation magnetization of the "free layer" ferromagnet. 
    
    Byx = 0.005 (Tesla)
        Coercive field along the y direction. Nominally equal to
        ms*(Nyy-Nxx), where Nxx+Nyy+Nzz = 1 are the elements of the (assumed
        diagonal) demagnetization tensor.
        
    Bzx = 0.8 (Tesla)
        Coercive field along the z-direction. Nominally equal to
        ms*(Nzz-Nxx).
        
    damping = 0.02 (unitless)
        Gradient damping parameter. The Gilbert damping parameter
        works here, but the exact direction of the torque may be slightly 
        different. This difference is small if damping << 1.

    resistivity_fm = 65.2 (units not needed)
        Resistivity of the ferromagnetic film. Only the 
        ratio of resistivities matters for determining the spin Hall torque.
        
    resistivity_nm = 21.9 (units not needed)
        Resistivity of the normal metal film. Only the 
        ratio of resistivities matters for determining the spin Hall torque.
            
    torque_type = 1 (flag)
        How to convert current (mA) into torque. 1=spin Hall (current along x), 
        2=sinusoidal (current along z).
        
    hall_angle = 0.05 (unitless)
        Spin Hall angle for torque_type=1.
        
    g_factor = 2.002319 (unitless)
        Electron g-factor, used only to calculate spin torque efficiency. 
    
    efficiency = 1.0 (unitless)
        Fraction of maximum per-electron spin transfer torque that is actually
        applied to the free layer.
    
    gyromagnetic_magnitude = 176084607648.0 (radians / second Tesla)
        Magnitude of the gyromagnetic ratio
    
    I0 = 0.0 (milliamps)
        DC bias.
        
    I1 = 0.0 (milliamps)
        RF current amplitude.
        
    f1 = 2.0 (Hertz scaled by time_scale)
        RF current frequency.
        
    Bx_per_mA = 0.0, By_per_mA = 0.0, Bz_per_mA = 0.0  (Tesla/milliamp)
        How much field in the x, y, and z directions are applied per milliamp
        of current.
    
    B0 = 0.02 (Tesla)
        Static applied field magnitude
        
    Bx_hat = 0, By_hat = 1, Bz_hat = 0 (unitless)
        Unit vector direction of applied field.
        
    mx = 1, my = 0, mz = 0 (unitless)
        Unit vector direction of the free layer magnetization.
    
    Mx = 0, My = 1, Mz = 0 (unitless)
        Unit vector direction of the fixed layer magnetization. Not used for
        spin Hall torque (torque_type=1).
    
    log_level = 0 (flag)
        How much information to include in the log output. This is mostly for
        Jack to mess with, and won't produce reliable results.
        
    steps = 10000
        How many time steps per simulation.
    
    t0 = 0 (seconds scaled by time_scale)
        Starting time of simulation. Used, e.g., to calculate the RF current.
        
    dt = 0.0005 ( scaled by time_scale)
        Time step.
        
    time_scale = 1e-9 (seconds per time value)
        Natural units for the time-domain. 
        1e-9 means the system works with nanoseconds and gigahertz.
        
    Notes
    -----
    Device dimensions and the saturation magnetization ms are used to estimate 
    the total angular momentum of the free layer, which affects the per-spin torque efficiency
    and langevin field. Shape anisotropy is specified by the coercive fields
    Byx and Bzx. The cross section thickness_fm*width is also used
    to calculate the current density for the spin Hall torque (torque_type=1).
    
    Applied current is I(t) = I0 + I1*cos(2*pi*f1*t).

    """
    
    def __init__(self, **kwargs):

        # Load the shared library / dll
        self._macrospin = _c.CDLL(_path_dll)

        # Create the settings structure, and set the default values.
        self.data = _data(
        
            T               = 0,
            thickness_fm    = 10,
            thickness_nm    = 10,
            resistivity_fm  = 65.2,
            resistivity_nm  = 21.9,
            length          = 6000,
            width           = 500,
            ms              = 0.8,
            Byx             = 0.005,
            Bzx             = 0.8,
            damping         = 0.02,
            
            torque_type     = 1,
            hall_angle      = 0.05,
            g_factor        = 2.002319,
            efficiency      = 1.0,
            gyromagnetic_magnitude = 176084607648.0,
            
            I0              = 0.0,
            I1              = 0.0,
            f1              = 2.0,
            Bx_per_mA       = 0.0,
            By_per_mA       = 0.0,
            Bz_per_mA       = 0.0,
            
            B0              = 0.02,
            Bx_hat          = 0,
            By_hat          = 1,
            Bz_hat          = 0,
            
            mx = 1,
            my = 0,
            mz = 0,
            Mx = 0,
            My = 1,
            Mz = 0,
            
            log_level       = 0,
            steps           = 10000,
            t0              = 0,
            dt              = 0.0005,
            time_scale      = 1e-9)
        
        # Solution arrays
        self.solution_mx = None
        self.solution_my = None
        self.solution_mz = None
        
        # Update settings
        self.set(**kwargs)
        
    def __repr__(self):
        """
        Returns the string that appears when you inspect the object.
        """
        s = "Macrospin settings:\n\n"
        
        s = s+"  T              = %0.5G K\n"  % self['T']
        s = s+"  thickness_fm   = %0.5G nm\n" % self['thickness_fm']
        s = s+"  thickness_nm   = %0.5G nm\n" % self['thickness_nm']
        s = s+"  resistivity_fm = %0.5G (no units needed)\n" % self['resistivity_fm']
        s = s+"  resistivity_nm = %0.5G (no units needed)\n" % self['resistivity_nm']
        s = s+"  length         = %0.5G nm\n" % self['length']
        s = s+"  width          = %0.5G nm\n" % self['width']
        s = s+"  ms             = %0.5G T\n" % self['ms']
        s = s+"  Byx            = %0.5G T\n" % self['Byx']
        s = s+"  Bzx            = %0.5G T\n" % self['Bzx']
        s = s+"  damping        = %0.5G (unitless)\n" % self['damping']
        s = s+"\n"
        
        s = s+"  torque_type            = %i\n"    % self['torque_type']
        s = s+"  hall_angle             = %0.5G\n" % self['hall_angle']
        s = s+"  g_factor               = %0.5G\n" % self['g_factor']
        s = s+"  efficiency             = %0.5G\n" % self['efficiency']
        s = s+"  gyromagnetic_magnitude = %0.5G\n" % self['gyromagnetic_magnitude']
        s = s+"\n"
        
        s = s+"  I0        = %0.5G mA\n" % self['I0']
        s = s+"  I1        = %0.5G mA\n" % self['I1']
        s = s+"  f1        = %0.5G / time_scale Hz\n" % self['f1']
        s = s+"  Bx_per_mA = %0.5G T/mA\n" % self['Bx_per_mA']
        s = s+"  By_per_mA = %0.5G T/mA\n" % self['By_per_mA']
        s = s+"  Bz_per_mA = %0.5G T/mA\n" % self['Bz_per_mA']
        s = s+"\n"
        
        s = s+"  B0        = %0.5G T\n" % self['B0']
        s = s+"  Bx_hat    = %0.5G\n" % self['Bx_hat']
        s = s+"  By_hat    = %0.5G\n" % self['By_hat']
        s = s+"  Bz_hat    = %0.5G\n" % self['Bz_hat']
        s = s+"\n"
        
        s = s+"  mx = %0.5G\n" % self['mx']
        s = s+"  my = %0.5G\n" % self['my']
        s = s+"  mz = %0.5G\n" % self['mz']
        s = s+"\n"
        
        s = s+"  Mx = %0.5G\n" % self['Mx']
        s = s+"  My = %0.5G\n" % self['My']
        s = s+"  Mz = %0.5G\n" % self['Mz']
        s = s+"\n"
        
        s = s+"  log_level  = %0.5G\n" % self['log_level']
        s = s+"  steps      = %0.5G\n" % self['steps']
        s = s+"  t0         = %0.5G * time_scale s\n" % self['t0']
        s = s+"  dt         = %0.5G * time_scale s\n" % self['dt']
        s = s+"  time_scale = %0.5G\n" % self['time_scale']
        
        return s
    
    
    def get(self, k='T'):
        """
        Returns the setting for the supplied keyword (string) or list of keywords.
        """
        # Simple call
        if type(k) == str: return eval("self.data."+k, dict(self=self))
        
        # List call
        out = []
        for x in k: out.append(eval("self.data."+x, dict(self=self)))
        return out
        
    def set(self, **kwargs):
        """
        Updates the simulation settings based on the supplied keywords.
        """
        for k in kwargs.keys(): exec("self.data."+k+"=kwargs[k]", dict(self=self, kwargs=kwargs, k=k))
        return self

    def __call__(self, *args, **kwargs):
        """
        If a string argument is set, returns returns the value of that parameter.
        If keyword arguments are provided, calls set() and returns self.
        """
        if len(args): return self.get(args[0])
        else:         return self.set(**kwargs)

    def __getitem__(self, k)   : return self.get(k)
    def __setitem__(self, k, v): return self.set(**{k:v})
    
    def set_angles(self, m_theta=None, m_phi=None, M_theta=None, M_phi=None, B_theta=None, B_phi=None, **kwargs):
        """
        Set's the free (m), fixed (M), and field (B) angles (degrees). Items 
        set to None are ignored or set to zero (if only one is specified).
        
        For example, set_angles(B_theta=90, t0=0) will set B_theta=90, B_phi=0 
        and reset the clock.
        
        **kwargs are sent to self.set() at the end.
        
        Parameters
        ----------
        m_theta, M_theta, B_theta:
            Angle away from the x-axis (degrees).
        m_phi, M_phi, B_phi:
            Angle around the x-axis (degrees), with 0 in the +y direction.
        """
        to_radians = _n.pi / 180.0

        # Free layer        
        if not m_theta == None or not m_phi == None:
            if m_theta == None: m_theta = 0
            if m_phi   == None: m_phi   = 0
            self['mx'] = _n.cos(m_theta*to_radians)
            self['my'] = _n.sin(m_theta*to_radians)*_n.cos(m_phi*to_radians)
            self['mz'] = _n.sin(m_theta*to_radians)*_n.sin(m_phi*to_radians)

        # Fixed layer        
        if not M_theta == None or not M_phi == None:
            if M_theta == None: M_theta = 0
            if M_phi   == None: M_phi   = 0
            self['Mx'] = _n.cos(M_theta*to_radians)
            self['My'] = _n.sin(M_theta*to_radians)*_n.cos(M_phi*to_radians)
            self['Mz'] = _n.sin(M_theta*to_radians)*_n.sin(M_phi*to_radians)

        # Field
        if not B_theta == None or not B_phi == None:
            if B_theta == None: B_theta = 0
            if B_phi   == None: B_phi   = 0
            self['Bx_hat'] = _n.cos(B_theta*to_radians)
            self['By_hat'] = _n.sin(B_theta*to_radians)*_n.cos(B_phi*to_radians)
            self['Bz_hat'] = _n.sin(B_theta*to_radians)*_n.sin(B_phi*to_radians)
        
        # Update other parameters
        self.set(**kwargs)
        
        return self
    
    def run(self, **kwargs):
        """
        Creates the solution arrays and runs the solver, which updates the 
        dynamic quantities, and stores the free layer trajectory in 
        self.solution_mx, self.solution_my, and self.solution_mz.
        
        **kwargs are sent to self.set() prior to execution.
        """
        
        self.set(**kwargs)
        
        # Create the solution arrays
        self.solution_mx = _n.zeros(self['steps'])
        self.solution_my = _n.zeros(self['steps'])
        self.solution_mz = _n.zeros(self['steps'])
        
        # Store them in c-ready form
        self.data.solution_mx = _n.ctypeslib.as_ctypes(self.solution_mx)
        self.data.solution_my = _n.ctypeslib.as_ctypes(self.solution_my)
        self.data.solution_mz = _n.ctypeslib.as_ctypes(self.solution_mz)
        
        # Solve it.
        self._macrospin.Solve_Huen.restype = None
        self._macrospin.Solve_Huen(_c.byref(self.data))
        
        return self



if __name__ == '__main__':

    import pylab
    
    # Create the solver object, and initialize the magnetization
    m = solver(steps=2000).set_angles(m_theta=45)
    
    # Get a time-domain array
    t = _n.linspace(m['t0'],m['t0']+(m['steps']-1)*m['dt'], m['steps'])
    
    # Run the solver.
    m.run()
    pylab.plot(t, m.solution_mx, label='mx')
    pylab.plot(t, m.solution_my, label='my')
    pylab.plot(t, m.solution_mz, label='mz')
    
    # Get a time-domain array
    t = _n.linspace(m['t0'],m['t0']+(m['steps']-1)*m['dt'], m['steps'])
    
    m.run()
    pylab.plot(t, m.solution_mx, label='mx')
    pylab.plot(t, m.solution_my, label='my')
    pylab.plot(t, m.solution_mz, label='mz')
    
    pylab.legend()
    pylab.xlabel('Time (ns)')
    pylab.ylabel('Magnetization')
    
