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
        ('Byx1',           _c.c_double), # In-plane anisotropy, element 1 (T)
        ('Bzx1',           _c.c_double), # Out-of-plane anisotropy, element 1 (T)
        ('Byx2',           _c.c_double), # In-plane anisotropy, element 2 (T)
        ('Bzx2',           _c.c_double), # Out-of-plane anisotropy, element 2 (T)
        ('Bex1',           _c.c_double), # Exchange-like field felt by element 1 from element 2 (T)
        ('Bex2',           _c.c_double), # Exchange-like field felt by elemnet 2 from element 1 (T)
        ('damping1',       _c.c_double), # Gradient damping, element 1 (unitless)
        ('damping2',       _c.c_double), # Gradient damping, element 2 (unitless)    
        ('g0',             _c.c_double), # Magnitude of the gyromagnetic ratio, nominally 8.794033700300592e10*g_factor (rad / s T)
            
        # Solver stuff
        ('steps',          _c.c_int),    # Number of steps
        ('dt',             _c.c_double), # Time step (seconds)
        
        # User-supplied arrays
        ("applied_Bx1", _c.POINTER(_c.c_double)),
        ("applied_By1", _c.POINTER(_c.c_double)),
        ("applied_Bz1", _c.POINTER(_c.c_double)),
        ("applied_Bx2", _c.POINTER(_c.c_double)),
        ("applied_By2", _c.POINTER(_c.c_double)),
        ("applied_Bz2", _c.POINTER(_c.c_double)),
        
        ("applied_tx1", _c.POINTER(_c.c_double)),
        ("applied_ty1", _c.POINTER(_c.c_double)),
        ("applied_tz1", _c.POINTER(_c.c_double)),
        ("applied_tx2", _c.POINTER(_c.c_double)),
        ("applied_ty2", _c.POINTER(_c.c_double)),
        ("applied_tz2", _c.POINTER(_c.c_double)),
        
        # Solution arrays
        ("Bx1", _c.POINTER(_c.c_double)),
        ("By1", _c.POINTER(_c.c_double)),
        ("Bz1", _c.POINTER(_c.c_double)),
        ("Bx2", _c.POINTER(_c.c_double)),
        ("By2", _c.POINTER(_c.c_double)),
        ("Bz2", _c.POINTER(_c.c_double)),
        
        ("tx1", _c.POINTER(_c.c_double)),
        ("ty1", _c.POINTER(_c.c_double)),
        ("tz1", _c.POINTER(_c.c_double)),
        ("tx2", _c.POINTER(_c.c_double)),
        ("ty2", _c.POINTER(_c.c_double)),
        ("tz2", _c.POINTER(_c.c_double)),
        
        ("mx1", _c.POINTER(_c.c_double)),
        ("my1", _c.POINTER(_c.c_double)),
        ("mz1", _c.POINTER(_c.c_double)),
        ("mx2", _c.POINTER(_c.c_double)),
        ("my2", _c.POINTER(_c.c_double)),
        ("mz2", _c.POINTER(_c.c_double)),
        
        # Internal values used by solver
        ('dmx1', _c.c_double), # Time step (seconds)
        ('dmy1', _c.c_double), # Time step (seconds)
        ('dmz1', _c.c_double), # Time step (seconds)
        ('dmx2', _c.c_double), # Time step (seconds)
        ('dmy2', _c.c_double), # Time step (seconds)
        ('dmz2', _c.c_double), # Time step (seconds)
        
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
    
    Byx1, Byx2 = 0.005 (Tesla)
        Coercive field along the y direction for each magnetic element. Nominally equal to
        ms*(Nyy-Nxx), where Nxx+Nyy+Nzz = 1 are the elements of the (assumed
        diagonal) demagnetization tensor.
        
    Bzx1, Bzx2 = 0.8 (Tesla)
        Coercive field along the z-direction for each magnetic volume. Nominally equal to
        ms*(Nzz-Nxx).
        
    Bex1, Bex2 = 0 (Tesla)
        Exchange-like field experienced by elements 1 and 2, respectively, 
        parallel to magnetizations 2 and 1, respectively.
        
    damping1, damping2 = 0.02 (unitless)
        Gradient damping parameter. The Gilbert damping parameter
        works here, but the exact direction of the torque may be slightly 
        different. This difference is small if damping << 1.

    g1, g2 = 176084607648.0 (radians / second Tesla)
        Magnitude of the gyromagnetic ratio for each magnetic element.
    
    steps = 10000
        How many time steps per simulation.
    
    dt = 0.0005e-9 (seconds)
        Time step.
        
    applied_Bx1, applied_By1, applied_Bz1 = 0 (Tesla)
    applied_Bx2, applied_By2, applied_Bz2 = 0 (Tesla)
        Externally applied field for each layer. Can specify a value or 
        an array of length "steps" defined above to create a time-dependent
        field.
    
    applied_tx1, applied_ty1, applied_tz1 = 0 (inverse seconds)
    applied_tx2, applied_ty2, applied_tz2 = 0 (inverse seconds)
        Externally applied torque for each layer in units of how much the 
        magnetization unit vectors change per unit time. Can specify a value or
        an array of length "steps" defined above to create a time-dependent
        torque.
        
    log_level = 0 (flag)
        How much information to include in the log output. This is mostly for
        Jack to mess with, and won't produce reliable results.
    """
    
    def __init__(self, **kwargs):

        # Load the shared library / dll
        self._macrospin = _c.CDLL(_path_dll)


#############################################

        # Create the settings structure, and set the default values.
        self.data = _data(
        
            Byx1            = 0.005,
            Bzx1            = 0.8,
            Byx2            = 0.005,
            Bzx2            = 0.8,
            
            Bex1            = 0,
            Bex2            = 0,
            
            damping1        = 0.01,
            damping2        = 0.01,
            
            g1              = 176084607648.0,
            g2              = 176084607648.0,
            
            log_level       = 0,
            steps           = 10000,
            dt              = 0.0005e-9)
        
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
        s = s+"  damping        = %0.5G K\n" % self['damping']
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
    
