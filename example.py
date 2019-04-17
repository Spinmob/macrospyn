import pylab
import numpy as np
import macrospyn
    
# Create the solver object, and initialize the magnetization
m = macrospyn.solver(steps=2000).set_angles(m_theta=45)

# Get a time-domain array
t = np.linspace(m['t0'],m['t0']+(m['steps']-1)*m['dt'], m['steps'])

# Run the solver.
m.run()
pylab.plot(t, m.solution_mx, label='mx')
pylab.plot(t, m.solution_my, label='my')
pylab.plot(t, m.solution_mz, label='mz')

# Get a time-domain array
t = np.linspace(m['t0'],m['t0']+(m['steps']-1)*m['dt'], m['steps'])

m.run()
pylab.plot(t, m.solution_mx, label='mx')
pylab.plot(t, m.solution_my, label='my')
pylab.plot(t, m.solution_mz, label='mz')

pylab.legend()
pylab.xlabel('Time (ns)')
pylab.ylabel('Magnetization')