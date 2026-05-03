import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class Constants:
    G = 6.67430e-11        # Gravitational const. (m^3 kg^-1 s^-2)
    M_EARTH = 5.972e24     # Earth mass (kg)
    R_EARTH = 6371e3       # Earth radius (m)
    MU = G * M_EARTH       # Gravitational parameter (m^3 s^-2)


class Satellite:
    def __init__(self, mass=3.743, area=0.9336, inertia=None):
        self.mass = mass       # kg
        self.area = area       # m^2
        self.inertia = inertia if inertia is None else print("Inertia matrix must be provided") # kg m^2
        self.inertia_inv = np.linalg.inv(self.inertia)
       
    def initial_state(self, altitude=550e3):
        r0 = np.array([Constants.R_EARTH + altitude, 0, 0])  # Initial position (m)
        v0 = np.array([0, np.sqrt(Constants.MU / np.linalg.norm(r0)), 0])  # Circular orbit velocity (m/s)
        return np.concatenate([r0, v0])
class Environment:
    def __init__(self):
        self.I_matrix = I_matrix

    def magnetic_torque(self, r):
        D=0.2 # D is the spacecraft’s residual dipole moment in A.m²
        M=7.8e15 # M is the magnetic moment of the Earth multiplied by the magnetic constant in tesla• m3
        l=2 # For worst-case disturbance analysis: Take lambda = 2. This corresponds to 
# the magnetic poles where the Earth's field is strongest,
# resulting in the highest possible residual torque on your CubeSat
        r_norm = np.linalg.norm(r)
        r_hat = r / r_norm
        B=(M*l)/(r_norm)**3 #Magnetic field in Tesla
        return D*B
I_matrix = np.array([
    [0.075, 0.0031, 0.0002],
    [0.0031, 0.1929, -3.2416e-5],
    [0.0002, -3.2416e-5, 0.1745],]
    ) # Inertia matrix calculated with m=3.743 kg
env = Environment()
altitudes = []
torque_values = []
for h in np.linspace(0, 550e3, 550):
    r = np.array([Constants.R_EARTH + h, 0, 0])
    torque_magnetism = env.magnetic_torque(r)
    altitudes.append(h / 1000)
    torque_values.append(np.linalg.norm(torque_magnetism))
plt.figure("Magnetic Torque")
plt.title("Magnetic Torque vs Altitude")
plt.plot(torque_values, altitudes)
plt.xlabel("Magnetic Torque (Nm)")
plt.ylabel("Altitud (km)")
plt.grid(True)
plt.show()