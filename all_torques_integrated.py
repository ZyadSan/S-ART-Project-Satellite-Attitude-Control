import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class Constants:
    G = 6.67430e-11        # Gravitational const. (m^3 kg^-1 s^-2)
    M_EARTH = 5.972e24     # Earth mass (kg)
    R_EARTH = 6371e3       # Earth radius (m)
    MU = G * M_EARTH       # Gravitational parameter (m^3 s^-2)
    # Additional constants for new torques
    P_SUN = 4.56e-6        # Solar radiation pressure at Earth (N/m^2) = Φ/c
    C = 2.998e8            # Speed of light (m/s)
    K_B = 1.380649e-23     # Boltzmann constant (J/K)
    M_H2O = 2.99e-26       # Mass of H2O molecule (kg)
    T_DEFAULT = 273.15 + 20  # Default temperature (K) ~20°C
    SOLAR_FLUX = 1361.0    # Solar flux at Earth (W/m^2)


class Satellite:
    def __init__(self, mass=3.743, area=0.9336, inertia=None, 
                 # SRP parameters - asymmetric geometry creates torque
                 reflectance_diff=0.1, reflectance_spec=0.5, 
                 srp_area=0.01,  # Effective area for SRP (m^2)
                 # Lever arm: CP offset from CM (asymmetric)
                 cp_s=np.array([0.05, 0.02, 0.01]),  # Center of pressure SRP (m)
                 cm=np.array([0, 0, 0]),     # Center of mass (m)
                 # Aero parameters
                 cd=2.2,                       # Drag coefficient
                 aero_area=0.01,               # Effective area for aero (m^2)
                 cp_aero=np.array([0.05, 0.02, 0.01]),  # Center of pressure aero (m)
                 # Mass expulsion parameters
                 mdot=1e-6,                    # Mass flow rate (kg/s)
                 T_temp=None,                  # Temperature for mass expulsion
                 # CubeSat dimensions (1U = 10cm)
                 dimension=0.1):               # CubeSat side length (m)
        self.mass = mass       # kg
        self.area = area       # m^2
        self.inertia = inertia if inertia is not None else print("Inertia matrix must be provided")
        self.inertia_inv = np.linalg.inv(self.inertia)
        
        # SRP parameters
        self.reflectance_diff = reflectance_diff
        self.reflectance_spec = reflectance_spec
        self.srp_area = srp_area
        self.cp_s = cp_s
        self.cm = cm
        
        # Aero parameters
        self.cd = cd
        self.aero_area = aero_area
        self.cp_aero = cp_aero
        
        # Mass expulsion parameters
        self.mdot = mdot
        self.T_temp = T_temp if T_temp is not None else Constants.T_DEFAULT
        
        # CubeSat dimensions
        self.dimension = dimension
        
    def initial_state(self, altitude=550e3):
        r0 = np.array([Constants.R_EARTH + altitude, 0, 0])  # Initial position (m)
        v0 = np.array([0, np.sqrt(Constants.MU / np.linalg.norm(r0)), 0])  # Circular orbit velocity (m/s)
        return np.concatenate([r0, v0])


class Environment:
    def __init__(self):
        self.I_matrix = I_matrix

    def atmospheric_density(self, altitude):
        """
        Improved atmospheric density model using exponential approximation
        with different scale heights for different altitude regimes.
        Based on standard atmosphere models for LEO.
        """
        if altitude < 0:
            return 1.225  # Sea level
        
        # Standard exponential atmosphere model
        if altitude < 100e3:
            rho0 = 1.225
            H = 8500.0  # m
            rho = rho0 * np.exp(-altitude / H)
        elif altitude < 200e3:
            rho0 = 5.6e-7
            H = 28000.0  # m
            h_ref = 100e3
            rho = rho0 * np.exp(-(altitude - h_ref) / H)
        elif altitude < 400e3:
            rho0 = 2.5e-10
            H = 45000.0  # m
            h_ref = 200e3
            rho = rho0 * np.exp(-(altitude - h_ref) / H)
        else:
            rho0 = 3e-12
            H = 60000.0  # m
            h_ref = 400e3
            rho = rho0 * np.exp(-(altitude - h_ref) / H)
        
        return max(rho, 1e-20)  # Minimum density to avoid underflow

    def magnetic_torque(self, r):
        D = 0.2  # Spacecraft residual dipole moment (A·m²)
        M = 7.8e15  # Earth's magnetic moment × magnetic constant (T·m³)
        l = 2  # Worst-case: magnetic poles (strongest field)
        r_norm = np.linalg.norm(r)
        B = (M * l) / (r_norm)**3  # Magnetic field (Tesla)
        return D * B

    # ============================================================
    # 1. SOLAR RADIATION PRESSURE TORQUE
    # ============================================================
    def solar_radiation_pressure_torque(self, sat, r, sun_direction=None):
        """
        Solar Radiation Pressure Torque
        
        Using the simplified formula from the image:
        L_SRP = (Φ/c) * A_s * (1 + q) * (cp_s - cm) * cos(φ)
        
        Where:
        - Φ: solar flux (W/m^2) ~ 1361 W/m^2 at Earth
        - c: speed of light
        - A_s: exposed surface area
        - q: reflectance factor (q = specular reflectance)
        - cp_s: center of solar pressure
        - cm: center of mass
        - φ: angle between sun direction and surface normal
        
        Note: SRP torque is INDEPENDENT of altitude (at Earth orbit)
        """
        solar_flux = Constants.SOLAR_FLUX  # W/m^2
        
        # Default sun direction (from +X direction, hitting the face)
        if sun_direction is None:
            sun_direction = np.array([1.0, 0.0, 0.0])
        sun_direction = sun_direction / np.linalg.norm(sun_direction)
        
        # Surface normal (panel normal points outward in +X)
        surface_normal = np.array([1.0, 0.0, 0.0])
        
        # Angle between sun direction and surface normal
        cos_phi = np.dot(sun_direction, surface_normal)
        cos_phi = max(cos_phi, 0.0)
        
        # Reflectance factor q (using specular reflectance)
        q = sat.reflectance_spec
        
        # Lever arm: center of pressure - center of mass
        r_lever = sat.cp_s - sat.cm
        
        # Using the full formula from the image with reflectance terms:
        # F_SRP = -P_⊙ * S_i * [2*(R_diff/3 + R_spec*cos(θ))*n_B + (1-R_spec)*s] * max(cos(θ),0)
        
        P_sun = Constants.P_SUN  # N/m^2
        R_diff = sat.reflectance_diff
        R_spec = sat.reflectance_spec
        
        # Force vector components
        # Component along surface normal
        normal_component = 2 * (R_diff/3 + R_spec * cos_phi) * surface_normal
        # Component along sun direction  
        sun_component = (1 - R_spec) * sun_direction
        
        force_vector = -P_sun * sat.srp_area * (normal_component + sun_component) * max(cos_phi, 0)
        
        # Torque = r_lever × F
        torque = np.cross(r_lever, force_vector)
        
        return torque

    # ============================================================
    # 2. AERODYNAMIC TORQUE
    # ============================================================
    def aerodynamic_torque(self, sat, r, v=None):
        """
        Aerodynamic Torque
        
        From the image:
        L_aero = Σ(r^i × F_aero^i)
        
        F_aero^i = -1/2 * ρ * C_D * ||v_rel|| * v_rel * S_i * max(cos(θ_aero^i), 0)
        
        Where:
        - ρ: atmospheric density
        - C_D: drag coefficient
        - v_rel: relative velocity (spacecraft velocity - atmosphere velocity)
        - S_i: surface area of panel i
        - θ_aero^i: angle between velocity and surface normal
        """
        r_norm = np.linalg.norm(r)
        altitude = r_norm - Constants.R_EARTH
        
        # Atmospheric density
        rho = self.atmospheric_density(altitude)
        
        # Orbital velocity (circular orbit assumption)
        if v is None:
            v_orbit = np.sqrt(Constants.MU / r_norm)
            v = np.array([0.0, v_orbit, 0.0])
        
        # Relative velocity (atmosphere rotates with Earth)
        omega_earth = 7.2921159e-5  # rad/s
        v_atm = omega_earth * np.cross(np.array([0, 0, 1]), r)
        v_rel = v - v_atm
        v_rel_mag = np.linalg.norm(v_rel)
        
        if v_rel_mag < 1e-10:
            return np.array([0.0, 0.0, 0.0])
        
        v_rel_hat = v_rel / v_rel_mag
        
        # Surface normal: the face pointing INTO the velocity direction
        # For orbital motion in +Y direction, the face perpendicular to velocity
        surface_normal = -v_rel_hat  # Face pointing into the flow
        
        # For a flat plate perpendicular to flow, cos(θ) = 1
        cos_theta = 1.0
        
        # Drag force magnitude
        # F = -1/2 * ρ * C_D * ||v||^2 * A * max(cos(θ), 0) * direction
        force_magnitude = 0.5 * rho * sat.cd * v_rel_mag**2 * sat.aero_area * cos_theta
        
        # Force vector (opposes velocity)
        force_vector = -force_magnitude * v_rel_hat
        
        # Lever arm from center of mass to center of pressure
        r_lever = sat.cp_aero - sat.cm
        
        # Torque = r × F
        torque = np.cross(r_lever, force_vector)
        
        return torque

    # ============================================================
    # 3. MASS EXPULSION TORQUE
    # ============================================================
    def mass_expulsion_torque(self, sat, r, v=None):
        """
        Mass Expulsion Torque
        
        From the image:
        L_mexp = r × F_mexp = -ṁ * r × v_rel
        
        v_rel = sqrt(2*k_B*T / m_H2O) = 370 m/s (for typical outgassing)
        
        Where:
        - ṁ: mass flow rate (kg/s)
        - r: position vector from CM to thruster/nozzle (lever arm)
        - v_rel: exhaust velocity relative to spacecraft
        
        Note: Mass expulsion torque is INDEPENDENT of altitude
        """
        # Calculate exhaust velocity from thermal model
        # v_rel = sqrt(2 * k_B * T / m_H2O)
        v_rel_mag = np.sqrt(2 * Constants.K_B * sat.T_temp / Constants.M_H2O)
        
        # For outgassing, assume asymmetric venting from one corner/face
        expulsion_direction = np.array([1.0, 0.0, 0.0])
        
        # Relative velocity vector of expelled mass
        v_rel = v_rel_mag * expulsion_direction
        
        # Lever arm: distance from CM to outgassing point
        r_lever = np.array([0.05, 0.03, 0.02])  # Asymmetric offset
        
        # Torque: L = -ṁ * (r × v_rel)
        torque = -sat.mdot * np.cross(r_lever, v_rel)
        
        return torque


# ============================================================
# INERTIA MATRIX (from original code)
# ============================================================
I_matrix = np.array([
    [0.075, 0.0031, 0.0002],
    [0.0031, 0.1929, -3.2416e-5],
    [0.0002, -3.2416e-5, 0.1745],
])  # Inertia matrix calculated with m=3.743 kg


# ============================================================
# MAIN SIMULATION
# ============================================================
env = Environment()
sat = Satellite(inertia=I_matrix)

altitudes = []
torque_magnetic = []
torque_srp = []
torque_aero = []
torque_mexp = []

for h in np.linspace(0, 550e3, 550):
    r = np.array([Constants.R_EARTH + h, 0, 0])
    
    # Magnetic torque
    torque_mag = env.magnetic_torque(r)
    
    # Solar Radiation Pressure torque
    torque_s = env.solar_radiation_pressure_torque(sat, r)
    
    # Aerodynamic torque
    torque_a = env.aerodynamic_torque(sat, r)
    
    # Mass expulsion torque
    torque_m = env.mass_expulsion_torque(sat, r)
    
    altitudes.append(h / 1000)
    torque_magnetic.append(np.linalg.norm(torque_mag))
    torque_srp.append(np.linalg.norm(torque_s))
    torque_aero.append(np.linalg.norm(torque_a))
    torque_mexp.append(np.linalg.norm(torque_m))


# ============================================================
# PLOTTING ALL FOUR TORQUES - Individual plots
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('CubeSat Disturbance Torques vs Altitude (SART Project)', fontsize=14, fontweight='bold')

# 1. Magnetic Torque
axes[0, 0].plot(torque_magnetic, altitudes, 'b-', linewidth=1.5)
axes[0, 0].set_title('Magnetic Torque')
axes[0, 0].set_xlabel('Torque (Nm)')
axes[0, 0].set_ylabel('Altitude (km)')
axes[0, 0].grid(True)
axes[0, 0].ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

# 2. Solar Radiation Pressure Torque
axes[0, 1].plot(torque_srp, altitudes, 'orange', linewidth=1.5)
axes[0, 1].set_title('Solar Radiation Pressure Torque')
axes[0, 1].set_xlabel('Torque (Nm)')
axes[0, 1].set_ylabel('Altitude (km)')
axes[0, 1].grid(True)
axes[0, 1].ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

# 3. Aerodynamic Torque
axes[1, 0].plot(torque_aero, altitudes, 'g-', linewidth=1.5)
axes[1, 0].set_title('Aerodynamic Torque')
axes[1, 0].set_xlabel('Torque (Nm)')
axes[1, 0].set_ylabel('Altitude (km)')
axes[1, 0].grid(True)
axes[1, 0].ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

# 4. Mass Expulsion Torque
axes[1, 1].plot(torque_mexp, altitudes, 'r-', linewidth=1.5)
axes[1, 1].set_title('Mass Expulsion Torque')
axes[1, 1].set_xlabel('Torque (Nm)')
axes[1, 1].set_ylabel('Altitude (km)')
axes[1, 1].grid(True)
axes[1, 1].ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

plt.tight_layout()
plt.savefig('sart_torques_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

# ============================================================
# COMBINED COMPARISON PLOT (log scale for better comparison)
# ============================================================
plt.figure("Combined Torque Comparison", figsize=(10, 7))
plt.title('All Disturbance Torques vs Altitude (SART Project)', fontweight='bold')

torque_magnetic_np = np.array(torque_magnetic)
torque_srp_np = np.array(torque_srp)
torque_aero_np = np.array(torque_aero)
torque_mexp_np = np.array(torque_mexp)
altitudes_np = np.array(altitudes)

plt.semilogy(torque_magnetic_np, altitudes_np, 'b-', label='Magnetic', linewidth=2)
plt.semilogy(np.maximum(torque_srp_np, 1e-15), altitudes_np, 'orange', label='Solar Radiation Pressure', linewidth=2)
plt.semilogy(np.maximum(torque_aero_np, 1e-15), altitudes_np, 'g-', label='Aerodynamic', linewidth=2)
plt.semilogy(np.maximum(torque_mexp_np, 1e-15), altitudes_np, 'r-', label='Mass Expulsion', linewidth=2)
plt.xlabel('Torque (Nm) - Log Scale')
plt.ylabel('Altitude (km)')
plt.legend(loc='best')
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('sart_torques_combined_log.png', dpi=150, bbox_inches='tight')
plt.show()

# Print summary values at target altitude (550 km)
idx = -1
print("=" * 60)
print(f"TORQUE VALUES AT TARGET ALTITUDE (550 km)")
print("=" * 60)
print(f"Magnetic Torque:           {torque_magnetic[idx]:.4e} Nm")
print(f"Solar Radiation Pressure:  {torque_srp[idx]:.4e} Nm")
print(f"Aerodynamic Torque:        {torque_aero[idx]:.4e} Nm")
print(f"Mass Expulsion Torque:     {torque_mexp[idx]:.4e} Nm")
print("=" * 60)

# Also print at key altitudes
for h_target in [100, 200, 300, 400, 500, 550]:
    idx_h = int(h_target)
    if idx_h < len(altitudes):
        print(f"\nAt {h_target} km:")
        print(f"  Magnetic: {torque_magnetic[idx_h]:.4e} Nm")
        print(f"  SRP:      {torque_srp[idx_h]:.4e} Nm")
        print(f"  Aero:     {torque_aero[idx_h]:.4e} Nm")
        print(f"  MExp:     {torque_mexp[idx_h]:.4e} Nm")
        print(f"  Density:  {env.atmospheric_density(h_target*1000):.4e} kg/m^3")