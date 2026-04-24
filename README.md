# S-ART-Project-Satellite-Attitude-Control
The following are all the files that are related to the attitude determination and control system of the S-ART (Student Astronomical Radio Telescope) project's cubesat. This repository contains a the file required to run a dynamic simulator between the following subsystems:

• Orbital Propagation: Handles eccentric orbits considering the Earth's J2 flattening (precession) and atmospheric
braking using the NRLMSISE-00 model (a standard orbital decay model).
• Draft of Operating Modes: A document with all the operating modes the satellite should have.
• Rigid Body Dynamics (Real Effects): Euler's equations solver. It is used to see how the satellite rotates on its own. It propagates the following external perturbation torques:
  • Gravity Gradient Torque
  • Magnetic Torque
  • Solar Radiation Pressure Torque
  • Aerodynamic Torque
  • Mass Expulsion Torque 
• Reference Kinematics: Calculating, using pure geometry, the angular velocity w(ref) that the satellite would need to remain fixed, pointing towards the Nadir, a ground station, or the Earth's limb. 
