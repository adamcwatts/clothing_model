# currently modeling wool with parameters from Neves et al. 2015 paper
# Default values for ODE.IVP are 1e-3 for rtol and 1e-6 for atol which are too small (sinusoidal behavior).
# Proposed tolerances are rtol=1e-5, atol=1e-7
TIME_INPUT_PARAMETERS = \
    {
        'start': 0,  # [time in seconds]
        'finish min': 30.0,  # [minutes]
    }

FABRIC_INPUT_PARAMETERS = \
    {'regain': 0.15,  # [fraction from 0 to 1]
     'porosity of air in fabric': 0.931,  # [fraction from 0 to 1]
     'dry fiber density': 1300,  # [kg/m^3]
     'fiber specific heat': 1360,  # [J/ Kg K ]
     # 'fabric thickness': 0.00857,  # [m]
     'fabric thickness': 0.01,  # [m]
     'R_ef': 23.4,  # [m^2 Pa / W]
     'thermal conductivity of fiber': 0.20  # [W /K m],
     }

BOUNDARY_INPUT_PARAMETERS = \
    {
        'air temp': 35.0,  # [Celsius]
        'plate temp': 35.0,  # [Celsius]
        'air rh': 0.4,  # [fraction from 0 to 1]
        'plate rh': 1.0,  # [fraction from 0 to 1]
    }

FABRIC_IC_INPUT = \
    {
        'initial clothing temp': 33.7,  # [Celsius]
        'initial clothing rh': 0.35,  # [fraction from 0 to 1]
    }  # assumes iso-humid and iso-thermo

PDE_PHYSICS_INPUT = \
    {
        'membrane_air_length_equiv': 5.0E-3,  # [m]
        # 5mm effective air layer thickness between plate and textile at plate.
        # water-vapor permeable membrane has an inherent vapor resistance that’s always there even during calibration
        'length_still_air': 3.9E-3,  # [m]
        # 3.9mm effective air layer thickness between outer textile and ambient air

        'eps_clothing': 0.95,
        'sigma': 5.67 * 10 ** (-8),  # [w /m^2 K]
        'h_convection': 7.75,  # [W/m^2 K]
        # Heat transfer coefficient of air convection assuming boundary air-layer between outer clothing and environment
    }
# TOLERANCE = 0
# RELATIVE_TOLERANCE = 1E-5
# ABSOLUTE_TOLERANCE = 1E-7
NUMBER_OF_NODES = 30

