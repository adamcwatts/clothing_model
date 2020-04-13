# currently modeling wool with parameters from Neves et al. 2015 paper

TIME_INPUT_PARAMETERS = \
    {
        'start': 0,  # [time in seconds]
        'finish [min]': 10.0,  # [minutes]
    }

FABRIC_INPUT_PARAMETERS = \
    {
     'fabric: density effective [kg/m^3]': 103.6,  # [kg/m^3]
     'fabric: thickness [m]': 0.00857,  # [m]
     'fabric: R_ef [m]': 23.209,  # [m^2 Pa / W]
     'fabric: radiation coefficient [-]': 0.95,  # radiation coefficient
     'fabric: initial temp [C]': 33.7,  # [Celsius]
     'fabric: initial RH [-]': 0.35,  # [fraction from 0 to 1]
     }

FIBER_INPUT_PARAMETERS = {
    'fiber: regain[-]': 0.15,  # [fraction from 0 to 1]
    'fiber: density dry [kg/m^3]': 1300,  # [kg/m^3]
    'fiber: specific heat [J/ Kg K]': 1360,  # [J/ Kg K ]
    'fiber: diffusion to diameter [1/s]': 4.88e-4,  # [1 / s]
    'fiber: thermal conductivity [W/ K m]': 0.20,  # [W /K m],
}

WATER_INPUT_PARAMETERS = {
    'water: density [Kg/ m^3]': 1000,  # [Kg/m^3]
    'water: specific heat [J/ Kg K]':  4179,  # [J / (kg K)]
    'water: thermal conductivity [W/ K m]': 0.613,
    'water: molecular weight [g/mol]': 18.01528,  # [g / mol]
    # T. Bergman and A. Lavine, Fundamentals of Heat and Mass Transfer, 8th ed. Wiley.
}

WATER_VAPOR_INPUT_PARAMETERS = {
    'water vapor: specific heat [J/ Kg K]': 1862,
    'water vapor: thermal conductivity [W/ K m]': 2.46e-2,
}

AIR_INPUT_PARAMETERS = {
    'air: dry density [Kg/m^3]': 1.161,  # [Kg/m^3] at 35 C
    'air: specific heat [J/ Kg K]': 1003,
    'air: thermal conductivity [W/ K m]': 2.56e-2,
    'air: still air length over fabric [m]': 3.9E-3,  # [m]
    'air: molecular weight [g/mol]': 28.97,
    # 3.9mm effective air layer thickness between outer textile and ambient air
}

BOUNDARY_INPUT_PARAMETERS = \
    {
        'air temp': 35.0,  # [Celsius]
        'plate temp': 35.0,  # [Celsius]
        'air rh': 0.4,  # [fraction from 0 to 1]
        'plate rh': 1.0,  # [fraction from 0 to 1]
    }


PDE_PHYSICS_INPUT = \
    {
        'membrane_air_length_equiv': 5.0E-3,  # [m]
        # 5mm effective air layer thickness between plate and textile at plate.
        # water-vapor permeable membrane has an inherent vapor resistance that’s always there even during calibration
        'sigma': 5.67 * 10 ** (-8),  # [w /m^2 K]
        'h_convection': 7.75,  # [W/m^2 K]
        # Heat transfer coefficient of air convection assuming boundary air-layer between outer clothing and environment
        'R [J/mol K]': 8.3145,  # [J/mol K]
        'H_VAPORIZATION [J/g]':  2418,  # J / g
    }

NUMBER_OF_NODES = 15
