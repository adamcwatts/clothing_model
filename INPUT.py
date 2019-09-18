TIME_INPUT_PARAMETERS = \
    {
        'start': 0,
        'discrete time step': 0.1,
        'finish_min': 1,
    }

FABRIC_INPUT_PARAMETERS = \
    {'regain': 0.15,  # %
     'air porosity': 0.931,  # %
     'dry density': 1300,  # [kg/m^3]
     'specific heat': 1360,  # [J/ Kg K ]
     'fabric thickness': 0.00857,  # [m]
     'R_ef': 23.4  # [m^2 Pa / W]
     }

BOUNDARY_INPUT_PARAMETERS = \
    {
        'air temp': 35,  # [Celsius]
        'plate temp': 35,  # [Celsius]
        'air rh': 0.4,
        'plate rh': 1,
    }

TOLERANCE = 0
NUMBER_OF_NODES = 21
