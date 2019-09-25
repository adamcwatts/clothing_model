TIME_INPUT_PARAMETERS = \
    {
        'start': 0,  # [time in seconds]
        'discrete time step': 0.1,  # [seconds]
        'finish min': 1,  # [minutes]
    }

FABRIC_INPUT_PARAMETERS = \
    {'regain': 0.15,  # [fraction from 0 to 1]
     '% porosity of air in fabric': 0.931,  # [fraction from 0 to 1]
     'dry fiber density': 1300,  # [kg/m^3]
     'fiber specific heat': 1360,  # [J/ Kg K ]
     'fabric thickness': 0.00857,  # [m]
     'R_ef': 23.4  # [m^2 Pa / W]
     }

BOUNDARY_INPUT_PARAMETERS = \
    {
        'air temp': 35,  # [Celsius]
        'plate temp': 35,  # [Celsius]
        'air rh': 0.4,  # [fraction from 0 to 1]
        'plate rh': 1,  # [fraction from 0 to 1]
    }

IC_INPUT_PARAMETERS = \
    {
        'initial clothing temp': 33.7,  # [Celsius]
        'initial clothing rh': 0.35,  # [fraction from 0 to 1]
    }


TOLERANCE = 0
NUMBER_OF_NODES = 21
