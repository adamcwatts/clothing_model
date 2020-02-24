import functions
import MODEL_BC_IC as BC_IC  # BC_IC - Boundary and Initial Conditions
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

PI = BC_IC.PDE_PHYSICS_INPUT
H_VAPORIZATION = 2418  # J / g


def convert_float_to_array(dictionary) -> 'dictionary with numpy array values':
    for key, value in dictionary.items():
        dictionary[key] = np.array([value])
    return dictionary


def fabric_initial_conditions_to_df(fabric_ic: 'temp and relative humidity',
                                    fabric_data_dict: 'fabric data dictionary',
                                    node_thickness: 'array of thickness of nodes') -> 'Pandas array':
    # assumes fabric begins in an isothermal and iso-humid state
    temp_dict = {}
    for key, value in fabric_ic.items():
        # value = value[0]  # turns 1D array into scalar :TODO might not need this
        multiplier = np.ones(node_count)  # multiplies values by node count

        if key == 'initial clothing temp':
            temp_dict['initial clothing temp [C]'] = value * multiplier
            temp_dict['initial clothing temp [K]'] = (value + 273.15) * multiplier
        else:
            temp_dict[key] = value * multiplier

    for key, value in fabric_data_dict.items():
        if ('fiber' or 'thickness') not in key:
            temp_dict[key] = value * multiplier

    temp_dict['R_ef [m^2 Pa / W]'] = temp_dict.pop('R_ef')  # adds units
    temp_dict['total fabric thickness [m]'] = temp_dict.pop('fabric thickness')

    fabric_water_concentration = fabric_data_dict['dry fiber density'] * fabric_data_dict[
        'regain'] * functions.absorption(np.array(fabric_ic['initial clothing rh'])) * fabric_data_dict[
                                     'fabric porosity']

    temp_dict['water concentration absorbed by fabric [kg/m^3]'] = fabric_water_concentration * multiplier

    water_concentration_in_trapped_air_and_fabric = functions.concentration_calc(fabric_ic['initial clothing temp'],
                                                                                 fabric_ic['initial clothing rh'])
    temp_dict['water concentration in air [g/m^3]'] = water_concentration_in_trapped_air_and_fabric * multiplier

    temp_dict['thickness [m]'] = node_thickness

    df = pd.DataFrame.from_dict(temp_dict)
    names = [f'Fabric Node {x}' for x in range(node_count)]
    df.index = names
    # df.rename(columns={'initial clothing rh': 'Initial clothing rh'}, inplace=True)
    return df


def boundary_conditions_with_kelvin(boundary_condition):
    updated_boundary_conditions = boundary_condition.copy()
    for key, value in boundary_condition.items():
        if 'temp' in key:
            new_key_celsius = key + ' [C]'
            new_key_kelvin = key + ' [K]'
            updated_boundary_conditions[new_key_celsius] = updated_boundary_conditions.pop(key)  # update key
            updated_boundary_conditions[new_key_kelvin] = value + 273.15  # kelvin
            if 'air' in key:
                updated_boundary_conditions['water concentration in ambient air [g/m^3]'] = \
                    functions.concentration_calc(boundary_condition[key], boundary_condition['air rh'])[0]
            if 'plate' in key:
                updated_boundary_conditions['water concentration at plate [g/m^3]'] = \
                    functions.concentration_calc(boundary_condition[key], boundary_condition['plate rh'])[0]

    # df = pd.DataFrame.from_dict(boundary_condition)
    return updated_boundary_conditions


def calculate_time_step():
    # fabric -> material properties
    # assumes dt is global even though alpha changes after each time step. Alpha is based off dry fabric property
    # larger alpha -> smaller time step . Air_alpha >> water_alpha or fiber_alpha

    dx = BC_IC.FABRIC_INPUT_PARAMETERS['fabric thickness'] / node_count  # meters dx size
    k_fiber = BC_IC.FABRIC_INPUT_PARAMETERS['thermal conductivity of fiber']
    p_fiber = BC_IC.FABRIC_INPUT_PARAMETERS['dry fiber density']
    c_fiber = BC_IC.FABRIC_INPUT_PARAMETERS['fiber specific heat']

    p_fab = (0.07) * (p_fiber) + (1 - 0.07) * (1.225)  # density  averaged [kg/m^3]
    c_fab = (0.07) * (c_fiber) + (1 - 0.07) * (1.000)  # specific heat capacity averaged [J/kg k]
    k_fab = (0.07) * (k_fiber) + (1 - 0.07) * (0.024)  # thermal conductivity averaged [w/m K]

    alpha = k_fab / (p_fab * c_fab)
    # water alpha = 1.43 × 10−7
    dt = dx ** 2 / (2 * alpha)  # less than this number
    return dt


def ode_plateflux_temp(t, y):
    #  y[0] refers to heat flux needed by plate
    #  y[n] refers to temperature at nth node : 1 <= n <= number of nodes

    # ode_params = {}
    # n numbers of nodes + 1 node for boundary at plate

    temp_at_clothing_boundary = y[PDE.node_count]  # last element of vector associated with temperatures

    h_radiation = PI['eps_clothing'] * PI['sigma'] * (
            temp_at_clothing_boundary + boundary_conditions['air temp [K]']) * (
                          temp_at_clothing_boundary ** 2 + boundary_conditions['air temp [K]'] ** 2)

    # r_out = 1 / (h_radiation + OPI['h_convection'])    U_out = 1 / r_out

    r_rad_convect = 1 / (h_radiation + PI['h_convection'])

    # delta_x = df_fabric_initial_conditions['thickness [m]'].to_numpy()

    K = df_wet_fabric['wet fabric thermal conductivity [W/mk]'].to_numpy()
    C = df_wet_fabric['wet fabric specific heat [J/kg K]'].to_numpy()
    P = df_wet_fabric['wet fabric density [kg/m^3]'].to_numpy()
    U = K / PDE.wool_data.delta_x

    # corrects last element for radiation, convection heat loss
    U[-1] = 1 / (r_rad_convect + 0.5 * (PDE.wool_data.delta_x[-1] / K[-1]))

    material = 1 / (PDE.wool_data.delta_x * np.multiply(C, P))

    temperature_array = np.append(y[1:], boundary_conditions['air temp [K]'])

    q_dry = np.diff(temperature_array) * -1 * U  # -1 due to diff goes a[i+1] - a[i], I need a[i] - a[i+1]
    q_dry = np.concatenate(([y[0]], q_dry), axis=0)  # adds heat flow from plate

    # flux at plate ODE
    q_gen = (boundary_conditions['plate temp [K]'] - y[1]) * (K[0] / (PDE.wool_data.delta_x[0] * 0.5))
    heat_flows = PDE.q_sorp_history[idx, :] - PDE.q_evap_history[idx, :] + PDE.q_conden_history[idx, :]

    # TEMPORARY NO OTHER HEAT FLOWS COUPLED TO TEST
    clothing_dydt = material * (np.diff(q_dry * -1) + heat_flows)
    # -1 due to diff goes a[i+1] - a[i], I need a[i] - a[i+1]

    dydt = np.concatenate(([q_gen], clothing_dydt), axis=0)

    # q_1_2 = (y[1] - y[2]) * U[0]
    # q_2_3 = (y[2] - y[3]) * U[1]
    # q_3_out = (y[3] - boundary_conditions['air temp [K]']) * U[2]  #
    #
    # dydt_1 = material[0] * (y[0] - q_1_2)  # delta temp at node ODE 1
    # dydt_2 = material[1] * (q_1_2 - q_2_3) # delta temp at node ODE 2
    # dydt3 = material[2] * (q_2_3 - q_3_out) # delta temp at node ODE 3

    # dydt = [q_gen, dydt_1, dydt_2, dydt3]

    return dydt


def ode_concentration(t, c):
    # delta_x = df_fabric_initial_conditions['thickness [m]'].to_numpy()  # thickness of each node
    # D_WF = df_fabric_initial_conditions['diffusivity of water though fabric [m^2 /s]'].to_numpy()
    #
    # delta_path = np.copy(delta_x)  # copy array to create actual path between nodes
    # delta_path[0] = delta_path[0] * 0.5 + OPI['membrane_air_length_equiv']  # node next to hot plate
    # delta_path[-1] = delta_path[-1] * 0.5 + OPI['length_still_air']
    #
    # material = D_WF / delta_x

    concentration_array = np.concatenate(
        ([boundary_conditions['water concentration at plate [g/m^3]']], c,
         [boundary_conditions['water concentration in ambient air [g/m^3]']]), axis=0)

    delta_c = np.diff(concentration_array) * -1

    # dcdt_0 = material[0] * (delta_c[0] / (delta_x[0] * 0.5 + membrane_air_length_equiv) - delta_c[1] / delta_x[0])
    # dcdt_1 = material[1] * (delta_c[1] / delta_x[1] - (delta_c[2]) / delta_x[1])
    # dcdt_2 = material[2] * (delta_c[2] / delta_x[2] - delta_c[3] / (delta_x[2] * 0.5 + still_air_length))

    dcdt = PDE.wool_data.concen_stiffness * \
           (delta_c[0:-1] / PDE.wool_data.delta_path[0:-1] - delta_c[1:] / PDE.wool_data.delta_path[1:])

    return dcdt


def solution_to_df(ode_class):
    ode_class.Solutions.temp_c_history = ode_class.Solutions.temp_k_history - 273.15  # append a C array
    ode_class.Solutions.rh_history *= 100  # convert fractions to  %
    ode_class.Solutions.rh_post_absorption_history *= 100  # convert fractions to  %

    # takes ODE solution and creates data frame that's easy to read solutions off of

    all_keys = ode_class.Solutions.__dict__.keys()
    export_keys = [key for key in all_keys if 'history' in key]  # return only keys with 'history' in them
    stripped_keys = [key.split('_history')[0] for key in export_keys]  # remove _history term at end of string
    cleaned_key = [key.replace("_", ' ') for key in stripped_keys if 'q' not in key]  # replace underscore but no q
    cleaned_key.extend(key for key in stripped_keys if 'q' in key)  # extend with q while keeping underscore

    unit_dict = {
        'concen': '[g/m^3]', 'temp k': '[K]', 'temp c': '[C]', 'q': '[W/m^2]', 'enthalpy': '[J/g]',
        'heat': '[W/m^2]', 'fabric': '[g/m^3]', 'rh': '[%]'
    }

    index_names = [f'{round(k, 5)} [s]' for k in tspan]
    column_names = []  # nested list of column names for each dataframe to export

    for key in cleaned_key:
        for sub_key in list(unit_dict.keys()):
            if sub_key in key:
                # print(key, unit_dict[sub_key])

                names = [f'{key.title()} at Node {h} {unit_dict[sub_key]}' for h in range(0, node_count)]
                names.insert(0, 'Time [minutes]')
                column_names.append(names)
    print('\n')

    for k, item in enumerate(export_keys):
        print(f'{k}:', f'Exporting {item}.csv file')
        # add column of minutes at the beginning of each data array
        my_data = np.concatenate((np.round(tspan[:, None], decimals=3) / 60, getattr(ode_class.Solutions, item)),
                                 axis=1)
        df = pd.DataFrame(data=my_data, columns=column_names[k], index=index_names)

        # TODO write function that checks to see if solutions already exist in folder,
        #  if so write to a different folder with a parameter of user choice
        file_path = './SOLUTION_CSV/' + f'{item}.csv'
        df.to_csv(file_path, encoding='utf-8')

    display_messages(2)
    # fabric_solutions = {f'Temp at Node {i + 1} [C]': ode_sol.y[i + 1] - 273.15 for i in range(ODE.node_count)}
    # fabric_solutions['Time (Minutes)'] = ode_sol.t / 60
    # fabric_solutions['Hot Plate Flux'] = ode_sol.y[0]


def display_messages(case):
    if case == 1:
        print('\n****************************** CODE STARTING  *******************************')
    elif case == 2:
        print('\n****************************** CODE COMPLETED  *******************************')


def print_progressbar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', print_end="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=print_end)
    # Print New Line on Complete
    if iteration == total:
        print()


class MaterialData:
    def __init__(self):
        self.delta_x = None
        self.delta_path = None
        self.D_WF = None
        self.concen_stiffness = None


class ArrayGen:
    def __init__(self, dict):
        for key, value in dict.items():
            setattr(self, key, np.zeros(value))


class SolutionTree:
    def __init__(self, my_pde, node_tuple):
        zeros_array = np.zeros(node_tuple)
        time_size = node_tuple[0]

        self.temp_k_history = np.zeros(node_tuple)  # [Kelvin]
        self.temp_c_history = np.zeros(node_tuple)  # [Celsius]
        self.temp_k_history[0, :] = my_pde.fabric_IC['initial clothing temp [K]'].to_numpy()

        self.raw_water_concen_history = np.copy(zeros_array)  # [g/m^3]
        self.raw_water_concen_history[0, :] = my_pde.fabric_IC[
            'water concentration in air [g/m^3]']

        self.water_vapor_concen_history = np.copy(zeros_array)  # [g/m^3]
        self.water_vapor_concen_history[0, :] = my_pde.fabric_IC[
            'water concentration in air [g/m^3]']

        self.conden_concen_history = np.copy(zeros_array)  # [(g/m^3]  NEED TO VERIFY
        self.fabric_sorption_history = np.copy(
            zeros_array)  # [g/m^3] :  water density absorbed or desorbed by fabric

        self.fiber_water_concen_history = np.copy(zeros_array)  # [(kg/m^3]
        self.fiber_water_concen_history[0, :] = my_pde.fabric_IC[
            'water concentration absorbed by fabric [kg/m^3]']

        self.rh_history = np.copy(zeros_array)
        self.rh_history[0, :] = my_pde.fabric_IC['initial clothing rh']

        self.rh_post_absorption_history = np.copy(zeros_array)
        self.rh_post_absorption_history[0, :] = my_pde.fabric_IC['initial clothing rh']

        self.enthalpy_sorption_history = np.copy(zeros_array)

        self.q_conden_history = np.copy(zeros_array)
        self.q_evap_history = np.copy(zeros_array)

        self.q_sorp_history = np.copy(zeros_array)  # sum of q_1 and q_L
        self.q_1_history = np.copy(zeros_array)
        self.q_L_history = np.copy(zeros_array)

        self.heat_flows = np.zeros((time_size, 1 + my_pde.node_count))

        return


class SolutionSetup:
    def __init__(self, fabric_ic, fabric_df_wet, allocated_time_size, dt, BC_kelvin):
        self.node_count = BC_IC.NUMBER_OF_NODES
        self.fabric_wet_df = fabric_df_wet  # Imports wet fabric dataframe, updates each time step
        self.fabric_IC = fabric_ic  # Imports fabric initial conditions dataframe,
        self.dt = dt
        self.dx = self.fabric_IC['thickness [m]'].values[1]
        self.BC_K = BC_kelvin  # imports kelvin BC conditions
        node_tuple = (allocated_time_size, self.node_count)

        self.Solutions = SolutionTree(self, node_tuple)

        self.wool_data = MaterialData()

    def zeta_gen(self):
        dx_array = np.zeros(self.node_count)  # Create array of zeros for dx_array
        dx_array.fill(self.dx)  # fill dx array with interior node which is dx_array

        zeta_i = self.dt / (self.fabric_wet_df['wet fabric density [kg/m^3]'].values *
                            dx_array ** 2 * self.fabric_wet_df['wet fabric specific heat [J/kg K]'].values)

        zeta_i[-1] = zeta_i[-1] * 2 * self.dx  # adjust last value for Mixed BC

        return zeta_i

    def radiation_transfer_coeff(self, t=0):
        # t is the current time step

        h = PI['eps_clothing'] * PI['sigma'] * (self.Solutions.temp_k_history[t, -1] + self.BC_K['air temp [K]']) * \
            (self.Solutions.temp_k_history[t, -1] ** 2 + self.BC_K['air temp [K]'] ** 2)
        return h

    def heat_matrix_ftcs_gen(self, t=0):
        zeta = self.zeta_gen()  # run zeta_gen to produce zeta array
        h_r = self.radiation_transfer_coeff(t)
        h_c = PI['h_convection']

        n = self.node_count + 1  # plus 1 for ambient boundary temperature
        a = np.zeros([n, n])
        k = self.fabric_wet_df['wet fabric thermal conductivity [W/mk]']

        # returns row vector of total heat generated at node from previous time step
        q_net = self.Solutions.q_conden_history[t, :] + self.Solutions.q_evap_history[t, :] + \
                self.Solutions.q_1_history[t, :] + self.Solutions.q_L_history[t, :]

        a[0, 0] = 1  # assumes dirichlet BC at far left Node
        a[n - 1, n - 1] = 1

        for row in range(n - 1):
            if row != 0 and row != n - 2:  # dont modify 0th row. Nor modify last or 2nd to last row
                # A[row, row-1:row+2] = fo
                a[row, row - 1], a[row, row + 1] = zeta[row] * k[row - 1], zeta[row] * k[row + 1]
                # a[row, row] = fo[row]
                a[row, row] = 1 - zeta[row] * (k[row - 1] + k[row + 1])

            elif row != 0 and row == n - 2:  # second to last row modify for radiation and convection transfer
                a[row, row - 1], a[row, row + 1] = zeta[-1] * (k[row - 1] / self.dx), zeta[-1] * (h_r + h_c)
                a[row, row] = 1 - zeta[-1] * (h_r + h_c + k[row - 1] / self.dx)

        q_net[0:-1] *= self.dx  # dont multiply last element with dx as zeta corrected for that
        return a, q_net*zeta

    def FTCS_heat(self, ft_matrix, phi, t=0):
        temp_old = np.append(self.Solutions.temp_k_history[t, :], self.BC_K['air temp [K]'])  # attach BC to array
        temp_new = np.matmul(ft_matrix, temp_old)

        phi[0] = 0  # boundary condition cannot heat up!
        return temp_new[0:-1] + phi  # dont return ambient temp hence 0 to -1

    def gamma_gen(self):

        dx_array = np.zeros(self.node_count)  # Create array of zeros for dx_array
        dx_array.fill(self.dx)  # fill dx array with interior node which is dx_array

        gamma_i = self.dt / (self.fabric_IC['diffusion resistance through fabric [s/m]'].values * dx_array)
        #
        gamma_i[-1] = gamma_i[-1] * self.dx  # adjust last value for Mixed BC

        return gamma_i

    def diffusion_matrix_ftcs_gen(self):
        gamma = self.gamma_gen()  # run zeta_gen to produce zeta array
        n = self.node_count + 1  # plus 1 for ambient boundary temperature
        a = np.zeros([n, n])
        a[0, 0] = 1  # TODO might want to preload these arrays so it doesnt need to do this every time step
        a[n - 1, n - 1] = 1

        for row in range(n - 1):
            if row != 0 and row != n - 2:  # dont modify 0th row. Nor modify last or 2nd to last row
                # A[row, row-1:row+2] = fo
                a[row, row - 1], a[row, row + 1] = gamma[row], gamma[row]
                # a[row, row] = fo[row]
                a[row, row] = 1 - 2 * gamma[row]

            elif row != 0 and row == n - 2:  # second to last row modify for radiation and convection transfer
                a[row, row - 1], a[row, row + 1] = gamma[-1] / self.dx, gamma[-1] / PI['length_still_air']
                a[row, row] = 1 - gamma[-1] * (1 / PI['length_still_air'] + 1 / self.dx)

        return a

    def FTCS_diffusion(self, ft_matrix, t=0):
        concentration_old = np.append(self.Solutions.raw_water_concen_history[t, :],
                                      self.BC_K['water concentration in ambient air [g/m^3]'])  # attach BC to array
        concentration_new = np.matmul(ft_matrix, concentration_old)

        return concentration_new[0:-1]  # dont return ambient air concentration

    def concentration_constants(self):
        # Constants for ODE concentration
        self.wool_data.delta_x = self.fabric_IC['thickness [m]'].to_numpy()  # thickness of each node
        self.wool_data.D_WF = self.fabric_IC['diffusivity of water though fabric [m^2 /s]'].to_numpy()

        self.wool_data.delta_path = np.copy(
            self.wool_data.delta_x)  # copy array to create actual path between nodes

        # node next to hot plate
        self.wool_data.delta_path[0] = self.wool_data.delta_path[0] * 0.5 + PI['membrane_air_length_equiv']

        self.wool_data.delta_path = np.append(self.wool_data.delta_path,
                                              self.wool_data.delta_path[-1] * 0.5 + PI['length_still_air'])

        self.wool_data.concen_stiffness = self.wool_data.D_WF / self.wool_data.delta_x

    def post_concentration(self, index):
        # REMEMBER INDEX ALREADY IS PLUS 1

        # under normal conditions where RH of air > RH fabric
        # fiber sorption occurs, reduces ambient RH while absorption water, giving off heat
        # else, inverse occurs

        # if near 100% ambient RH, hard to tell if sorption occurs first then condensation or inverse

        current_concentration = self.Solutions.raw_water_concen_history[index, :]

        # RH before any sorption of water in fabric
        self.Solutions.rh_history[index, :] = \
            functions.relative_humidity_calc(current_concentration, self.fabric_IC['initial clothing temp [K]'])

        # create passing tuple to feed into rh_equilibrium function
        passing_parameters = (
            self.fabric_IC, current_concentration, self.Solutions.temp_k_history[index - 1, :],
            self.Solutions.rh_post_absorption_history[index - 1, :])

        # Use Relative humidity equilibrium function
        relative_humidity_equilibrium, sorption, eq_air_concentration = functions.rh_equilibrium(*passing_parameters)

        # assumes sorption priority before condensation
        self.Solutions.rh_post_absorption_history[index, :] = relative_humidity_equilibrium

        # sorption
        self.Solutions.fabric_sorption_history[index, :] = sorption  # [g/m^3]

        # Updated water concentration inside the fabric [kg/m^3], scale sorption from grams to kg   / m^3
        self.Solutions.fiber_water_concen_history[index, :] = self.Solutions.fiber_water_concen_history[index - 1, :] + \
                                                              (sorption / 1000)

        # equilibrium air concentration due to equilibrium
        self.Solutions.water_vapor_concen_history[index, :] = eq_air_concentration

        # create passing tuple to feed into condensation_checker
        passing_parameters = (self.Solutions.rh_post_absorption_history[index, :],
                              self.Solutions.water_vapor_concen_history[index, :],
                              self.Solutions.temp_k_history[index - 1, :],
                              self.Solutions.conden_concen_history[index, :])

        # Use condensation_checker to check for condensation using previous time steps temperatures
        self.Solutions.rh_post_absorption_history[index, :], \
        self.Solutions.water_vapor_concen_history[index, :], \
        self.Solutions.conden_concen_history[index, :] = functions.condensation_checker(*passing_parameters)

    def post_concentration_sorption_flows(self, index):
        # REMEMBER INDEX ALREADY IS PLUS 1

        node_thickness = self.fabric_IC['thickness [m]']

        # heat of sorption due to regain in [J/g]
        self.Solutions.enthalpy_sorption_history[index, :] = functions.vectorized_h_sorp_cal(
            self.Solutions.rh_post_absorption_history[index - 1, :],
            self.Solutions.rh_post_absorption_history[index, :])

        flux = self.Solutions.fabric_sorption_history[index, :] * node_thickness / self.dt

        self.Solutions.q_1_history[index, :] = flux * self.Solutions.enthalpy_sorption_history[index, :]
        self.Solutions.q_L_history[index, :] = flux * H_VAPORIZATION
        self.Solutions.q_sorp_history[index, :] = self.Solutions.q_1_history[index, :] + self.Solutions.q_L_history[
                                                                                         index, :]

        self.Solutions.q_conden_history[index, :] = functions.q_condensation(
            self.Solutions.conden_concen_history[index, :],
            node_thickness, self.dt)

        self.Solutions.q_evap_history[index, :], self.Solutions.conden_concen_history[index, :] = \
            functions.q_evaporation(self.Solutions.rh_post_absorption_history[index, :],
                                    self.Solutions.conden_concen_history[index, :],
                                    node_thickness, self.dt)

    def post_plateflux_temp(self, ode_temp_solution, index):
        current_temperature = ode_temp_solution.y[1:, -1]
        self.temp_k_history[index, :] = current_temperature

        self.heat_flows[index, 0] = ode_temp_solution.y[0, -1]


if __name__ == '__main__':
    node_count = BC_IC.NUMBER_OF_NODES
    start = BC_IC.TIME_INPUT_PARAMETERS['start']
    finish_seconds = BC_IC.TIME_INPUT_PARAMETERS['finish min'] * 60

    dt_min = calculate_time_step()
    time_step_lower_bound = functions.dt_generator(dt_min)  # base 10 time step
    time_step = time_step_lower_bound * (dt_min // time_step_lower_bound)  # even rounds to lower decimal place

    tspan = np.arange(start, finish_seconds + time_step, time_step)
    t_size = tspan.shape[0]

    fabric_data = functions.fabric_parameters(BC_IC.FABRIC_INPUT_PARAMETERS)
    fabric_node_dimensions = functions.node_generator(fabric_data, node_count)

    boundary_conditions = boundary_conditions_with_kelvin(BC_IC.BOUNDARY_INPUT_PARAMETERS)
    fix_floats_to_numpy_array = convert_float_to_array(BC_IC.FABRIC_IC_INPUT)

    df_fabric_initial_conditions = fabric_initial_conditions_to_df(fix_floats_to_numpy_array, fabric_data,
                                                                   fabric_node_dimensions)

    df_wet_fabric = functions.wet_fabric_calc(df_fabric_initial_conditions,
                                              df_fabric_initial_conditions['initial clothing rh'].to_numpy())

    # initialize PDE CLASS arrays and key parameters
    PDE = SolutionSetup(df_fabric_initial_conditions, df_wet_fabric, t_size, time_step, boundary_conditions)

    # First time time through SCHEME
    A, Phi = PDE.heat_matrix_ftcs_gen()
    first_temp = np.append(PDE.Solutions.temp_k_history[0, :], PDE.BC_K['air temp [K]'])
    first_temp[0] = PDE.BC_K['plate temp [K]']
    Phi[0] = 0
    PDE.Solutions.temp_k_history[0 + 1, :] = np.matmul(A, first_temp)[0:-1] + Phi  # dont include last value, since its T_ambient

    B = PDE.diffusion_matrix_ftcs_gen()
    first_concen = np.append(PDE.Solutions.raw_water_concen_history[0, :],
                             PDE.BC_K['water concentration in ambient air [g/m^3]'])
    first_concen[0] = PDE.BC_K['water concentration at plate [g/m^3]']
    PDE.Solutions.raw_water_concen_history[0 + 1, :] = np.matmul(B, first_concen)[0: -1]  # dont include C_ambient
    PDE.post_concentration(0 + 1)
    PDE.post_concentration_sorption_flows(0 + 1)
    PDE.fabric_wet_df = functions.wet_fabric_calc(df_fabric_initial_conditions,
                                                  PDE.Solutions.rh_post_absorption_history[0 + 1, :])

    end_value = t_size - 1
    for i in range(1, end_value):
        print_progressbar(i, end_value, prefix='Progress', suffix='Complete', length=50)
        A, Phi = PDE.heat_matrix_ftcs_gen(t=i)
        PDE.Solutions.temp_k_history[i + 1, :] = PDE.FTCS_heat(A, Phi, t=i)

        B = PDE.diffusion_matrix_ftcs_gen()
        PDE.Solutions.raw_water_concen_history[i + 1, :] = PDE.FTCS_diffusion(B, t=i)

        PDE.post_concentration(i + 1)
        PDE.post_concentration_sorption_flows(i + 1)
        PDE.fabric_wet_df = functions.wet_fabric_calc(df_fabric_initial_conditions,
                                                      PDE.Solutions.rh_post_absorption_history[i + 1, :])

        # update fabric properties due to regain
        # df_wet_fabric = functions.wet_fabric_calc(df_fabric_initial_conditions, PDE.rh_post_absorption_history[1, :])

    print('DONE')
    # PDE.concentration_constants(df_fabric_initial_conditions)  # pre-calculate constants for concentration ODE
    #
    # # ODE STUFF BELOW
    # plate_flux_temp_init = np.insert(df_fabric_initial_conditions['initial clothing temp [K]'].to_numpy(), 0, 0)
    # concentration_init = df_fabric_initial_conditions['water concentration in air [g/m^3]'].to_numpy()
    #
    # idx = 1
    # end_value = tspan.shape[0] - 1
    # display_messages(1)
    # print_progressbar(0, end_value, prefix='Progress', suffix='Complete', length=50)
    #
    # # start concentration ODE
    # concentration_solution = solve_ivp(ode_concentration, [tspan[0], tspan[idx]], concentration_init, rtol=RTOL,
    #                                    atol=ATOL)
    # PDE.post_concentration(concentration_solution, idx)
    # PDE.post_concentration_sorption_flows(df_fabric_initial_conditions, idx)
    #
    # # update fabric properties due to regain
    # df_wet_fabric = functions.wet_fabric_calc(df_fabric_initial_conditions, PDE.rh_post_absorption_history[1, :])
    #
    # # start temp/flux ODE
    # PDE.matrix_ftcs_gen(df_wet_fabric, node_count)
    #
    # plate_flux_temp_solution = solve_ivp(ode_plateflux_temp, [tspan[0], tspan[idx]], plate_flux_temp_init,
    #                                      rtol=RTOL,
    #                                      atol=ATOL)
    # PDE.post_plateflux_temp(plate_flux_temp_solution, idx)
    #
    # for i in range(1, end_value):
    #     # display_messages(3)
    #     idx = i
    #     print_progressbar(i, end_value, prefix='Progress', suffix='Complete', length=50)
    #
    #     c_sol = solve_ivp(ode_concentration, [tspan[i], tspan[i + 1]], PDE.water_vapor_concen_history[i, :],
    #                       rtol=RTOL,
    #                       atol=ATOL)
    #     PDE.post_concentration(c_sol, i + 1)
    #     PDE.post_concentration_sorption_flows(c_sol, i + 1)
    #
    #     # is it  ODE.rh_post_absorption_history[i, :] or  ODE.rh_post_absorption_history[i + 1, :] ?????????
    #     df_wet_fabric = functions.wet_fabric_calc(df_fabric_initial_conditions,
    #                                               PDE.rh_post_absorption_history[i, :])
    #
    #     new_ic = np.insert(PDE.temp_k_history[i, :], 0, PDE.heat_flows[i, 0])
    #
    #     t_sol = solve_ivp(ode_plateflux_temp, [tspan[i], tspan[i + 1]], new_ic, rtol=RTOL, atol=ATOL)
    #     PDE.post_plateflux_temp(t_sol, i + 1)
    #
    solution_to_df(PDE)

    # concentration_solution = solve_ivp(ode_concentration, [tspan[0], tspan[-1]], concentration_init, t_eval=tspan)
    # concentration_sol_df = solution_to_df(concentration_solution, 'c')
    # concentration_sol_df.to_csv('concentration_solution_summary.csv', sep='\t', encoding='utf-8')
    #
    # plate_flux_temp_solution = solve_ivp(ode_plateflux_temp, [tspan[0], tspan[-1]], plate_flux_temp_init, t_eval=tspan)
    # temperature_sol_df = solution_to_df(plate_flux_temp_solution, 't')
    # temperature_sol_df.to_csv('plate_flux_and_temp_solution_summary.csv', sep='\t', encoding='utf-8')

    # solution_plots(sol_df)
