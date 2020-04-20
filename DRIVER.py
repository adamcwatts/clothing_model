import functions
import MODEL_BC_IC as BC_IC  # BC_IC - Boundary and Initial Conditions
import math
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

PI = BC_IC.PDE_PHYSICS_INPUT
H_VAPORIZATION = PI['H_VAPORIZATION [J/g]']


def boundary_conditions_with_kelvin(boundary_condition):
    updated_boundary_conditions = boundary_condition.copy()
    for key, value in boundary_condition.items():
        if 'temp' in key:
            new_key_celsius = key + ' [C]'
            new_key_kelvin = key + ' [K]'
            updated_boundary_conditions[new_key_celsius] = updated_boundary_conditions.pop(key)  # update key
            updated_boundary_conditions[new_key_kelvin] = value + 273.15  # kelvin
            if 'air' in key:
                updated_boundary_conditions['water: concentration in ambient air [g/m^3]'] = \
                    functions.concentration_calc(boundary_condition[key], boundary_condition['air: rh'])[0]

                updated_boundary_conditions['water: concentration in ambient air [kg/m^3]'] = \
                    updated_boundary_conditions['water: concentration in ambient air [g/m^3]'] / 1000

            if 'plate' in key:
                updated_boundary_conditions['water: concentration at plate [g/m^3]'] = \
                    functions.concentration_calc(boundary_condition[key], boundary_condition['plate: rh'])[0]

                updated_boundary_conditions['water: concentration at plate [kg/m^3]'] = \
                    updated_boundary_conditions['water: concentration at plate [g/m^3]'] / 1000

    # df = pd.DataFrame.from_dict(boundary_condition)
    return updated_boundary_conditions


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


class ArrayGen:
    def __init__(self, dict):
        for key, value in dict.items():
            setattr(self, key, np.zeros(value))


class SolutionTree:
    def __init__(self, my_pde, node_tuple):
        zeros_array = np.zeros(node_tuple)
        time_size = node_tuple[0]

        self.temp_k_history = np.zeros(node_tuple)  # [Kelvin]
        self.temp_k_history[0, :] = my_pde.fabric_df['Temp [K]']

        self.temp_c_history = np.zeros(node_tuple)  # [Celsius]
        self.temp_c_history[0, :] = my_pde.fabric_df['Temp [C]']

        # self.raw_water_concen_history = np.copy(zeros_array)  # [g/m^3]
        # self.raw_water_concen_history[0, :] = my_pde.fabric_df[
        #     'water concentration in air [g/m^3]']
        #
        self.water_vapor_concen_history = np.copy(zeros_array)  # [kg/m^3]
        self.water_vapor_concen_history[0, :] = my_pde.materials.VAPOR['water vapor: density [Kg/m^3]']

        self.conden_concen_history = np.copy(zeros_array)  # [(kg/m^3]  NEED TO VERIFY
        self.fabric_sorption_history = np.copy(zeros_array)  # [kg/m^3] :  water density absorbed or desorbed by fabric

        # self.fiber_water_concen_history = np.copy(zeros_array)  # [(kg/m^3]
        # self.fiber_water_concen_history[0, :] = my_pde.fabric_df[
        #     'water concentration absorbed by fabric [kg/m^3]']

        self.rh_history = np.copy(zeros_array)
        self.rh_history[0, :] = my_pde.fabric_df['RH [-]']

        self.rh_post_absorption_history = np.copy(zeros_array)
        self.rh_post_absorption_history[0, :] = my_pde.fabric_df['RH [-]']

        self.enthalpy_sorption_history = np.copy(zeros_array)

        self.q_conden_history = np.copy(zeros_array)
        self.q_evap_history = np.copy(zeros_array)

        self.q_sorp_history = np.copy(zeros_array)  # sum of q_1 and q_L
        self.q_1_history = np.copy(zeros_array)
        self.q_L_history = np.copy(zeros_array)

        self.heat_flows = np.zeros((time_size, 1 + my_pde.node_count))

        self.volume_fraction_gas_history = np.copy(zeros_array)
        self.volume_fraction_gas_history[0, :] = my_pde.fabric_df["VF: gas [-]"]

        self.volume_fraction_water_history = np.copy(zeros_array)
        self.volume_fraction_water_history[0, :] = my_pde.fabric_df['VF: water [-]']
        return


class Material:
    def __init__(self, material_dict):
        # material_names = ['FABRIC', 'FIBER', 'WATER', 'AIR', 'VAPOR', 'GAS']
        # for i, material_name in enumerate(material_names):
        #     setattr(self, material_name, material_dict[i])
        self.FABRIC = material_dict[0]
        self.FIBER = material_dict[1]
        self.WATER = material_dict[2]
        self.AIR = material_dict[3]
        self.VAPOR = material_dict[4]
        self.GAS = material_dict[5]


class SolutionSetup:
    def __init__(self, material_list, grid_dimensions, allocated_time_size, dt, BC_kelvin):
        self.counter = int(0)
        self.count_to = allocated_time_size
        self.node_count = BC_IC.NUMBER_OF_NODES
        # self.fabric_wet_df = fabric_df_wet  # Imports wet fabric dataframe, updates each time step
        # self.fabric_IC = fabric_ic  # Imports fabric initial conditions dataframe,
        self.materials = Material(material_list)
        self.dt = dt
        self.dx = grid_dimensions[1]
        self.BC_K = BC_kelvin  # imports kelvin BC conditions
        node_tuple = (allocated_time_size, self.node_count)
        self.fabric_df = self.fabric_time_dependent()
        self.Solutions = SolutionTree(self, node_tuple)
        self.Coupling_Option = {'Fully Decoupled': 1, 'Coupled Equilibrium': 2, 'Coupled Double Exponential': 3, }

    @staticmethod
    def material_parameters() -> 'updated fabric dictionary':
        FABRIC = BC_IC.FABRIC_INPUT_PARAMETERS
        FIBER = BC_IC.FIBER_INPUT_PARAMETERS
        WATER = BC_IC.WATER_INPUT_PARAMETERS
        AIR = BC_IC.AIR_INPUT_PARAMETERS
        VAPOR = BC_IC.WATER_VAPOR_INPUT_PARAMETERS
        GAS = {}

        # FABRIC PROPERTIES
        plate_temp = BC_IC.BOUNDARY_INPUT_PARAMETERS['plate: temp']
        R_ef = FABRIC['fabric: R_ef [m]']
        diffusion_resistance = functions.evap_res_to_diffusion_res(plate_temp, R_ef)  # effective for fabric
        diffusivity_water_though_fabric = FABRIC['fabric: thickness [m]'] / diffusion_resistance

        FABRIC['fabric: diffusion resistance effective [s/m]'] = diffusion_resistance
        FABRIC['fabric: diffusivity effective [m^2 /s]'] = diffusivity_water_though_fabric
        rho_effective = FABRIC['fabric: density effective [kg/m^3]']
        rho_water = WATER['water: density [Kg/ m^3]']
        rho_fiber = FIBER['fiber: density dry [kg/m^3]']

        RH_initial = FABRIC['fabric: initial RH [-]']
        M_h20 = WATER['water: molecular weight [g/mol]'] * (1 / 1000)  # convert g/ mol to Kg/mol
        FABRIC['fabric: initial temp [K]'] = FABRIC['fabric: initial temp [C]'] + 273.15
        R = PI['R [J/mol K]']

        # convert kPa to Pa
        partial_vapor = RH_initial * functions.saturated_vapor_pressure(FABRIC['fabric: initial temp [C]']) * 1000

        rho_vapor = (partial_vapor * M_h20 / (R * FABRIC['fabric: initial temp [K]']))[0]

        VAPOR['water vapor: density [Kg/m^3]'] = rho_vapor  # AKA Water vapor concentration
        GAS['gas: density effective [Kg/m^3]'] = rho_vapor + AIR['air: dry density [Kg/m^3]']
        rho_gas = GAS['gas: density effective [Kg/m^3]']
        regain = FIBER['fiber: regain[-]'] * functions.absorption(FABRIC['fabric: initial RH [-]'])

        pass_args = (rho_effective, rho_water, rho_fiber, rho_gas, regain)

        eps_guess = np.array([0.9, 0.09, 0.01])
        epsilon_values = fsolve(functions.epsilon_equations, eps_guess, args=pass_args)
        epsilon_gas, epsilon_dry_fiber, epsilon_bound_water = epsilon_values

        GAS['gas: volume fraction [-]'] = epsilon_gas
        FIBER['fiber: volume fraction [-]'] = epsilon_dry_fiber
        WATER['water: volume fraction bound [-]'] = epsilon_bound_water
        GAS['gas: specific heat [J/ Kg K]'] = functions.cp_gas_calc(AIR, VAPOR, GAS)
        GAS['gas: thermal conductivity [W/ K m]'] = functions.k_gas_calc(AIR, VAPOR, GAS)

        k_solid = functions.k_solid_calc(WATER, FIBER)
        k_ef = functions.k_ef_calc(GAS, WATER, FIBER, k_solid)
        FABRIC['fabric: thermal conductivity effective [W/ K m]'] = k_ef
        C_p_eff = functions.cp_ef_calc(WATER, GAS, FIBER, AIR, VAPOR, FABRIC)
        FABRIC['fabric: specific heat effective [J/ Kg K]'] = C_p_eff

        materials = [FABRIC, FIBER, WATER, AIR, VAPOR, GAS]
        return materials

    @staticmethod
    def calculate_time_step(fabric_dict, dx):
        # fabric -> material properties
        # assumes dt is global even though alpha changes after each time step. Alpha is based off dry fabric property
        # larger alpha -> smaller time step . Air_alpha >> water_alpha or fiber_alpha

        # dx = BC_IC.FABRIC_INPUT_PARAMETERS['fabric: thickness [m]'] / BC_IC.NUMBER_OF_NODES  # meters dx size

        k_ef = fabric_dict['fabric: thermal conductivity effective [W/ K m]']
        C_p_eff = fabric_dict['fabric: specific heat effective [J/ Kg K]']
        rho_ef = fabric_dict['fabric: density effective [kg/m^3]']

        alpha = k_ef / (rho_ef * C_p_eff)
        dt = dx ** 2 / (2 * alpha)  # less than this number
        return dt

    @staticmethod
    def dt_generator(min_dt):
        base = math.ceil(math.fabs(math.log10(min_dt)))
        # new_dt = math.floor(min_dt*10**base) / 10**base
        return 10 ** (-base)

    def fabric_time_dependent(self):
        multiplier = np.ones(self.node_count)
        rho_water = self.materials.WATER['water: density [Kg/ m^3]']
        epsilon_water = self.materials.WATER['water: volume fraction bound [-]']

        rho_fiber = self.materials.FIBER['fiber: density dry [kg/m^3]']
        epsilon_fiber = self.materials.FIBER['fiber: volume fraction [-]']

        # time dependent fabric properties
        fabric_TD = \
            {'Temp [C]': self.materials.FABRIC['fabric: initial temp [C]'] * multiplier,
             'Temp [K]': self.materials.FABRIC['fabric: initial temp [K]'] * multiplier,
             'Specific Heat Effective [J/ Kg K]': self.materials.FABRIC[
                                                      'fabric: specific heat effective [J/ Kg K]'] * multiplier,
             'Thermal Conductivity Effective [W/ K m]': self.materials.FABRIC[
                                                            'fabric: thermal conductivity effective [W/ K m]'] * multiplier,
             'Density Effective [kg/m^3]': self.materials.FABRIC['fabric: density effective [kg/m^3]'] * multiplier,
             'VF: fiber [-]': self.materials.FIBER['fiber: volume fraction [-]'] * multiplier,
             'VF: water [-]': self.materials.WATER['water: volume fraction bound [-]'] * multiplier,
             'VF: gas [-]': self.materials.GAS['gas: volume fraction [-]'] * multiplier,
             'RH [-]': self.materials.FABRIC['fabric: initial RH [-]'] * multiplier,
             }

        fabric_TD['Diffusivity Effective [m^2/s]'] = self.materials.FABRIC[
                                                         'fabric: diffusivity effective [m^2 /s]'] * multiplier

        fabric_TD['Regain Factor: Instantaneous [-]'] = functions.regain_instant(epsilon_water,
                                                                                 epsilon_water,
                                                                                 rho_water,
                                                                                 rho_fiber) * multiplier

        fabric_DF = pd.DataFrame.from_dict(fabric_TD)
        names = [f'Fabric Node {x}' for x in range(self.node_count)]
        fabric_DF.index = names

        # fabric_DF = fabric_DF.reindex(sorted(fabric_DF.columns), axis=1)
        return fabric_DF

    def coupling_method(self, couple_choice):
        if self.Coupling_Option['Fully Decoupled'] == couple_choice:  # 1 -> DECOUPLED
            self.Solutions.water_vapor_concen_history[1, :] = self.Solutions.raw_water_concen_history[1, :]
            self.Solutions.rh_history[1, :] = \
                functions.relative_humidity_calc(self.Solutions.water_vapor_concen_history[1, :],
                                                 self.Solutions.temp_k_history[1, :])
            self.Solutions.rh_post_absorption_history[1, :] = self.Solutions.rh_history[1, :]

            for i in range(1, end_value):
                print_progressbar(i, end_value, prefix='Progress', suffix='Complete', length=50)
                # print(i)
                A, Phi = self.heat_matrix_ftcs_gen(t=i)
                self.Solutions.temp_k_history[i + 1, :] = self.FTCS_heat(A, Phi, t=i)

                B = self.diffusion_matrix_ftcs_gen()
                self.Solutions.raw_water_concen_history[i + 1, :] = self.FTCS_diffusion(B, t=i)
                self.Solutions.water_vapor_concen_history[i + 1, :] = self.Solutions.raw_water_concen_history[i + 1, :]
                self.Solutions.rh_history[i + 1, :] = \
                    functions.relative_humidity_calc(self.Solutions.water_vapor_concen_history[i + 1, :],
                                                     self.Solutions.temp_k_history[i + 1, :])
                self.Solutions.rh_post_absorption_history[i + 1, :] = self.Solutions.rh_history[i + 1, :]

        if self.Coupling_Option['Coupled Equilibrium'] == couple_choice:  # 2 -> COUPLED EQ

            # First time time through SCHEME
            A, Phi = PDE.heat_matrix_ftcs_gen()
            first_temp = np.append(PDE.Solutions.temp_k_history[0, :], PDE.BC_K['air: temp [K]'])
            first_temp[0] = PDE.BC_K['plate: temp [K]']

            # dont include last value, since its T_ambient
            PDE.Solutions.temp_k_history[0 + 1, :] = np.matmul(A, first_temp)[0:-1] + Phi

            B = PDE.diffusion_matrix_ftcs_gen()
            first_concen = np.append(PDE.Solutions.water_vapor_concen_history[0, :],
                                     PDE.BC_K['water: concentration in ambient air [kg/m^3]'])
            first_concen[0] = PDE.BC_K['water: concentration at plate [kg/m^3]']
            self.Solutions.rh_history[1, 0] = 1.00  # Fabric Surface Instantly at 100% RH

            PDE.solve_coupling(B, first_concen)
            self.post_concentration_sorption_flows()
            for i in range(1, self.count_to):
                print_progressbar(i, self.count_to, prefix='Progress', suffix='Complete', length=50)

                # update fabric properties due to regain

    def zeta_gen(self):
        # Heat transfer material coefficients based off dt, density, dx, and specific Heat

        dx_array = np.zeros(self.node_count)  # Create array of zeros for dx_array
        dx_array.fill(self.dx)  # fill dx array with interior node which is dx_array

        zeta_i = self.dt / (self.fabric_df['Density Effective [kg/m^3]'].values *
                            dx_array ** 2 * self.fabric_df['Specific Heat Effective [J/ Kg K]'].values)

        zeta_i[-1] = zeta_i[-1] * 2 * self.dx  # adjust last value for Mixed BC

        return zeta_i

    def radiation_transfer_coeff(self, t=0):
        # t is the current time step

        h = self.materials.FABRIC['fabric: radiation coefficient [-]'] * PI['sigma'] * \
            (self.Solutions.temp_k_history[t, -1] + self.BC_K['air: temp [K]']) * \
            (self.Solutions.temp_k_history[t, -1] ** 2 + self.BC_K['air: temp [K]'] ** 2)
        return h

    def heat_matrix_ftcs_gen(self, t=0):
        zeta = self.zeta_gen()  # run zeta_gen to produce zeta array
        h_r = self.radiation_transfer_coeff(t)
        h_c = PI['h_convection']

        n = self.node_count + 1  # plus 1 for ambient boundary temperature
        a = np.zeros([n, n])
        k = self.fabric_df['Thermal Conductivity Effective [W/ K m]']

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

        phi = q_net * zeta
        return a, phi

    def FTCS_heat(self, ft_matrix, phi, t=0):
        temp_old = np.append(self.Solutions.temp_k_history[t, :], self.BC_K['air: temp [K]'])  # attach BC to array
        temp_new = np.matmul(ft_matrix, temp_old)

        phi[0] = 0  # boundary condition cannot heat up!
        return temp_new[0:-1] + phi  # dont return ambient temp hence 0 to -1

    def gamma_gen(self):
        # Mass transfer material coefficients based off dt and evaporative diffusion resistance

        dx_array = np.zeros(self.node_count)  # Create array of zeros for dx_array
        dx_array.fill(self.dx)  # fill dx array with interior node which is dx_array
        gamma_i = self.dt * self.fabric_df['Diffusivity Effective [m^2/s]'].values / dx_array ** 2

        return gamma_i

    def diffusion_matrix_ftcs_gen(self):
        gamma = self.gamma_gen()  # run zeta_gen to produce zeta array
        still_air = self.materials.AIR['air: still air length over fabric [m]']
        n = self.node_count + 1  # plus 1 for ambient boundary temperature
        a = np.zeros([n, n])
        epsilon_gas = self.Solutions.volume_fraction_gas_history[self.counter, :]

        # TODO might want to preload these arrays so it doesnt need to do this every time step
        a[0, 0] = 1  # Concentration at Plate BC is invariant with time: steady state
        a[n - 1, n - 1] = 1  # Concentration at boundary air  is invariant with time: steady state

        for row in range(n - 1):
            if row != 0 and row != n - 2:  # dont modify 0th row. Nor modify last or 2nd to last row
                # A[row, row-1:row+2] = fo
                a[row, row - 1], a[row, row + 1] = gamma[row], gamma[row]
                # a[row, row] = fo[row]
                a[row, row] = epsilon_gas[row] - 2 * gamma[row]

            elif row != 0 and row == n - 2:  # second to last row modify for radiation and convection transfer
                a[row, row - 1], a[row, row + 1] = gamma[-1], gamma[-1] * self.dx / still_air
                a[row, row] = epsilon_gas[row] - gamma[-1] * (self.dx / still_air + 1)

        return a

    def m_dot(self, rho_vapor, node):
        # TODO Redundent code but still needed, FUNCTION?
        rho_water = self.materials.WATER['water: density [Kg/ m^3]']
        eps_water_prev = self.Solutions.volume_fraction_water_history[self.counter, node]
        eps_fiber = self.materials.FIBER['fiber: volume fraction [-]']
        rho_fiber = self.materials.FIBER['fiber: density dry [kg/m^3]']
        sorp_rate_factor = 16 * self.materials.FIBER['fiber: diffusion to diameter [1/s]']
        regain_inst = ((eps_water_prev * rho_water) / (eps_fiber * rho_fiber))

        temp_at_node_K = PDE.Solutions.temp_k_history[self.counter, node]  # Kelvin
        temp_at_node_C = temp_at_node_K - 273.15

        P_sat = functions.saturated_vapor_pressure(temp_at_node_C, None)[0]  # kPa
        P_sat *= 1000  # Pa

        # Pa
        Pa_v = rho_vapor * PI['R [J/mol K]'] * temp_at_node_K / self.materials.WATER['water: molecular weight [kg/mol]']
        RH = functions.RH_calc_basic(Pa_v, P_sat)
        regain_eq = self.materials.FIBER['fiber: regain[-]'] * functions.regain_function(RH)

        m_vs = sorp_rate_factor * eps_fiber * rho_fiber * (regain_eq - regain_inst)  # [kg / m^3]
        # m_dot_vs = 16 * self.materials.FIBER['fiber: diffusion to diameter [1/s]'] * \
        #            self.materials.FIBER['fiber: volume fraction [-]'] * \
        #            self.materials.FIBER['fiber: density dry [kg/m^3]']
        #
        # regain_eq = self.materials.FIBER['fiber: regain[-]'] * functions.absorption(self.Solutions.rh_history)
        #
        # epsilon_water = self.fabric_df['VF: water [-]']
        # epsilon_fiber = self.fabric_df['VF: fiber [-]']
        #
        # ones = np.ones(self.node_count)
        #
        # rho_water = self.materials.WATER['water: density [Kg/ m^3]'] * ones
        # rho_fiber = self.materials.FIBER['fiber: density [Kg/ m^3]'] * ones
        #
        # regain_instant = functions.regain_instant(epsilon_water, epsilon_fiber, rho_water, rho_fiber)
        # delta_regain = regain_eq - regain_instant
        # m_dot_vs *= delta_regain

        return m_vs, RH

    def coupling_equations(self, variables, *extra_args):
        eps_gas, eps_water, rho_vapor = variables  # solve these variables for next time step
        B_matrix, node, water_vapor_concen = extra_args

        rho_water = self.materials.WATER['water: density [Kg/ m^3]']
        eps_water_prev = self.Solutions.volume_fraction_water_history[self.counter, node]
        eps_gas_prev = self.Solutions.volume_fraction_gas_history[self.counter, node]
        eps_fiber = self.materials.FIBER['fiber: volume fraction [-]']
        rho_fiber = self.materials.FIBER['fiber: density dry [kg/m^3]']
        sorp_rate_factor = 16 * self.materials.FIBER['fiber: diffusion to diameter [1/s]']
        regain_inst = ((eps_water_prev * rho_water) / (eps_fiber * rho_fiber))

        temp_at_node_K = PDE.Solutions.temp_k_history[self.counter, node]  # Kelvin
        temp_at_node_C = temp_at_node_K - 273.15

        P_sat = functions.saturated_vapor_pressure(temp_at_node_C, None)[0]  # kPa
        P_sat *= 1000  # Pa

        # Pa
        Pa_v = rho_vapor * PI['R [J/mol K]'] * temp_at_node_K / self.materials.WATER['water: molecular weight [kg/mol]']
        RH = functions.RH_calc_basic(Pa_v, P_sat)
        regain_eq = self.materials.FIBER['fiber: regain[-]'] * functions.regain_function(RH)

        m_vs = sorp_rate_factor * eps_fiber * rho_fiber * (regain_eq - regain_inst)  # [kg / m^3]

        # continuity equation for volume fraction for future time step
        EQ1 = 1 - eps_gas - eps_water - eps_fiber
        EQ2 = rho_water * ((eps_water - eps_water_prev) / self.dt) - m_vs
        EQ3 = eps_gas * rho_vapor - np.dot(B_matrix[node, :], eps_gas_prev * water_vapor_concen) + EQ2 * self.dt

        # continuity equation for volume fraction for future time step
        # EQ1 = 1 - eps_gas - eps_water - eps_fiber
        #
        # term_1 = rho_water * (eps_water - eps_water_prev) / self.dt
        # term_2 = sorp_rate_factor * eps_fiber * rho_fiber * (regain_eq - regain_inst)
        # EQ2 = term_1 - term_2
        #
        # EQ3 = eps_gas * rho_vapor - np.dot(B_matrix[node, :], eps_gas_prev*water_vapor_concen) + EQ2

        return [EQ1, EQ2, EQ3]

    def solve_coupling(self, b_matrix, concen_array):
        # TODO FIGURE OUT WHY 2nd Node onward actually decrease their rho_vapor
        for node in range(1, self.node_count):
            x_0 = np.array([0.9, 0.007, 0.013])
            sol = fsolve(self.coupling_equations, x_0, args=(b_matrix, node, concen_array))
            # prior_values = np.array([PDE.fabric_df["VF: gas [-]"][1],
            #                          PDE.fabric_df["VF: water [-]"][1],
            #                          PDE.Solutions.water_vapor_concen_history[0, 1]])
            # print('Previous Values: ', prior_values)
            print('Solution Values', sol, '  Time: ', self.dt * self.counter)
            # self.coupling_equations(sol, b_matrix, node, concen_array)

            m_dot_vs, new_RH = self.m_dot(sol[2], node)
            PDE.Solutions.volume_fraction_gas_history[self.counter + 1, node] = sol[0]
            PDE.Solutions.volume_fraction_water_history[self.counter + 1, node] = sol[1]
            PDE.Solutions.water_vapor_concen_history[self.counter + 1, node] = sol[2]
            PDE.Solutions.rh_history[self.counter + 1, node] = new_RH
            PDE.Solutions.fabric_sorption_history[self.counter + 1, node] = m_dot_vs
            # RH = functions.relative_humidity_calc(sol[2] * 1000, PDE.Solutions.temp_k_history[self.counter, node])
            # print('\n' + f'Water Vapor Concen [kg/m^3] = {sol[2]}')
            # print(f'RH = {RH}')
            # rho_water = functions.concentration_calc([], RH, PDE.Solutions.temp_k_history[self.counter, node]) / 1000
            # print(f'My Function: {rho_water}')
            # print()
            # print('test')

    def FTCS_diffusion(self, ft_matrix, t=0):
        concentration_old = np.append(self.Solutions.water_vapor_concen_history[t, :],
                                      self.BC_K['water concentration in ambient air [g/m^3]'])  # attach BC to array

        # concentration_old[0] = PDE.BC_K['water concentration at plate [g/m^3]']  # always at max concen, slow heat
        concentration_new = np.matmul(ft_matrix, concentration_old)

        return concentration_new[0:-1]  # dont return ambient air concentration

    def dirichlet_step(self):
        #  Massive amount of parameters change once the fabric touches the plate

        delta_absorp = functions.absorption(self.BC_K['plate rh'] - self.fabric_IC['initial clothing rh'].array[0])[0]
        gain = self.fabric_IC['dry fabric density [g/m^3]'].array[0] * self.fabric_IC['regain'].array[0] * delta_absorp
        self.Solutions.fabric_sorption_history[1, 0] = gain  # increase of mass due to step change!

        # fiber_water_concen_history is in kg/m^3
        self.Solutions.fiber_water_concen_history[1:, 0] = self.Solutions.fiber_water_concen_history[0, 0] + gain / 1000

        self.Solutions.raw_water_concen_history[1:, 0] = self.BC_K['water concentration at plate [g/m^3]']
        self.Solutions.rh_history[1:, 0] = 1.0
        self.Solutions.rh_post_absorption_history[1:, 0] = 1.0  # dirichlet BC due to plate

        self.Solutions.water_vapor_concen_history[1:, 0] = self.BC_K['water concentration at plate [g/m^3]']
        self.Solutions.conden_concen_history[1:, 0] = 0

    def post_concentration(self, index):
        # REMEMBER INDEX ALREADY IS PLUS 1 for newest time step
        # temps at index are at current time step. Water heat flows 1 step behind

        # under normal conditions where RH of air > RH fabric
        # fiber sorption occurs, reduces ambient RH while absorption water, giving off heat
        # else, inverse occurs

        # if near 100% ambient RH, hard to tell if sorption occurs first then condensation or inverse

        current_concentration = self.Solutions.raw_water_concen_history[index, :]  # grams / m^3

        # RH before any sorption of water in fabric using temperature at current time step
        self.Solutions.rh_history[index, :] = \
            functions.relative_humidity_calc(current_concentration, self.Solutions.temp_k_history[index, :])

        # create passing tuple to feed into rh_equilibrium function
        # skips first parameter as the 1st param cannot change after step change
        passing_parameters = (
            self.fabric_IC, current_concentration[1:],
            self.Solutions.temp_k_history[index, 1:],
            self.Solutions.rh_post_absorption_history[index - 1, 1:]  # previous relative humidity
        )

        # Use Relative humidity equilibrium function
        relative_humidity_equilibrium, sorption, eq_air_concentration = functions.rh_equilibrium(*passing_parameters)

        # assumes sorption priority before condensation
        self.Solutions.rh_post_absorption_history[index, 1:] = relative_humidity_equilibrium

        # sorption
        self.Solutions.fabric_sorption_history[index, 1:] = sorption  # [g/m^3]

        # Updated water concentration inside the fabric [kg/m^3]
        # scale sorption from grams to kg / m^3
        self.Solutions.fiber_water_concen_history[index, 1:] = self.Solutions.fiber_water_concen_history[index - 1,
                                                               1:] + \
                                                               (sorption / 1000)

        # equilibrium air concentration due to equilibrium
        self.Solutions.water_vapor_concen_history[index, 1:] = eq_air_concentration

        # create passing tuple to feed into condensation_checker
        passing_parameters = (self.Solutions.rh_post_absorption_history[index, :],
                              self.Solutions.water_vapor_concen_history[index, :],
                              self.Solutions.temp_k_history[index, :],
                              self.Solutions.conden_concen_history[index, :])

        # Use condensation_checker to check for condensation using previous time steps temperatures
        self.Solutions.rh_post_absorption_history[index, :], \
        self.Solutions.water_vapor_concen_history[index, :], \
        self.Solutions.conden_concen_history[index, :] = functions.condensation_checker(*passing_parameters)

        return None

    def post_concentration_sorption_flows(self):
        # REMEMBER INDEX ALREADY IS PLUS 1
        index = self.counter + 1
        # TODO fix section for new method
        node_thickness = self.materials.FABRIC['fabric: thickness [m]']

        # heat of sorption due to regain in [J/g]
        self.Solutions.enthalpy_sorption_history[index, :] = functions.vectorized_h_sorp_cal(
            self.Solutions.rh_history[index - 1, :],
            self.Solutions.rh_history[index, :])

        flux = self.Solutions.fabric_sorption_history[index, :] * node_thickness / self.dt  # g /(m^2s)

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

        return None

    def post_plateflux_temp(self, ode_temp_solution, index):
        current_temperature = ode_temp_solution.y[1:, -1]
        self.temp_k_history[index, :] = current_temperature

        self.heat_flows[index, 0] = ode_temp_solution.y[0, -1]


if __name__ == '__main__':
    COUPLE_VAL = 2
    start = BC_IC.TIME_INPUT_PARAMETERS['start']
    finish_seconds = BC_IC.TIME_INPUT_PARAMETERS['finish [min]'] * 60
    material_pack = SolutionSetup.material_parameters()
    fabric_node_dimensions = functions.node_generator(True)

    dt_minimum = SolutionSetup.calculate_time_step(material_pack[0], min(fabric_node_dimensions))  # seconds
    time_step_lower_bound = SolutionSetup.dt_generator(dt_minimum)  # base 10 time step
    time_step = time_step_lower_bound * (dt_minimum // time_step_lower_bound)  # even rounds to lower decimal place

    tspan = np.arange(start, finish_seconds + time_step, time_step)
    t_size = tspan.shape[0]

    boundary_conditions = boundary_conditions_with_kelvin(BC_IC.BOUNDARY_INPUT_PARAMETERS)

    time_step = 0.005
    # initialize PDE CLASS arrays and key parameters
    PDE = SolutionSetup(material_pack, fabric_node_dimensions, t_size, time_step, boundary_conditions)

    PDE.coupling_method(COUPLE_VAL)

    # First time time through SCHEME
    # A, Phi = PDE.heat_matrix_ftcs_gen()
    # first_temp = np.append(PDE.Solutions.temp_k_history[0, :], PDE.BC_K['air: temp [K]'])
    # first_temp[0] = PDE.BC_K['plate: temp [K]']
    #
    # # dont include last value, since its T_ambient
    # PDE.Solutions.temp_k_history[0 + 1, :] = np.matmul(A, first_temp)[0:-1] + Phi
    #
    # B = PDE.diffusion_matrix_ftcs_gen()
    # first_concen = np.append(PDE.Solutions.water_vapor_concen_history[0, :],
    #                          PDE.BC_K['water: concentration in ambient air [kg/m^3]'])
    # first_concen[0] = PDE.BC_K['water: concentration at plate [kg/m^3]']
    #
    # PDE.solve_coupling(B, first_concen)

    # x_0 = np.array([0.9, 0.007, 0.013])
    # test_sol = fsolve(PDE.coupling_equations, x_0, args=(B, 1, first_concen))
    # print(test_sol - prior_values)

    # PDE.Solutions.raw_water_concen_history[0 + 1, :] = np.matmul(B, first_concen)[0: -1]  # dont include air boundary
    # PDE.dirichlet_step()
    # end_value = t_size - 1

    print('DONE')
    solution_to_df(PDE)
