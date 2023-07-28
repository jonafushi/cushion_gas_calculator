import CoolProp.CoolProp
import pandas as pd
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from scipy.stats import linregress
import math
import os
import sys


sys.stdout = open('console_output.txt', 'w')

# working directory
folder = os.getcwd()
# INPUT DATA ##########################################################################################################
# well parameters
re = 4920  # effective radius [ft] - equiv to 1500 metres
# rw = 0.25 # wellbore radius [ft] 6" tubing
S = 1  # skin factor - assumed as 1
delta_P = 1  # change in pressure between each row of data [psi]
Area = np.pi * re ** 2  # effective area [ft^2]
# Area_wb = np.pi * rw ** 2  # wellbore cross-sectional area [ft^2]
well_angle = 0  # degrees inclination of the well (0 = vertical)
# D_wb = rw * 2 * 12  # diameter of tubing in inches
epsilon = 0.0006  # pipe roughness assumed 0.0006 inches
dd_factor = 0.5  # wellhead pressure [psi] assumed to be 0.3 times the reservoir pressure - so drawdown factor is 0.3

# reservoir paraemters
rf = 0.9  # gas recovery factor
Ca = 31.62  # dietz shape factor

# standard conditions and surface temp
Tsc = 491.67  # temperature [Rankine]
Psc = 14.5038  # pressure [psia]
T_surface_C = 15  # surface temperature [celcius] - assumed to be 15 C (518.67 in rankine)
T_surface_R = T_surface_C * (
            9 / 5) + 491.67  # surface temperature [RANKINE] - assumed to be 15 C (518.67 in rankine)

# conversion factors
psi_to_pascal = 6894.76  # multiply for pascal / divide for psi
PaS_to_cP = 1000  # multiply for cP / divide for Pa S
kg_mol_to_lb_lbmol = 1000  # multiply for lb/lbmol / divide for kg/mol
kg_m3_to_lb_ft3 = 0.062428  # multiply for lb/ft3 / divide for kg/m3
secs_to_days = 86400  # multiply for days / divide for seconds
m3_to_ft3 = 35.3147  # multiply for feet / divide for metres
inches_to_metres = 0.0254
metres_to_feet = 3.28084  # multiply for feet / divide for metres
rankine_to_kelvin = 5 / 9  # multiply for kelvin / divide for rankine

# fluid parameters

HHV_H2 = 39.4  # kWh/kg
HHV_CH4 = 15.4  # kWh/kg
# function for specific gravity


def specific_gravity(gas):
    return CoolProp.CoolProp.PropsSI('D', 'T', 273.15, 'P', 100e5, gas) / CoolProp.CoolProp.PropsSI('D', 'T', 273.15, 'P', 100e5, 'AIR')

# Erosional velocity variables
# from API RP 14E (1984)
# c is a constant 100 for continuous service, 150 for intermittent service (assume intermittent for this well)
c = 150

### FUNCTIONS ##########################################################################################################

# calculate Z factor and make function for later use
def zfactorSI(temp, pressure, gas):
    return CoolProp.CoolProp.PropsSI('Z', 'T', temp, 'P', pressure, gas)

def pressure_col(Pr):
    Pwf = 1
    li = []
    while Pwf <= Pr + 1:
        li.append(Pwf)
        Pwf += delta_P
    return li

# calculate viscosity [cP] and make function for later use
def VISCcP(temp, pressure, gas):
    return CoolProp.CoolProp.PropsSI('V', 'T', temp, 'P', pressure, gas) * 1000

# calculate density and make function for later use [lb/ft3]
def dens_field(temp, pressure, gas):
    return CoolProp.CoolProp.PropsSI('D', 'T', temp, 'P', pressure, gas) * kg_m3_to_lb_ft3

# calculate density and make function for later use [kg/m3]
def dens_SI(temp, pressure, gas):
    return CoolProp.CoolProp.PropsSI('D', 'T', temp, 'P', pressure, gas)

# Pseduopressure - define function for 2(p/uz)
def two_puz(p, u, z):
    return 2 * (p / (u * z))

### FLOW RATE FUNCTIONS ###############################################################################################
# define flow rate equation- from gas reservoir engineering Lee & Wattenbarger
# Ppr - Ppwf = A*qg + B*qg^2

# define function for A from flow rate equation. Parameters are:
# temp [Kelvin], perm, thickness, Area, Dietz shape factor, wellbore radius, Skin
def A_flow(T, k, h, Area, Ca, rw, S):
    return (1.422 * (10 ** 6) * T / (k * h)) * (1.151 * np.log(10.06 * Area / (Ca * rw * rw)) - 0.75 + S)

# define function for B from flow rate equation. Parameters are:
# temp [Kelvin], perm, thickness, D (from flow rate equation)
def B(T, k, h, D):
    return 1.422 * (10 ** 6) * T * D / (k * h)

# define D from flow rate equation. Parameters are:
# beta, perm, molecular weight of gas/gas mixture, pressure @ STP, thickness, viscosity, wellbore radius, Temp@STP [Kelvin]
def D(beta, k, M, Psc, h, visc, rw, Tsc):
    return 2.715 * (10 ** -12) * beta * k * M * Psc / (h * visc * rw * Tsc)

# define beta from equation for D. Parameters are:
# perm, porosity
def beta(k, por):
    return 1.88 * (10 ** 10) * (k ** -1.47) * (por ** -0.53)

# define reaaranged flow rate equation to give qg. Parameters are:
# A, B, average reservoir pseudopressure, well flowing pseduopressure
def Forcheimer_qg(a, b, mr, mbh):
    return (-a + (((a ** 2) + 4 * b * (mr - mbh)) ** 0.5)) / (2 * b)

### TPR FUNCTIONS ######################################################################################################
# BHFP estimate : ptf is tubing head pressure, L is length of tubing, well_angle is angle of tubing (vertical = 0 degrees)
def BHFP_estimate(ptf, L, well_angle):
    return ptf + 0.25 * (ptf / 100) * ((L * np.cos(well_angle)) / 100)


# Calculate the Reynold's number: qg is flow rate [MSCF/d] , sg is specific gravity, d is diameter of the wellbore [inches], ug is viscosity in cP
def Nre(qg, sg, d, ug):
    nre = (20 * sg * qg) / (ug * d)
    if nre > 0:
        return nre
    else:
        return 100

# calculate the friction factor: Nre_tpr is the Reynold's number, epsilon is the pipe roughness [inches], D_wb wellbore diameter [inches]
def ff(Nre_tpr, epsilon, D_wb):
    if Nre_tpr > 2000:
        return 4 * (2.28 - 4 * np.log(
            (epsilon / D_wb) + (21.25 / Nre_tpr ** 0.9))) ** -2  # Jain and Swamee ref 13 CH 4
    else:
        return 64 / Nre_tpr

# calculate the variable s for the Pwf equation

def s(sg, L, well_angle, Zav_tpr, Tav_tpr):
    return 0.0375 * sg * L * np.cos(well_angle) / (Zav_tpr * Tav_tpr)

# pwf equation

def pwf_tpr(Pwh, s, qg, f_tpr, Tav_tpr, Zav_tpr, D_wb, well_angle):
    return (Pwh ** 2 * math.exp(s) + ((6.67e-4 * qg**2 * f_tpr * Tav_tpr**2 * Zav_tpr**2) / (D_wb**5 * np.cos(well_angle))) * (math.exp(s) - 1))**(1/2)

def TPR(qg, Pwh, L, Tres_rankine, sg, gas, D_wb):
        x = 1  # to keep an eye on the number of interations the while loop does later
        Z_over_Z = 1  # arbitrary number to kick off the loop
        V_over_V = 1  # arbitrary number to kick off the loop

        # step 1. estimate BHFP to start
        BHFP_est = BHFP_estimate(Pwh, L, well_angle)

        # step 2. use estimate to calculate average wellbore temp and pressure
        Tav_tpr = (Tres_rankine + T_surface_R) / 2
        print(str(Tav_tpr) + ' Rankine')
        Pav_tpr = (BHFP_est + Pwh) / 2
        print(str(Pav_tpr) + ' psi')

        # step 3. use av temp and pressure to calculate gas average viscosity and compressibility
        Zav_tpr = zfactorSI(Tav_tpr * rankine_to_kelvin, Pav_tpr * psi_to_pascal, gas)
        Viscav_tpr = VISCcP(Tav_tpr * rankine_to_kelvin, Pav_tpr * psi_to_pascal, gas)

        # step 4. calculate the Reynold's number
        Nre_tpr = Nre(qg, sg, D_wb, Viscav_tpr)
        print('Reynold number = ' + str(Nre_tpr))
        if Nre_tpr <= 2000:
            print('laminar flow')
        elif 2000 < Nre_tpr < 4000:
            print('unstable flow')
        else:
            print('turbulent flow')

        # step 5. calculate the friction factor
        f_tpr = ff(Nre_tpr, epsilon, D_wb)
        print('friction factor = ' + str(f_tpr))

        # set the tolerance of 0.5 % for our z and viscosity and begin the loop
        while ((round(Z_over_Z, 3) > 0.005) and (round(V_over_V, 3) > 0.005)):
            print('interation ' + str(x))  # keep an eye on the number of iterations it takes
            print('BHFP est ' + str(x) + ': ' + str(BHFP_est) + ' psi')

            # step 6. calculate s
            s_tpr = s(sg, L, well_angle, Zav_tpr, Tav_tpr)
            #s_tpr = s_tpr.astype(float)

            # step 7. calculate Pwf
            pwf_tpr_2 = pwf_tpr(Pwh, s_tpr, qg, f_tpr, Tav_tpr * rankine_to_kelvin, Zav_tpr, D_wb, well_angle)
            print('calculated BHFP est ' + str(x) + ': ' + str(pwf_tpr_2) + ' psi')

            # caculate z and visc new average pressure
            Pav_tpr_2 = (pwf_tpr_2 + Pwh) / 2

            Zav_tpr_2 = zfactorSI(Tav_tpr * rankine_to_kelvin, Pav_tpr_2 * psi_to_pascal, gas)
            Viscav_tpr_2 = VISCcP(Tav_tpr * rankine_to_kelvin, Pav_tpr_2 * psi_to_pascal, gas)

            # compare to original estimate
            Z_over_Z = (Zav_tpr_2 - Zav_tpr) / Zav_tpr_2
            V_over_V = (Viscav_tpr_2 - Viscav_tpr) / Viscav_tpr_2

            # reset variables
            Zav_tpr = Zav_tpr_2
            Viscav_tpr = Viscav_tpr_2
            x= x + 1


        return pwf_tpr_2


def TPR_Katz(sg, D, pwf, pwh, L, Tres_rankine, gas):

    Tav_tpr = (Tres_rankine + T_surface_R) / 2
    print(str(Tav_tpr) + ' Rankine')
    Pav_tpr = (pwf + pwh) / 2
    print(str(Pav_tpr) + ' psi')
    Zav_tpr = zfactorSI(Tav_tpr * rankine_to_kelvin, Pav_tpr * psi_to_pascal, gas)
    print(str(Zav_tpr) + ' Z factor')
    s_tpr = s(sg, L, well_angle, Zav_tpr, Tav_tpr)
    print(str(s_tpr) + ' s')
    f_tpr = (2*np.log(3.71/(epsilon/D)))**-2
    print('friction factor = ' + str(f_tpr))
    try:
        qg = (200000 *((s_tpr * D**5* (pwf**2 - math.exp(s_tpr) * pwh**2))/( sg * Tav_tpr * Zav_tpr * L * f_tpr * (math.exp(s_tpr)-1)))**0.5)/1e6
        return qg
    except:
        return 0

### EROSIONAL VELOCITY CALCULATION ####################################################################################

# define erosional velocity function from API RP 14E (1984)
# ve = c/(density)^0.5 where density will change with pwf and units are ft/s
# c is a constant 100 for continuous service, 150 for intermittent service
# ve * x-sectional area of wellbore = flow rate

# define function for ve -  [m/s]. Parameters are:
# constant c, well flowing pressure, reservoir temperature
def ve(c, density):
    return 1.22 * c / (density ** 0.5)

# calculate the velocity of the gas flow as per equation in Chapter 6, p. 62 from:
# Standard Handbook of Petroleum and Natural Gas Engineering doi: 10.1016/B978-0-12-383846-9.00006-0.
# ft/s
def gas_velocity_field(Q, Pb, Tf, Z, D, Pm, Tb):
    return (127.3 * 10 ** 3 * Q * Pb * Tf * Z / (D ** 2 * Pm * Tb)) / 60

## OPERATING POINT FUNCTIONS

def op_point_psi(qg_in, qg_out, pwf):
    try:
        first_line = LineString(np.column_stack((qg_in, pwf)))
        second_line = LineString(np.column_stack((qg_out, pwf)))
        intersection = first_line.intersection(second_line)
        op_point = intersection.x
        return op_point
    except:
        return 0

def op_point_pwf(qg_in, qg_out, pwf):
    try:
        first_line = LineString(np.column_stack((qg_in, pwf)))
        second_line = LineString(np.column_stack((qg_out, pwf)))
        intersection = first_line.intersection(second_line)
        op_point_pwf = intersection.y
        return op_point_pwf
    except:
        return 0

## EXPANSION FACTOR
def e_factor(z, T, P):
    return 35.3 * P / (z * T)

## PLATEAU RATE CALCULATION
def Qp(Qi, Tp, GIIPl):
    return 1 / ((1 / Qi) + (Tp / GIIPl))

## WELL POWER RATING EQUATION
def well_power_MW(gas, intersection_x):
    if gas == 'H2':
        return HHV_H2 * (intersection_x / m3_to_ft3 * 1000000) / 24 * PropsSI('D', 'T', Tsc * rankine_to_kelvin, 'P',
                                                                              Psc * psi_to_pascal, gas) / 1000
    elif gas == 'CH4':
        return HHV_CH4 * (intersection_x / m3_to_ft3 * 1000000) / 24 * PropsSI('D', 'T', Tsc * rankine_to_kelvin, 'P',
                                                                               Psc * psi_to_pascal, gas) / 1000
    else:
        return 0



########################################################################################################################
### MAIN CALCULATOR ####################################################################################################
def IPR_CALCULATOR(field_name, gas, Pr, k, poro, h, L, GIIP, Tres_rankine, Tp, D_wb, rw):

    # file setup field_name = 'XXX' gas = 'CH4'  # used to call fluid type from CoolProp - list:
    # http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
    title = gas + ' IPR calculation' + ' for the ' + str(field_name) + ' gas field'
    print(title)

    # some internal variables
    Tres_kelvin = Tres_rankine * 5 / 9  # reservoir temperature [Kelvin]
    Pwh = Pr * dd_factor  # wellhead pressure [psi] assumed to be 0.3 times the reservoir pressure
    sg = specific_gravity(gas)  # specific gravity of gas - air = 1
    print(str(gas) + ' SG = ' + str(round(sg, 5)))
    M = CoolProp.CoolProp.PropsSI('M', 'T', Tres_kelvin, 'P', (Psc * psi_to_pascal), gas) * kg_mol_to_lb_lbmol
    print(str(gas) + ' M = ' + str(round(M, 5)) + ' lb/lbmol')

## CREATING AND POPULATING THE DATAFRAME ##############################################################################

    # create a dataframe. will be added to later
    header = ['Z']
    df = pd.DataFrame(columns=header)

    # define pressure in previously defined increments (variable delta_P) up to reservoir pressure & add to dataframe
    df.insert(0, "pressure [psia]", pressure_col(Pr))
    print(df)

## COMPRESSIBILITY (Z), VISCOSITY, AND DENSITY CALCULATIONS ############################################################
    # add z factor to df
    df['Z'] = df.apply(lambda x: zfactorSI(Tres_kelvin, x['pressure [psia]'] * psi_to_pascal, gas), axis=1)

    # add viscosity to df
    df['VISC [cP]'] = df.apply(lambda x: VISCcP(Tres_kelvin, x['pressure [psia]'] * psi_to_pascal, gas), axis=1)

    # add density to df
    df['density [lb/ft3]'] = df.apply(lambda x: dens_field(Tres_kelvin, x['pressure [psia]'] * psi_to_pascal, gas), axis=1)

    # add density to df
    df['density [kg/m3]'] = df.apply(lambda x: dens_field(Tres_kelvin, x['pressure [psia]'] * psi_to_pascal, gas), axis=1)

    ## PSEUDOPRESSURE CALCULATION #######################################################################################
    # how to calculate pseudopressure: do the following steps as per Al-Husseiny trapezoid method
    # '2(p/uz)', 'delta P [psia]', '2(p/uz) * delta P', 'sum'

    # how to calculate pseudopressure: do the following steps as per Al-Husseiny trapezoid method
    # '2(p/uz)', 'delta P [psia]', '2(p/uz) * delta P', 'sum'

    df['2(p/uz) * delta P'] = df.apply(lambda x: two_puz(x['pressure [psia]'], x['VISC [cP]'], x['Z']),
                                       axis=1) * delta_P

    df['pseudopressure'] = df['2(p/uz) * delta P'].cumsum()
    # define the average reservoir pseudopressure
    Pres_pseudo = df['pseudopressure'].iloc[-1]

    ## FLOW RATE EQUATION - RESERVOIR INTO WELLBORE ######################################################################

    # make A a variable as it does not vary with pressure increments
    A = A_flow(Tres_rankine, k, h, Area, Ca, rw, S)

    # calculate the flow rate for each pressure increment and add it to dataframe
    df['qg [MMSCF/d]'] = df.apply(
        lambda x: Forcheimer_qg(A, B(Tres_rankine, k, h, D(beta(k, poro), k, M, Psc, h, x['VISC [cP]'], rw, Tsc)),
                                Pres_pseudo, x['pseudopressure']), axis=1)

    # define the AOFP/theoretical maximum rate qg @ pwf = 1 psi
    AOFP = df['qg [MMSCF/d]'].iloc[0]

    ## BHFP CALCULATION - FLOW RATE AT BOTTOM OF WELLBORE - KATZ METHOD

    df['TP_Qg [MMSCF/d]'] = df.apply(lambda x: TPR_Katz(sg, D_wb, x['pressure [psia]'], Pwh, L, Tres_rankine, gas), axis = 1)

    ## BHFP CALCULATION - PRESSURE AT BOTTOM OF WELLBORE


    df['TP_BHP [psia]'] = df.apply(lambda x: TPR(x['qg [MMSCF/d]'] * 1000, Pwh, L, Tres_rankine, sg, gas, D_wb), axis=1)

    ## EROSIONAL VELOCITY
    df['Vg [ft/s]'] = gas_velocity_field(df['qg [MMSCF/d]'], Pwh, Tres_rankine, df['Z'], D_wb, df['pressure [psia]'], T_surface_R)

    ## WELL DIAMETER
    df['Well diameter [inches]'] = D_wb

    # save df results ##################################################################################################
    print(df)
    csv_name =  field_name + '_' + str(Pr) + '_Pr_' + gas + '_' + str(D_wb) + '_inches_IPR.csv'
    print(csv_name)
    df.to_csv(csv_name)
    ####################################################################################################################

    ## CHECK INTERSECTION OF IPR AND TPR - needs to have try statement otherwise no intersect will throw errors

    first_line = LineString(np.column_stack((df['qg [MMSCF/d]'], df['pressure [psia]'])))
    #first_line = first_line.buffer(0)
    second_line = LineString(np.column_stack((df['qg [MMSCF/d]'], df['TP_BHP [psia]'])))
    #second_line = second_line.buffer(0)
    #try:
    intersection = first_line.intersection(second_line)


    ## OPERATING POINT COORDINATES
    operating_point = op_point_psi(df['qg [MMSCF/d]'], df['pressure [psia]'], df['TP_BHP [psia]'])
    operating_point_pwf = op_point_pwf(df['qg [MMSCF/d]'], df['pressure [psia]'], df['TP_BHP [psia]'])
    #except:
        #operating_point = None
        #operating_point_pwf = None


    ## GAS PLATEAU RATE CALCULATION ########################################################################################
    # lookup the number of wells in the field
    n_wells = wells_file.loc[(wells_file['field_name'] == field_name, 'n_wells')].iloc[0]


    # to calculate the GIIP for hydrogen we need to find the expansion factor ratio between methane and hydrogen at res conditions
    # expansion factor function z is compressibility factor, T is temperature in degrees rankine, P is pressure in psia
    CH4_e = e_factor(
        CoolProp.CoolProp.PropsSI('Z', 'T', Tres_kelvin, 'P', gas_fields_props['Pr_psia'].iloc[0] * psi_to_pascal, 'CH4'), Tres_rankine, gas_fields_props['Pr_psia'].iloc[0])
    H2_e = e_factor(
        CoolProp.CoolProp.PropsSI('Z', 'T', Tres_kelvin, 'P', gas_fields_props['Pr_psia'].iloc[0] * psi_to_pascal, 'H2'), Tres_rankine, gas_fields_props['Pr_psia'].iloc[0])
    e_ratio = H2_e / CH4_e
    # adjust giip for hydrogen
    if gas == 'H2':
        GIIP = GIIP * e_ratio
    # perform plateau rate calculation

    Q_plat = Qp((intersection.x * n_wells), Tp, GIIP)

## DRAWDOWN ############################################################################################################
    DD = Pr - intersection.y

## WORKING GAS AND CUSHION GAS VOLUMES #################################################################################
# assumes that there is 10% irrecoverable gas in the reservoir i.e. a recovery factor of 90%
    WGV = Q_plat * Tp
    CGV = (GIIP * rf) - WGV
    WG_CG_ratio = WGV / CGV
    CG_percentage = (CGV / (WGV + CGV)) * 100

## WELL POWER CALCULATION - how much energy does the well deliver
    well_power = well_power_MW(gas, intersection.x)

## EROSIONAL VELOCITY LIMITS

    Pwf_int = intersection.y
    dens_int_SI = dens_SI(Tres_kelvin, Pwf_int * psi_to_pascal, gas)
    ve_int_ft_s = ve(100, dens_int_SI) * metres_to_feet
    # convert to equivalent flow rate in MMSCF/d - do all in SI units then convert at end
    # to flow rate by multiplying velocity by area of wellbore and conversion to MMSCF/d
    V_limit_q = (ve_int_ft_s * (np.pi * D_wb ** 2) * secs_to_days) / 1000000
    # except:
    #     Q_plat = None
    #     DD = None
    #     WGV = 0
    #     CGV = (GIIP * rf) - WGV
    #     WG_CG_ratio = WGV / CGV
    #     CG_percentage = (CGV / (WGV + CGV)) * 100
    #     well_power = 0
    #     Pwf_int = None
    #     dens_int_SI = None
    #     ve_int_ft_s = 0
    #     # convert to equivalent flow rate in MMSCF/d - do all in SI units then convert at end
    #     # to flow rate by multiplying velocity by area of wellbore and conversion to MMSCF/d
    #     V_limit_q = None

    to_append = [field_name, n_wells, D_wb, gas, GIIP, Pr, DD, Tp, Q_plat, WGV, CGV, WG_CG_ratio, CG_percentage, ve_int_ft_s,
                 V_limit_q, well_power, operating_point, operating_point_pwf, AOFP]
    df_length = len(calc_df)
    calc_df.loc[df_length] = to_append
    calc_name = calc_df['field'].iloc[0] + '_' + calc_df['gas'].iloc[0] + '_' + str(
            calc_df['Well diameter [inches]'].iloc[0]) + '_inches'+ '_calculated_storage_props.csv'
    calc_df.to_csv(calc_name)

    # timesteps dataframe

    # add P/Z column for well forecasting
    p_over_z_init = gas_fields_props['Pr_psia'].iloc[0] / zfactorSI(gas_fields_props['Tres_R'].iloc[0] * rankine_to_kelvin, gas_fields_props['Pr_psia'].iloc[0] * psi_to_pascal, gas_fields_props['gas'].iloc[0])
    # p/z slope and intercept for Gp column
    x1 = 0
    x2 = GIIP  # MMSCF
    y1 = p_over_z_init
    y2 = 0
    slope, intercept, r_value, p_value, std_err = linregress([x1, x2], [y1, y2])

    # p over z each timestep
    p_over_z_timestep = df['pressure [psia]'].iloc[-1] / df['Z'].iloc[-1]

    # gas produced [MMSCF]
    Gp = (p_over_z_timestep - intercept) / slope

    # rate
    Qg = operating_point
    Pwf_op = operating_point_pwf
    to_append = [field_name, n_wells, D_wb, gas, GIIP, Pr, p_over_z_timestep, Gp, Qg, Pwf_op, e_ratio]
    df_length = len(df_steps)
    df_steps.loc[df_length] = to_append



    ## ## PLOT THE CURVES ####################################################################################################
    # plot the output (IPR & TPR curves) on the same graph
    # tells plot to share the axis i.e. add to one graph gca = get current axis
    ax = plt.gca()

    # IPR curve with label for gas type
    ipr_label = 'IPR - ' + str(gas)
    df.plot(kind='line', x='qg [MMSCF/d]', y='pressure [psia]', ax=ax, label=ipr_label, color='black')
    plt.xlim([0, df['qg [MMSCF/d]'].iloc[0] * 1.1])
    plt.ylim([0, Pr * 1.1])
    tpr_label = 'TPR ' + str(round(D_wb, 2)) + '" tubing ' + str(round(Pwh)) + ' psia'
    try:
        df.plot(x= 'qg [MMSCF/d]', y='TP_BHP [psia]', ax=ax, label=tpr_label, color='black', linestyle='--')
    except:
        pass

    ## INTERSECTION (IF IT INTERSECTS)
    #try:
    # define and show the intersection
    intersect_label = 'qg = ' + str(round(intersection.x, 4)) + ' MMSCF/d \nPwf = ' + str(round(intersection.y)) + ' psia'
    if intersection.geom_type == 'MultiPoint':
        plt.plot(*LineString(intersection).xy, 'o', label='intersect_label')
    elif intersection.geom_type == 'Point':
        plt.plot(*intersection.xy, 'o', label='operating rate')
        plt.annotate(intersect_label, (intersection.x, intersection.y * 1.1))
    #except:
      #  pass
    # sort out axis labels
    plt.xlabel('qg [MMSCF/d]')
    plt.ylabel('Pwf [psia]')

    # erosional velocity for intersection Pwf
    try:
        Ve_label = 'Ve limit: ' + str(round(ve_int_ft_s, 1)) + ' ft/s'
        plt.axvline(x=V_limit_q, label=Ve_label, linestyle='dashed', c='black', lw=0.75)
    except:
        pass



    # title and legend
    plt.legend(loc='lower left')
    plt.title(title)
    png_name = 'IPR_plot_' + field_name + '_' + str(Pr) + '_Pr_' + gas + str(D_wb) + '_inches.png'
    plt.savefig(png_name)
    plt.cla()
    plt.clf()



# import the wells file that contains the number of wells in each field
wells_file = pd.read_csv('INPUT/n_wells_single.csv')
print(wells_file)

# to automatically find all the input files
li = []
for root, dirs, files in os.walk(folder, topdown=False):
    for name in files:
        if "_timesteps" in name:
            print(os.path.join(root, name))
            x = os.path.join(root, name)
            li.append(x)
print(li)
list = pd.DataFrame(li)
list.to_csv('file_paths.csv')

for i in li:

## INPUT SETUP ##########################################################################################################

    # create a df to store key data calculated by model
    calcs_col_names = ["field", "n_wells", "Well diameter [inches]", "gas", "GIIP [MMSCF]", "initial reservoir pressure [psi]",
                       "drawdown [psi]", "plateau time [days]", "plateau rate [MMSCF/D]",
                       "WGV [MMSCF]", "CGV [MMSCF]", "WGV:CGV [ratio]", "Cushion gas requirement [%]",
                       "Ve limit [ft/s]", "Qg Ve limit [MMSCF/D]", "well power [MW]", "Qg Operating point [MMSCF/D]", "Pwf Operating point [psi]", "AOFP [MMSCF/d]"]
    calc_df = pd.DataFrame(columns=calcs_col_names)

    steps_col_names = ['field_name', 'n_wells', 'Well diameter [inches]', 'gas', 'GIIP [MMSCF]', 'Pr [psia]', 'p_over_z', 'Gp [MMSCF]',
                       'Qg [MMSCF/D]', 'Pwf [psia]', 'expansion ratio']

    df_steps = pd.DataFrame(columns=steps_col_names)

    print(df_steps)
    # x is used later to loop the naming function function
    x = 1

    # reads file to get input data
    gas_fields_props = pd.read_csv(i)

    # make a new folder to save shizzle in
    parent_dir = folder +'\OUTPUT'
    directory = str(gas_fields_props['field_name'].iloc[0]) + '_' + str(gas_fields_props['gas'].iloc[0]) + '_' + str(
            gas_fields_props['D_wb_inches'].iloc[0]) + '_inches'
    path = os.path.join(parent_dir, directory)
    os.mkdir(path)
    os.chdir(path)

    y = 1
    for i,r in gas_fields_props.iterrows():
        #plt.figure(y)
        IPR_CALCULATOR(r['field_name'], r['gas'], r['Pr_psia'],r['k_mD'],r['poro_frac'],r['h_ft'],r['L_ft'],r['GIIP_MMSCF'], r['Tres_R'], r['plateau'], r['D_wb_inches'], r['rw_ft'])
        print(str(y) + ' IPR curves constructed')
        # plt.show(block=False)
        # plt.pause(0.5)
        # plt.close(fig=int(y))
        y = y + 1

    # timestep calculations
    # average gas produced
    df_steps['delta Gp [MMSCF]'] = df_steps['Gp [MMSCF]'].diff()

    # average rate
    df_steps['Qg_av [MMSCF/D]'] = df_steps['Qg [MMSCF/D]'].rolling(window=2).mean()

    # number of days between each rate
    df_steps['delta_t [days]'] = df_steps['delta Gp [MMSCF]'] / df_steps['Qg_av [MMSCF/D]']

    # cumulative days - production time
    df_steps['t [days]'] = df_steps['delta_t [days]'].cumsum()
    timesteps_name = df_steps['field_name'].iloc[0] + '_timesteps_' + df_steps['gas'].iloc[0] + '_' + str(df_steps['Well diameter [inches]'].iloc[0]) + '_inches.csv'
    df_steps.to_csv(timesteps_name)

    df_steps.plot(kind='line', x='t [days]', y='Qg_av [MMSCF/D]', marker='s', legend=False, color='k', linewidth=1)
    plt.ylabel('Qg [MMSCF/D]')
    plt.xlim(0, )
    plt.ylim(0, )
    png_name = df_steps['field_name'].iloc[0] +'_forecast_' + df_steps['gas'].iloc[0] + '_' + str(df_steps['Well diameter [inches]'].iloc[0]) + '_inches.png'
    plt.savefig(png_name)
    plt.cla()
    plt.clf()

print('Finshed! ' + str(y) + ' IPR curves constructed')
sys.stdout.close()

