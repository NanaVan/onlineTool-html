#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import numpy as np

# three mode in CSRe
'''
vertical acceptance: A_nu [pi mm mrad] can be found in ring setting
average beta function: average_beta [m]
minimum of beta_y / beta_x: min_beta_xy can be found in ring setting
momentum acceptance: momentum_acceptance = delta P / P, can be found in ring setting
emittance: varepsilon [pi mm mrad], always < 10
'''
ring_modes = {
        'inter_target': {
            'A_nu': 75,
            'average_beta': 7.90744811,
            'min_beta_xy': 8.7/25.7,
            'momentum_acceptance': 0.02,
            'emittance': 6
            },
        'isochronous': {
            'A_nu': 20,
            'average_beta': 7.4971677,
            'min_beta_xy': 12.2/28.1,
            'momentum_acceptance': 0.007, 
            'emittance': 6
            },
        'normal': {
            'A_nu': 80,
            'average_beta': 7.9073749,
            'min_beta_xy': 8.2/17.6,
            'momentum_acceptance': 0.026, 
            'emittance': 6
            }
        }

# gas components in CSRe
'''
charge of gas: Z_t
partial voltage ratio: quantity 
'''
gas_CSRe_detail = {
        'H': {
            'Z_t': 1,
            'quantity': 0.827379
            },
        'C': {
            'Z_t': 6,
            'quantity': 0.043096
            },
        'N': {
            'Z_t': 7,
            'quantity': 0.081668
            },
        'O': {
            'Z_t': 8,
            'quantity': 0.047371
            },
        'He': {
            'Z_t': 2,
            'quantity': 3.75e-5
            },
        'Ar': {
            'Z_t': 18,
            'quantity': 9.25e-6
            }
        }
gas_CSRe = {
        'H': {
            'Z_t': 1,
            'quantity': 0.827379
            },
        'C': {
            'Z_t': 6,
            'quantity': 0.043096
            },
        'N': {
            'Z_t': 7,
            'quantity': 0.081668
            },
        'O': {
            'Z_t': 8,
            'quantity': 0.047371
            }
        }

# EC setting in CSRe
'''
vertical temperature of eletron beam: T_v, [eV] default 1
effective length of ECooling: len_cool, [m]
CSRe length: len_C, [m]
E-beam current: I_e, [A]
E-beam radius: r_0, [m]
E-beam energy: E_e, [eV]
'''
EC_CSRe = {
        'T_v': 1,
        'len_cool': 3.4,
        'len_C': 128.8,
        'I_e': 0.2,
        'r_0': 0.02,
        'E_e': 13.7145e3,
        }

# Binding Energy of ions from NIST Atomic Spectra Database Ionization Energies Data (limited to Li-like)
# format: 'element': {Z-charge: ionization energy [eV]}
bindEnergy = {
        'H':{0: 13.598434599702}, 
        'He': {0: 24.587389011, 1: 54.417765486}, 
        'Li':{3: 5.391714996, 2: 75.6400970, 1: 122.45435913}, 
        'Be': {3: 18.21115, 2: 153.896205, 1: 217.71858459}, 
        'B': {3: 37.93059, 2: 259.3715, 1: 340.2260225}, 
        'C': {3: 64.49352, 2: 392.090518, 1: 489.99320779}, 
        'N': {3: 97.8901, 2: 552.06733, '1': 667.0461377}, 
        'O': {3: 138.1189, 2: 739.32683, 1: 871.409913}, 
        'F': {3: 185.1868, 2: 953.89805, 1: 1103.1175302}, 
        'Ne': {3: 239.0970, 2: 1195.80784, 1: 1362.199256}, 
        'Na': {3: 299.856, 2: 1465.134502, 1: 1648.702285}, 
        'Mg': {3: 367.489, 2: 1761.80488, 1: 1962.663889}, 
        'Al': {3: 442.005, 2: 2085.97702, 1: 2304.140359}, 
        'Si': {3: 523.415, 2: 2437.65815, 1: 2673.177958}, 
        'P': {3: 611.741, 2: 2816.90879, 1: 3069.842145}, 
        'S': {3: 706.994, 2: 3223.7807, 1: 3494.188518}, 
        'Cl': {3: 809.198, 2: 3658.3438, 1: 3946.29179}, 
        'Ar': {3: 918.375, 2: 4120.6657, 1: 4426.22407}, 
        'K': {3: 1034.542, 2: 4610.87018, 1: 4934.04979}, 
        'Ca': {3: 1157.726, 2: 5128.8578, 1: 5469.86358}, 
        'Sc': {3: 1287.957, 2: 5674.9037, 1: 6033.75643}, 
        'Ti': {3: 1425.257, 2: 6249.0226, 1: 6625.81023}, 
        'V': {3: 1569.656, 2: 6851.3109, 1: 7246.12624}, 
        'Cr': {3: 1721.183, 2: 7481.8628, 1: 7894.80289}, 
        'Mn': {3: 1879.873, 2: 8140.7872, 1: 8571.95438}, 
        'Fe': {3: 2045.759, 2: 8828.1879, 1: 9277.6886}, 
        'Co': {3: 2218.876, 2: 9544.1833, 1: 10012.1297}, 
        'Ni': {3: 2399.259, 2: 10288.8862, 1: 10775.3948}, 
        'Cu': {3: 2586.954, 2: 11062.4313, 1: 11567.6237}, 
        'Zn': {3: 2781.996, 2: 11864.9399, 1: 12388.9427}, 
        'Ga': {3: 2984.426, 2: 12696.5575, 1: 13239.5029}, 
        'Ge': {3: 3194.293, 2: 13557.4208, 1: 14119.4457}, 
        'As': {3: 3411.643, 2: 14447.678, 1: 15028.9251}, 
        'Se': {3: 3636.526, 2: 15367.491, 1: 15968.1075}, 
        'Br': {3: 3868.986, 2: 16317.011, 1: 16937.1497}, 
        'Kr': {3: 4109.083, 2: 17296.424, 1: 17936.2405}, 
        'Rb': {3: 4356.865, 2: 18305.884, 1: 18965.5484}, 
        'Sr': {3: 4612.397, 2: 19345.588, 1: 20025.2673}, 
        'Y': {3: 4875.731, 2: 20415.717, 1: 21115.588}, 
        'Zr': {3: 5146.935, 2: 21516.469, 1: 22236.712}, 
        'Nb': {3: 5426.066, 2: 22648.046, 1: 23388.850}, 
        'Mo': {3: 5713.194, 2: 23810.654, 1: 24572.213}, 
        'Tc': {3: 6008.391, 2: 25004.533, 1: 25787.047}, 
        'Ru': {3: 6311.721, 2: 26229.895, 1: 27033.564}, 
        'Rh': {3: 6623.262, 2: 27486.983, 1: 28312.031}, 
        'Pd': {3: 6943.097, 2: 28776.034, 1: 29622.678}, 
        'Ag': {3: 7271.298, 2: 30097.318, 1: 30965.780}, 
        'Cd': {3: 7607.95, 2: 31451.062, 1: 32341.587}, 
        'In': {3: 7953.14, 2: 32837.592, 1: 33750.404}, 
        'Sn': {3: 8306.95, 2: 34257.143, 1: 35192.501}, 
        'Sb': {3: 8669.48, 2: 35710.028, 1: 36668.183}, 
        'Te': {3: 9040.83, 2: 37196.522, 1: 38177.740}, 
        'I': {3: 9421.10, 2: 38716.996, 1: 39721.549}, 
        'Xe': {3: 9810.37, 2: 40271.724, 1: 41299.892}, 
        'Cs': {3: 10208.78, 2: 41861.075, 1: 42913.144}, 
        'Ba': {3: 10616.42, 2: 43485.366, 1: 44561.633}, 
        'La': {3: 11033.40, 2: 45144.996, 1: 46245.77}, 
        'Ce': {3: 11459.85, 2: 46840.306, 1: 47965.89}, 
        'Pr': {3: 11895.89, 2: 48571.71, 1: 49722.44}, 
        'Nd': {3: 12341.66, 2: 50339.59, 1: 51515.78}, 
        'Pm': {3: 12797.26, 2: 52144.29, 1: 53346.31}, 
        'Sm': {3: 13262.85, 2: 53986.12, 1: 55214.30}, 
        'Eu': {3: 13738.58, 2: 55865.92, 1: 57120.64}, 
        'Gd': {3: 14224.57, 2: 57783.90, 1: 59065.54}, 
        'Tb': {3: 14721.02, 2: 59739.3, 1: 61050.1}, 
        'Dy': {3: 15228.06, 2: 61736.56, 1: 63073.23}, 
        'Ho': {3: 15745.77, 2: 63772.43, 1: 65137.13}, 
        'Er': {3: 16274.56, 2: 65848.24, 1: 67241.48}, 
        'Tm': {3: 16814.34, 2: 67965.26, 1: 69387.45}, 
        'Yb': {3: 17365.44, 2: 70123.04, 1: 71574.63}, 
        'Lu': {3: 17928.05, 2: 72322.91, 1: 73804.35}, 
        'Hf': {3: 18502.32, 2: 74565.93, 1: 76077.70}, 
        'Ta': {3: 19088.51, 2: 76852.03, 1: 78394.63}, 
        'W': {3: 19686.74, 2: 79181.94, 1: 80755.91}, 
        'Re': {3: 20297.40, 2: 81556.90, 1: 83162.41}, 
        'Os': {3: 20920.60, 2: 83976.21, 1: 85614.42}, 
        'Ir': {3: 21556.60, 2: 86438.9, 1: 88113.6}, 
        'Pt': {3: 22205.7, 2: 88955.18, 1: 90659.84}, 
        'Au': {3: 22868.1, 2: 91515.82, 1: 93254.62}, 
        'Hg': {3: 23544.1, 2: 94124.70, 1: 95898.19}, 
        'Tl': {3: 24234.1, 2: 96783.21, 1: 98592.12}, 
        'Pb': {3: 24938.2, 2: 99491.85, 1: 101336.7}, 
        'Bi': {3: 25656.9, 2: 102251.76, 1: 104133.4}, 
        'Po': {3: 26390.4, 2: 105064.3, 1: 106983.4}, 
        'At': {3: 27139.0, 2: 107923.4, 1: 109887.2}, 
        'Rn': {3: 27903.1, 2: 110842.0, 1: 112842.2}, 
        'Fr': {3: 28683.4, 2: 113817.2, 1: 115857.5}, 
        'Ra': {3: 29479.8, 2: 116848.7, 1: 118929.5}, 
        'Ac': {3: 30293.1, 2: 119938.6, 1: 122063.1}, 
        'Th': {3: 31122.8, 2: 123086.4, 1: 125250.3}, 
        'Pa': {3: 31971.6, 2: 126296.6, 1: 128507}, 
        'U': {3: 32836.5, 2: 129570.3, 1: 131816.2}, 
        'Np': {3: 33722.2, 2: 132901.8, 1: 135202}, 
        'Pu': {3: 34625.8, 2: 136299.2, 1: 138640.2}, 
        'Am': {3: 35549.4, 2: 139769.5, 1: 142153.5}, 
        'Cm': {3: 36493.0, 2: 143299.6, 1: 145740.1}, 
        'Bk': {3: 37457.6, 2: 146904.7, 1: 149398}, 
        'Cf': {3: 38443.5, 2: 150579.3, 1: 153124}, 
        'Es': {3: 39451.4, 2: 154328.1, 1: 156927}, 
        'Fm': {3: 40482.2, 2: 158152.5, 1: 160808}, 
        'Md': {3: 41548, 1: 164764}, 
        'No': {3: 42632, 1: 168804}, 
        'Lr': {3: 43759, 1: 172928}, 
        'Rf': {1: 177142}, 
        'Db': {1: 181445}, 
        'Sg': {1: 185835}, 
        'Bh': {1: 190329}, 
        'Hs': {1: 194911}, 
        'Mt': {1: 199605}, 
        'Ds': {1: 204394}
        }

class life_estimation(object):
    '''
    estimate 1/e life time of the ion

    ion_mass:       the mass of ion to be calculate, [amu]
    ion_charge:     the charge of the ion, e.g. 26 For Fe26+
    atomic_number:  the atomic number of the ion, e.g. 26 For Fe 
    ion_gamma:      the Lorentzian factor of the ion, e.g. 1.1049
    gas_pressure:   the gas pressure in ring, [mbar], e.g. 1e-11
    ring_mode:  selected from ring_modes, e.g. ring_modes['inter_target']
    EC_setting: setting for ECooling, e.g. EC_CSRe
    EC_on:      True for ECooling on, False for off
    gas_in_ring:gas components of the ring, e.g. gas_CSRe
    normal_temperature: temperature to calculate, default 20 centigrade
    '''

    r_e = 2.8179403262e-15 # [m] classical electron radius
    k_b = 1.380649e-23 # [J/K] Boltzmann's constant
    N_A = 6.02214076e23 # [mol^(-1)] Avogadro constant
    c = 299792458 # [m/s] speed of light in vacuum
    e = 1.602176634e-19 # [C] elementary charge
    me = 9.1093837015e-31 # [kg] electron mass
    eV2J = 1.602176634e-19 # [J/eV]
    u = 1.66053906660e-27 # [kg] atomic mass unit 
    MeV2u = 1.07354410233e-3 # [u/MeV] amount of u per 1 MeV
    elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds']

    def __init__(self, ion_mass, ion_charge, atomic_number, ion_gamma, gas_pressure, ring_mode, EC_setting, EC_on=False, gas_temperature=25, gas_in_ring=gas_CSRe, bindEnergy=bindEnergy):
        self.m = ion_mass
        self.ion_element = self.elements[atomic_number-1]
        self.r_i = self.me * self.r_e / (self.m * self.u)
        self.Q = ion_charge
        self.Z = atomic_number
        self.gamma = ion_gamma
        self.beta = np.sqrt(1-1/self.gamma**2)
        self.gas_pressure = gas_pressure # [mbar]
        self.T_normal = 273.15 + gas_temperature # [K]
        self.gas_in_ring = gas_in_ring
        self.EC_on = EC_on
        self.ring_mode = ring_mode
        self.bindEnergy = bindEnergy
        if self.EC_on:
            self.EC_setting = EC_setting

    def result_show(self, result_print=True, mode='automatic', EC_method='Dmitriev', EL_method='Habs'):
        '''
        show 1/e life time result
        '''
        lambda_total = 0
        if self.EC_on:
            lambda_rec = self.calc_life_recombination(self.EC_setting['T_v'], self.EC_setting['len_cool'], self.EC_setting['len_C'], self.EC_setting['I_e'], self.EC_setting['r_0'])
            if result_print: print('tau_rec: {:.4E} [s] / {:.4E} [min] / {:.4E} [h]'.format(1/lambda_rec, 1/lambda_rec/60, 1/lambda_rec/3600))
            lambda_total += lambda_rec
        lambda_ss, lambda_ms = self.calc_life_Scattering(self.ring_mode['A_nu'], self.ring_mode['average_beta'], self.ring_mode['emittance'], self.gas_in_ring)
        if result_print:
            print('tau_ss: {:.4E} [s] / {:.4E} [min] / {:.4E} [h]'.format(1/lambda_ss, 1/lambda_ss/60, 1/lambda_ss/3600))
            print('tau_ms: {:.4E} [s] / {:.4E} [min] / {:.4E} [h]'.format(1/lambda_ms, 1/lambda_ms/60, 1/lambda_ms/3600))
        lambda_total += (lambda_ss + lambda_ms)
        if mode == 'manual':
            lambda_ec = self.calc_life_eletronCapture_manual(self.gas_in_ring, method=EC_method)
            if result_print: print('tau_ec: {:.4E} [s] / {:.4E} [min] / {:.4E} [h]'.format(1/lambda_ec, 1/lambda_ec/60, 1/lambda_ec/3600))
            lambda_total += lambda_ec
            if self.Z > self.Q:
                lambda_el = self.calc_life_eletronLoss_manual(self.gas_in_ring, method=EL_method)
                if result_print: print('tau_el: {:.4E} [s] / {:.4E} [min] / {:.4E} [h]'.format(1/lambda_el, 1/lambda_el/60, 1/lambda_el/3600))
                lambda_total += lambda_el
        else:
            lambda_ec = self.calc_life_eletronCapture_automatic(self.gas_in_ring)
            if result_print: print('tau_ec: {:.4E} [s] / {:.4E} [min] / {:.4E} [h]'.format(1/lambda_ec, 1/lambda_ec/60, 1/lambda_ec/3600))
            lambda_total += lambda_ec
            if self.Z > self.Q:
                lambda_el = self.calc_life_eletronLoss_automatic(self.gas_in_ring)
                if result_print: print('tau_el: {:.4E} [s] / {:.4E} [min] / {:.4E} [h]'.format(1/lambda_el, 1/lambda_el/60, 1/lambda_el/3600))
                lambda_total += lambda_el
        lambda_nr = self.calc_life_nuclearReaction(self.gas_in_ring)
        if result_print: print('tau_nr: {:.4E} [s] / {:.4E} [min] / {:.4E} [h]'.format(1/lambda_nr, 1/lambda_nr/60, 1/lambda_nr/3600))
        lambda_total += lambda_nr
        if result_print:
            print('1/e lifetime:\ntau: {:.2E} [s] / {:.2E} [min] / {:.2E} [h]'.format(1/lambda_total, 1/lambda_total/60, 1/lambda_total/3600))
            print("ES: {:.2f} %, EC: {:.2f} %, NR: {:.2f} %".format((lambda_ss+lambda_ms)/lambda_total*100, lambda_ec/lambda_total*100, lambda_nr/lambda_total*100))
        if self.EC_on and result_print:
            print("REC: {:.2f} %".format(lambda_rec/lambda_total*100))
        if self.Z > self.Q and result_print:
            print("EL: {:.2f} %".format(lambda_el/lambda_total*100))
        return 1/lambda_total # [s]
    
    def calc_life_recombination(self, T_v, len_cool, len_C, I_e, r_0):
        '''
        estimate 1/e life time of the ion from recombination
        recombination@ H.Poth. Electron cooling: Theory, experiment, application. Physics Reports. 196 (1990): 135-297

        T_v:        vertical temperature of eletron beam, [eV]
                    usually set to be 1eV, since the original temperature from E-gun is around 0.1~0.2 eV
        len_cool:   effective length of ECooling, [m]
        len_C:      CSRe length, [m]
        I_e:        E-beam current, [A]
        r_0:        E-beam radius, [m]
        '''
        alpha_rec = 3.02e-13 * self.Q**2 / T_v * (np.log(11.32 * self.Q / np.sqrt(T_v)) + 0.14 * (T_v / self.Q**2)**(1/3)) # [cm^3/s]
        #e_gamma = 1+ E_e * self.eV2J / self.me / self.c**2
        e_gamma = self.gamma
        rho = I_e / np.pi / r_0**2 / self.e / self.c / np.sqrt(1 - 1 / e_gamma**2) * 1e-6 # [cm^(-3)]
        return len_cool * alpha_rec * rho / self.gamma**2 / len_C 

    def calc_life_Scattering(self, A_nu, average_beta, emittance, gas_in_ring):
        '''
        estimate 1/e life time of the ion from Scattering
        with assumption: round vacuum chamber, beta_x = beta_z = <beta_y>, varepsilon_x = varepsilon_z = varepsilon
        single-scattering@ H.Poth. Electron cooling: Theory, experiment, application. Physics Reports. 196 (1990): 135-297
        multiple-scattering@ B.Franzke. Interaction of stored ion beams with the residual gas. CAS-CERN Accelerator School: 4th advanced accelerator physics course. (1992): 100-119

        A_nu:                   vertical acceptance, can be found in ring setting, [pi mm mrad]
        average_beta:           average beta function, [m]
        emittance:              initial emittance, [pi mm mrad]
        '''
        lambda_ss, lambda_ms = 0, 0
        sqr_angle_acceptance = A_nu * 1e-6 / average_beta * (1 - emittance / A_nu)**2
        for gas_component in gas_in_ring.keys():
            Z_t = gas_in_ring[gas_component]['Z_t']
            quantity = gas_in_ring[gas_component]['quantity']
            rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T 
            tau_ss_temp = self.beta**3 * self.gamma**2 * sqr_angle_acceptance / (4 * np.pi * Z_t**2 * self.Q**2 * self.r_i**2 * rho * self.c)
            tau_ms_temp = self.beta**3 * self.gamma**2 * (A_nu - emittance) * 1e-6 / (32 * average_beta * Z_t**2 * self.Q**2 * self.r_i**2 * rho * self.c) / np.log(204 / Z_t**(1/3))
            lambda_ss += 1/tau_ss_temp
            lambda_ms += 1/tau_ms_temp
        return lambda_ss, lambda_ms 

    def calc_life_eletronCapture_automatic(self, gas_in_ring):
        '''
        estimate 1/e life time of the ion from Electron Capture (EC)
        EC is one of the main mechanism of the beam loss in the ring.
	   	Ion's magnetic rigidity will be changed when capturing the electron from the residual gas.

	    default: this article gives a good result when Z_t < 18, ion_energy from 40-1000 MeV/u
	    single-EC@ I.S. Dmitriev, et al. On the target thickness to attain equilibrium charge distribution in a beam of fast ions. Nucl. Instr. Meth. Phys. Res. B. 14. (1986): 515-526
	   	option: this article gives a good result when ion_energy[keV/amu]/Z_t^1.25/ion_charge^0.7 from 10-1000, ion_charge >= 3
	   	EC@ A.S. Schlacher, et al. Electron capture for fast highly charged ions in gas targets: an empirical scaling rule. Phys. Rev. A. 27. (1983) 11: 3372
	   	option: this article gives a good result when ion_Z >= 36; gamma = (1 + T/931.5) <= 1.1
		single-EC@ B. Franzke. Vacuum requirements for heavy ion synchrotrons. IEEE Trans. Nucl. Sci. NS-28 (1981): 2116-2118        
        '''
        lambda_ec = 0
        if ((self.gamma - 1) / self.MeV2u * 1e3 / self.Q**0.7 / 18**1.25 > 10 and (self.gamma - 1) / self.MeV2u * 1e3 /self.Q**0.7 < 1000 and self.Q >= 3):
            print("EC: Schlacher Method")
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                E_temp = (self.gamma - 1) / self.MeV2u * 1e3 / self.Z**0.7 / Z_t**1.25 # [keV]
                lambda_ec += 1.1e-8 / E_temp**4.8 * (1 - np.exp(-0.037 * E_temp**2.2)) * (1 - np.exp(-2.44e-5 * E_temp**2.6)) * 1e-4 * rho * self.beta * self.c * self.Q**0.5 / Z_t**1.8 # [s^(-1)]
        elif (self.Z >= 36 and self.gamma <= 1.1):
            print("EC: Franzke Method")
            q_bar = self.Z * (1 - np.exp(-137 * self.beta / self.Z**0.67))
            if self.Q < q_bar:
                b = 4
            else:
                b = 2
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                lambda_ec += 2e-28 * self.Z**0.5 * q_bar**2 * Z_t * (1 - np.exp(-137 * self.beta / Z_t**0.67)) * (self.gamma**2 - 1)**(-2) * (self.Q / q_bar)**b * rho * self.beta * self.c # [s^(-1)]
        else:
            default_state = True
            default_state *= ((self.gamma - 1) / self.MeV2u) >= 40 and ((self.gamma - 1) / self.MeV2u) <= 1000
            if self.Q <= 18:
                K = 1
            else:
                K = 1.20 - 0.01 * self.Z
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                lambda_ec += K * 1e-33 * Z_t * self.Z**5 * self.beta**(-9) * 1e-4 * rho * self.c
                default_state *= (Z_t < 18)
            if default_state:
                print("EC: Dmitriev Method")
            else:
                print("EC: (Dmitriev Method)")
        return lambda_ec

    def calc_life_eletronCapture_manual(self, gas_in_ring, method='Dmitriev'):
        '''
        estimate 1/e life time of the ion from Electron Capture (EC)
        EC is one of the main mechanism of the beam loss in the ring.
	   	Ion's magnetic rigidity will be changed when capturing the electron from the residual gas.

	    default: this article gives a good result when Z_t < 18, ion_energy from 40-1000 MeV/u
	    single-EC@ I.S. Dmitriev, et al. On the target thickness to attain equilibrium charge distribution in a beam of fast ions. Nucl. Instr. Meth. Phys. Res. B. 14. (1986): 515-526
	   	option: this article gives a good result when ion_energy[keV/amu]/Z_t^1.25/ion_charge^0.7 from 10-1000, ion_charge >= 3
	   	EC@ A.S. Schlacher, et al. Electron capture for fast highly charged ions in gas targets: an empirical scaling rule. Phys. Rev. A. 27. (1983) 11: 3372
	   	option: this article gives a good result when ion_Z >= 36; gamma = (1 + T/931.5) <= 1.1
		single-EC@ B. Franzke. Vacuum requirements for heavy ion synchrotrons. IEEE Trans. Nucl. Sci. NS-28 (1981): 2116-2118   
        '''
        lambda_ec = 0
        if method == 'Schlacher':
            print("(manual) EC: Schlacher Method")
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                E_temp = (self.gamma - 1) / self.MeV2u * 1e3 / self.Z**0.7 / Z_t**1.25 # [keV]
                lambda_ec += 1.1e-8 / E_temp**4.8 * (1 - np.exp(-0.037 * E_temp**2.2)) * (1 - np.exp(-2.44e-5 * E_temp**2.6)) * 1e-4 * rho * self.beta * self.c * self.Q**0.5 / Z_t**1.8 # [s^(-1)]
        elif method == 'Franzke':
            print("(manual) EC: Franzke Method")
            q_bar = self.Z * (1 - np.exp(-137 * self.beta / self.Z**0.67))
            if self.Q < q_bar:
                b = 4
            else:
                b = 2
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                lambda_ec += 2e-28 * self.Z**0.5 * q_bar**2 * Z_t * (1 - np.exp(-137 * self.beta / Z_t**0.67)) * (self.gamma**2 - 1)**(-2) * (self.Q / q_bar)**b * rho * self.beta * self.c # [s^(-1)]
        else:
            default_state = True
            default_state *= ((self.gamma - 1) / self.MeV2u) >= 40 and ((self.gamma - 1) / self.MeV2u) <= 1000
            if self.Q <= 18:
                K = 1
            else:
                K = 1.20 - 0.01 * self.Z
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                lambda_ec += K * 1e-33 * Z_t * self.Z**5 * self.beta**(-9) * 1e-4 * rho * self.c
                default_state *= (Z_t < 18)
            if default_state:
                print("(manual) EC: Dmitriev Method")
            else:
                print("(manual) EC: (Dmitriev Method)")
        return lambda_ec

    def calc_life_eletronLoss_automatic(self, gas_in_ring):
        '''
        estimate 1/e lifetime of the ion from Electron Stripping (EL)
		EL is one of the main mechanism of the beam loss in the ring.
		EL is forbidden for the bare ion.
		with assumption: ion_beta / alpha > ion_Z

		default: ion velocity beta > Z/137
		EL@ D. Habs, et al. First experiments with the Heidelberg test storage ring TSR. Nucl. Instr. Meth. Phys. Res. B. 43. (1989): 390-410
		option: this article gives a good result when ion_Z >= 36; gamma = (1 + T/931.5) <= 1.1
		single-EL@ B. Franzke. Vacuum requirements for heavy ion synchrotrons. IEEE Trans. Nucl. Sci. NS-28 (1981): 2116-2118
        '''
        alpha_0 = 5.29177210903e-11 # [m], Bohr radius
        alpha = 7.2973525693e-3 # fine-structure constant

        ion_num = self.Z - self.Q
        lambda_el = 0
        if (self.Z >=36 and self.gamma <= 1.1):
            print("EL: Franzke Method")
            q_bar = self.Z * (1 - np.exp(-137 * self.beta * self.Z**0.67))
            if self.Q < q_bar:
                b = -2.3
            else:
                b = -4
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                lambda_el += 3.5 * 10**(-18 + (0.71 * np.log(self.Q))**1.5) * 1e-4 * q_bar**(-2) * Z_t * (1 - np.exp(-137 * self.beta / Z_t**0.67)) * (self.gamma**2 - 1)**(-0.5) * (self.Q / q_bar)**b
        else:
            if self.beta > self.Z / 137:
                print("EL: Habs Method")
            else:
                print("EL: (Habs Method)")
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                if self.Z >= Z_t:
                    I_ratio = 0
                    for j in range(1,ion_num+1):
                        I_ratio += self.bindEnergy['H'][0] / self.bindEnergy[self.ion_element][j]
                    lambda_el += 4 * np.pi * alpha_0**2 * (alpha / self.beta)**2 * (Z_t + 1) * Z_t * I_ratio * rho * self.beta * self.c
                elif self.Z <= Z_t**(1/3):
                    I_ratio = 0
                    for j in range(1,ion_num+1):
                        I_ratio += np.sqrt(self.bindEnergy['H'][0] / self.bindEnergy[self.ion_element][j])
                    lambda_el += np.pi * alpha_0**2 * alpha / self.beta * Z_t**(2/3) * I_ratio * rho * self.beta * self.c
                else:
                    pass
        return lambda_el

    def calc_life_eletronLoss_manual(self, gas_in_ring, method='Habs'):
        '''
        estimate 1/e lifetime of the ion from Electron Stripping (EL)
		EL is one of the main mechanism of the beam loss in the ring.
		EL is forbidden for the bare ion.
		with assumption: ion_beta / alpha > ion_Z

		default: ion velocity beta > Z/137
		EL@ D. Habs, et al. First experiments with the Heidelberg test storage ring TSR. Nucl. Instr. Meth. Phys. Res. B. 43. (1989): 390-410
		option: this article gives a good result when ion_Z >= 36; gamma = (1 + T/931.5) <= 1.1
		single-EL@ B. Franzke. Vacuum requirements for heavy ion synchrotrons. IEEE Trans. Nucl. Sci. NS-28 (1981): 2116-2118
        '''
        alpha_0 = 5.29177210903e-11 # [m], Bohr radius
        alpha = 7.2973525693e-3 # fine-structure constant

        ion_num = self.Z - self.Q
        lambda_el = 0
        if method=='Habs':
            print("(manual) EL: Franzke Method")
            q_bar = self.Z * (1 - np.exp(-137 * self.beta * self.Z**0.67))
            if self.Q < q_bar:
                b = -2.3
            else:
                b = -4
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                lambda_el += 3.5 * 10**(-18 + (0.71 * np.log(self.Q))**1.5) * 1e-4 * q_bar**(-2) * Z_t * (1 - np.exp(-137 * self.beta / Z_t**0.67)) * (self.gamma**2 - 1)**(-0.5) * (self.Q / q_bar)**b
        else:
            if self.beta > self.Z / 137:
                print("(manual) EL: Habs Method")
            else:
                print("(manual) EL: (Habs Method)")
            for gas_component in gas_in_ring.keys():
                Z_t = gas_in_ring[gas_component]['Z_t']
                quantity = gas_in_ring[gas_component]['quantity']
                rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T
                if self.Z >= Z_t:
                    I_ratio = 0
                    for j in range(1,ion_num+1):
                        I_ratio += self.bindEnergy['H'][0] / self.bindEnergy[self.ion_element][j]
                    lambda_el += 4 * np.pi * alpha_0**2 * (alpha / self.beta)**2 * (Z_t + 1) * Z_t * I_ratio * rho * self.beta * self.c
                elif self.Z <= Z_t**(1/3):
                    I_ratio = 0
                    for j in range(1,ion_num+1):
                        I_ratio += np.sqrt(self.bindEnergy['H'][0] / self.bindEnergy[self.ion_element][j])
                    lambda_el += np.pi * alpha_0**2 * alpha / self.beta * Z_t**(2/3) * I_ratio * rho * self.beta * self.c
                else:
                    pass
        return lambda_el

    def calc_life_nuclearReaction(self, gas_in_ring):
        '''
        estimate 1/e lifetime of the ion from Nuclear Reaction (NR)
		
		NR @ X.Y. Zhang. Study of Lifetime of Light Ion Beam stored in CSRm for Internal Target Experiment. (2005) MSc thesis
        '''
        ion_A = np.round(self.m)
        lambda_nr = 0
        for gas_component in gas_in_ring.keys():
            Z_t = gas_in_ring[gas_component]['Z_t']
            quantity = gas_in_ring[gas_component]['quantity']
            rho = self.gas_pressure * 1e2 * quantity / self.k_b / self.T_normal # from General Gas Equation P = rho*k_b*T 
            lambda_nr += np.pi * (ion_A**(1/3) + (2*Z_t)**(1/3))**2 * 1e-30 / self.beta * rho * self.c
        return lambda_nr
