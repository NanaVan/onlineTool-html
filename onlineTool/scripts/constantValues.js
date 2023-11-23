// parameter from CODATA 2018
const speed_c = 299792458; // [m/s] speed of light in vacuum
const elementary_charge = 1.602176634e-19; // [C] elementary charge
const me = 5.48579909065e-4; // [u] electron mass in u
const me_kg = 9.1093837015e-31; // [kg] electron mass in kg
const u2kg = 1.66053906660e-27; // [kg/u] amount of kg per 1 u
const eV2J = 1.602176634e-19; // [J/eV]
const MeV2u = 1.07354410233e-3; // [u/MeV] amount of u per 1 MeV

// elements
const elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds'];

// binding Energy of ions from NIST Atomic Spectra Database Ionization Energies Data
// format: 'element': {Z-charge: ionization energy [eV]}
var bindEnergy = {'H':{0: 13.598434599702}, 'He': {0: 24.587389011, 1: 54.417765486}, 'Li':{3: 5.391714996, 2: 75.6400970, 1: 122.45435913}, 'Be': {3: 18.21115, 2: 153.896205, 1: 217.71858459}, 'B': {3: 37.93059, 2: 259.3715, 1: 340.2260225}, 'C': {3: 64.49352, 2: 392.090518, 1: 489.99320779}, 'N': {3: 97.8901, 2: 552.06733, 1: 667.0461377}, 'O': {3: 138.1189, 2: 739.32683, 1: 871.409913}, 'F': {3: 185.1868, 2: 953.89805, 1: 1103.1175302}, 'Ne': {3: 239.0970, 2: 1195.80784, 1: 1362.199256}, 'Na': {3: 299.856, 2: 1465.134502, 1: 1648.702285}, 'Mg': {3: 367.489, 2: 1761.80488, 1: 1962.663889}, 'Al': {3: 442.005, 2: 2085.97702, 1: 2304.140359}, 'Si': {3: 523.415, 2: 2437.65815, 1: 2673.177958}, 'P': {3: 611.741, 2: 2816.90879, 1: 3069.842145}, 'S': {3: 706.994, 2: 3223.7807, 1: 3494.188518}, 'Cl': {3: 809.198, 2: 3658.3438, 1: 3946.29179}, 'Ar': {3: 918.375, 2: 4120.6657, 1: 4426.22407}, 'K': {3: 1034.542, 2: 4610.87018, 1: 4934.04979}, 'Ca': {3: 1157.726, 2: 5128.8578, 1: 5469.86358}, 'Sc': {3: 1287.957, 2: 5674.9037, 1: 6033.75643}, 'Ti': {3: 1425.257, 2: 6249.0226, 1: 6625.81023}, 'V': {3: 1569.656, 2: 6851.3109, 1: 7246.12624}, 'Cr': {3: 1721.183, 2: 7481.8628, 1: 7894.80289}, 'Mn': {3: 1879.873, 2: 8140.7872, 1: 8571.95438}, 'Fe': {3: 2045.759, 2: 8828.1879, 1: 9277.6886}, 'Co': {3: 2218.876, 2: 9544.1833, 1: 10012.1297}, 'Ni': {3: 2399.259, 2: 10288.8862, 1: 10775.3948}, 'Cu': {3: 2586.954, 2: 11062.4313, 1: 11567.6237}, 'Zn': {3: 2781.996, 2: 11864.9399, 1: 12388.9427}, 'Ga': {3: 2984.426, 2: 12696.5575, 1: 13239.5029}, 'Ge': {3: 3194.293, 2: 13557.4208, 1: 14119.4457}, 'As': {3: 3411.643, 2: 14447.678, 1: 15028.9251}, 'Se': {3: 3636.526, 2: 15367.491, 1: 15968.1075}, 'Br': {3: 3868.986, 2: 16317.011, 1: 16937.1497}, 'Kr': {3: 4109.083, 2: 17296.424, 1: 17936.2405}, 'Rb': {3: 4356.865, 2: 18305.884, 1: 18965.5484}, 'Sr': {3: 4612.397, 2: 19345.588, 1: 20025.2673}, 'Y': {3: 4875.731, 2: 20415.717, 1: 21115.588}, 'Zr': {3: 5146.935, 2: 21516.469, 1: 22236.712}, 'Nb': {3: 5426.066, 2: 22648.046, 1: 23388.850}, 'Mo': {3: 5713.194, 2: 23810.654, 1: 24572.213}, 'Tc': {3: 6008.391, 2: 25004.533, 1: 25787.047}, 'Ru': {3: 6311.721, 2: 26229.895, 1: 27033.564}, 'Rh': {3: 6623.262, 2: 27486.983, 1: 28312.031}, 'Pd': {3: 6943.097, 2: 28776.034, 1: 29622.678}, 'Ag': {3: 7271.298, 2: 30097.318, 1: 30965.780}, 'Cd': {3: 7607.95, 2: 31451.062, 1: 32341.587}, 'In': {3: 7953.14, 2: 32837.592, 1: 33750.404}, 'Sn': {3: 8306.95, 2: 34257.143, 1: 35192.501}, 'Sb': {3: 8669.48, 2: 35710.028, 1: 36668.183}, 'Te': {3: 9040.83, 2: 37196.522, 1: 38177.740}, 'I': {3: 9421.10, 2: 38716.996, 1: 39721.549}, 'Xe': {3: 9810.37, 2: 40271.724, 1: 41299.892}, 'Cs': {3: 10208.78, 2: 41861.075, 1: 42913.144}, 'Ba': {3: 10616.42, 2: 43485.366, 1: 44561.633}, 'La': {3: 11033.40, 2: 45144.996, 1: 46245.77}, 'Ce': {3: 11459.85, 2: 46840.306, 1: 47965.89}, 'Pr': {3: 11895.89, 2: 48571.71, 1: 49722.44}, 'Nd': {3: 12341.66, 2: 50339.59, 1: 51515.78}, 'Pm': {3: 12797.26, 2: 52144.29, 1: 53346.31}, 'Sm': {3: 13262.85, 2: 53986.12, 1: 55214.30}, 'Eu': {3: 13738.58, 2: 55865.92, 1: 57120.64}, 'Gd': {3: 14224.57, 2: 57783.90, 1: 59065.54}, 'Tb': {3: 14721.02, 2: 59739.3, 1: 61050.1}, 'Dy': {3: 15228.06, 2: 61736.56, 1: 63073.23}, 'Ho': {3: 15745.77, 2: 63772.43, 1: 65137.13}, 'Er': {3: 16274.56, 2: 65848.24, 1: 67241.48}, 'Tm': {3: 16814.34, 2: 67965.26, 1: 69387.45}, 'Yb': {3: 17365.44, 2: 70123.04, 1: 71574.63}, 'Lu': {3: 17928.05, 2: 72322.91, 1: 73804.35}, 'Hf': {3: 18502.32, 2: 74565.93, 1: 76077.70}, 'Ta': {3: 19088.51, 2: 76852.03, 1: 78394.63}, 'W': {3: 19686.74, 2: 79181.94, 1: 80755.91}, 'Re': {3: 20297.40, 2: 81556.90, 1: 83162.41}, 'Os': {3: 20920.60, 2: 83976.21, 1: 85614.42}, 'Ir': {3: 21556.60, 2: 86438.9, 1: 88113.6}, 'Pt': {3: 22205.7, 2: 88955.18, 1: 90659.84}, 'Au': {3: 22868.1, 2: 91515.82, 1: 93254.62}, 'Hg': {3: 23544.1, 2: 94124.70, 1: 95898.19}, 'Tl': {3: 24234.1, 2: 96783.21, 1: 98592.12}, 'Pb': {3: 24938.2, 2: 99491.85, 1: 101336.7}, 'Bi': {3: 25656.9, 2: 102251.76, 1: 104133.4}, 'Po': {3: 26390.4, 2: 105064.3, 1: 106983.4}, 'At': {3: 27139.0, 2: 107923.4, 1: 109887.2}, 'Rn': {3: 27903.1, 2: 110842.0, 1: 112842.2}, 'Fr': {3: 28683.4, 2: 113817.2, 1: 115857.5}, 'Ra': {3: 29479.8, 2: 116848.7, 1: 118929.5}, 'Ac': {3: 30293.1, 2: 119938.6, 1: 122063.1}, 'Th': {3: 31122.8, 2: 123086.4, 1: 125250.3}, 'Pa': {3: 31971.6, 2: 126296.6, 1: 128507}, 'U': {3: 32836.5, 2: 129570.3, 1: 131816.2}, 'Np': {3: 33722.2, 2: 132901.8, 1: 135202}, 'Pu': {3: 34625.8, 2: 136299.2, 1: 138640.2}, 'Am': {3: 35549.4, 2: 139769.5, 1: 142153.5}, 'Cm': {3: 36493.0, 2: 143299.6, 1: 145740.1}, 'Bk': {3: 37457.6, 2: 146904.7, 1: 149398}, 'Cf': {3: 38443.5, 2: 150579.3, 1: 153124}, 'Es': {3: 39451.4, 2: 154328.1, 1: 156927}, 'Fm': {3: 40482.2, 2: 158152.5, 1: 160808}, 'Md': {3: 41548, 1: 164764}, 'No': {3: 42632, 1: 168804}, 'Lr': {3: 43759, 1: 172928}, 'Rf': {1: 177142}, 'Db': {1: 181445}, 'Sg': {1: 185835}, 'Bh': {1: 190329}, 'Hs': {1: 194911}, 'Mt': {1: 199605}, 'Ds': {1: 204394}};

function toSci(xx, digi){
  //return xx.toExponential(digi).replace(/e\+?/, ' x 10^');
  return xx.toExponential(digi);
}

