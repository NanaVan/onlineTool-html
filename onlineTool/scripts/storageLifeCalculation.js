function calc_U_e(){
	var U_e = (window.gamma - 1) * me_kg * Math.pow(speed_c, 2) / elementary_charge * 1e-3; // [kV]
	document.getElementById('U_e').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + U_e.toFixed(3) + '</span>';
}

function LifeTimeREC(){
	/*
		estimate 1/e lifetime of the ion from Radioactive Electron Capture (REC)
		REC is the main mechanism of the beam loss under the ECooling

		REC@ H.Poth. Electron cooling: Theory, experiment, application. Physics Reports. 196 (1990): 135-297
    */
	let Q = Number(document.getElementById('ion_charge').value);
	let T_v = Number(document.getElementById('T_v').value);
	let e_gamma  = Number(document.getElementById('ion_gamma').value);
	let I_e = Number(document.getElementById('I_e').value);
	let r_0 = Number(document.getElementById('r_0').value);
	let len_C = Number(document.getElementById('len_C').value);
	let len_cool = Number(document.getElementById('len_cool').value);

	var alpha_rec = 3.02e-13 * Math.pow(Q, 2) / T_v * (Math.log(11.32 * Q / Math.pow(T_v, 1/2)) + 0.14 * Math.pow((T_v / Math.pow(Q, 2)), 1/3)); // [cm^3/s]
	var rho = I_e / Math.PI	/ Math.pow(r_0, 2) / elementary_charge / speed_c / Math.pow(1-1/Math.pow(e_gamma, 2), 1/2) * 1e-6; // [cm^-3]
	
	var tau_rec = Math.pow(e_gamma, 2) * len_C / len_cool / alpha_rec / rho; // [s]
	var result = document.getElementById('result_lifetime_rec');
	result.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_rec,3) + '</span>';
}


function LifeTimeES(){
	/*
		estimate 1/e lifetime of the ion from Elastic Scattering (ES)
	   	ES is one of the main mechanism of the beam loss with the reaction with the residual gas in the storage ring

	   	with assumption: round vacuum chamber, beta_x = beta_z = <beta_y>, varepsilon_x = varepsilon_z = varepsilon
        single-scattering@ H.Poth. Electron cooling: Theory, experiment, application. Physics Reports. 196 (1990): 135-297
        multiple-scattering@ B.Franzke. Interaction of stored ion beams with the residual gas. CAS-CERN Accelerator School: 4th advanced accelerator physics course. (1992): 100-119
	*/
	const k_b = 1.380649e-23; // [J/K] Boltzmann's constant
	const r_e = 2.8179403262e-15; // [m] classical electron radius
	const u = 1.66053906660e-27; // [kg] atomic mass unit

	let gas_pressure = Number(document.getElementById('gas_pressure').value); // [mbar]
	let gas_temperture = Number(document.getElementById('gas_temperture').value) + 273.15; // [K]
	var gas_rho = gas_pressure * 1e2 / k_b / gas_temperture; // [m^-3]

	let ion_Z = document.getElementById('ion_Z').value;
	var r_i = me_kg * r_e / window.ion_mass / u; // [m]
	
	var mode_index = document.getElementById('selectMode').selectedIndex;
	var A_nu = 0;
	var average_beta = 0;
	switch (mode_index) {
		case 0:
			A_nu = Number(document.getElementById('interTG_A_nu').value); // [πmm mrad]
			average_beta = Number(document.getElementById('interTG_beta').value); // [m]
			break;
		case 1:
			A_nu = Number(document.getElementById('isochronous_A_nu').value);
			average_beta = Number(document.getElementById('isochronous_beta').value);
			break;
		case 2:
			A_nu = Number(document.getElementById('normal_A_nu').value);
			average_beta = Number(document.getElementById('normal_beta').value);
			break;
	}
	
	let emittance = Number(document.getElementById('emittance').value); // [m]
	var sqr_angle_acceptance = A_nu * 1e-6 / average_beta * Math.pow(1 - emittance/A_nu, 2); // only value, no unit

	var lambda_ss_part = 0;
	var lambda_ms_part = 0;
	for (var i=0; i<6; i++) {
		let Z_t = Number(document.getElementById(`Gas${i}_Z_t`).value);
		let quantity = Number(document.getElementById(`Gas${i}_quantity_t`).value);
		var lambda_ss_temp = Math.pow(Z_t, 2) * quantity;  
		var lambda_ms_temp = Math.pow(Z_t, 2) * quantity * Math.log(204 / Math.pow(Z_t, 1/3)); 
		lambda_ss_part += lambda_ss_temp;
		lambda_ms_part += lambda_ms_temp;
	}
	var tau_ss = Math.pow(window.beta, 3) * Math.pow(window.gamma, 2) * sqr_angle_acceptance / 4 / Math.PI / Math.pow(ion_Z, 2) / Math.pow(r_i, 2) / speed_c / gas_rho / lambda_ss_part ;	
	var tau_ms = Math.pow(window.beta, 3) * Math.pow(window.gamma, 2) * (A_nu - emittance) * 1e-6 / average_beta / 32 / Math.pow(ion_Z, 2) / Math.pow(r_i, 2) / speed_c / gas_rho / lambda_ms_part;
	var tau_es = 1/(1/tau_ss + 1/tau_ms);

	var result_ss = document.getElementById('result_lifetime_ss');
	result_ss.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_ss,3) + '</span>';
	var result_ms = document.getElementById('result_lifetime_ms');
	result_ms.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_ms,3) + '</span>';
	var result_es = document.getElementById('result_lifetime_es');
	result_es.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_es,3) + '</span>';
}

function LifeTimeEC(){
	/*
		estimate 1/e lifetime of the ion from Electron Capture (EC)
	   	EC is one of the main mechanism of the beam loss in the ring.
	   	Ion's magnetic rigidity will be changed when capturing the electron from the residual gas.

	    default: this article gives a good result when Z_t < 18, ion_energy from 40-1000 MeV/u
	    single-EC@ I.S. Dmitriev, et al. On the target thickness to attain equilibrium charge distribution in a beam of fast ions. Nucl. Instr. Meth. Phys. Res. B. 14. (1986): 515-526
	   	option: this article gives a good result when ion_energy[keV/amu]/Z_t^1.25/ion_charge^0.7 from 10-1000, ion_charge >= 3
	   	EC@ A.S. Schlacher, et al. Electron capture for fast highly charged ions in gas targets: an empirical scaling rule. Phys. Rev. A. 27. (1983) 11: 3372
	   	option: this article gives a good result when ion_Z >= 36; gamma = (1 + T/931.5) <= 1.1
		single-EC@ B. Franzke. Vacuum requirements for heavy ion synchrotrons. IEEE Trans. Nucl. Sci. NS-28 (1981): 2116-2118
	*/			
	const k_b = 1.380649e-23; // [J/K] Boltzmann's constant

	let gas_pressure = Number(document.getElementById('gas_pressure').value); // [mbar]
	let gas_temperture = Number(document.getElementById('gas_temperture').value) + 273.15; // [K]
	var gas_rho = gas_pressure * 1e2 / k_b / gas_temperture; // [m^-3]

	let ion_Z = Number(document.getElementById('ion_Z').value);
	let ion_charge = Number(document.getElementById('ion_charge').value);
	
	var lambda_ec = 0;
	if ((window.gamma-1)/MeV2u*1e3/Math.pow(ion_charge, 0.7)/Math.pow(18,1.25) > 10 && (window.gamma-1)/MeV2u*1e3/Math.pow(ion_charge, 0.7) < 1000 && ion_charge >=3) {
	    document.getElementById('EC_method').innerHTML = "Schlacher Method";
		for (var i=0; i<6; i++){
			let Z_t = Number(document.getElementById(`Gas${i}_Z_t`).value);
			let quantity = Number(document.getElementById(`Gas${i}_quantity_t`).value);
			var rho = gas_rho * quantity; // [m^-3]
			var E_temp = (window.gamma-1)/MeV2u*1e3/Math.pow(ion_Z, 0.7)/Math.pow(Z_t,1.25); // [keV]
			var lambda_ec_temp = 1.1e-8/Math.pow(E_temp, 4.8) * (1 - Math.exp(-0.037 * Math.pow(E_temp, 2.2))) * (1 - Math.exp(-2.44e-5 * Math.pow(E_temp, 2.6))) * 1e-4 * rho * window.beta * speed_c * Math.pow(ion_charge, 0.5) / Math.pow(Z_t, 1.8); // [s^-1]
			lambda_ec += lambda_ec_temp;
		}
	}else if (ion_Z >= 36 && window.gamma <= 1.1){
	    document.getElementById('EC_method').innerHTML = "Franzke Method";
		var q_bar = ion_Z * (1 - Math.exp(-137*window.beta/Math.pow(ion_Z, 0.67)));
		if (ion_charge < q_bar){
			var b = 4;
		}else{
			var b = 2;
		}
		for (var i=0; i<6; i++){
			let Z_t = Number(document.getElementById(`Gas${i}_Z_t`).value);
			let quantity = Number(document.getElementById(`Gas${i}_quantity_t`).value);
			var rho = gas_rho * quantity; // [m^-3]
			var lambda_ec_temp = 2e-24 * 1e-4 * Math.pow(ion_Z, 1/2) * Math.pow(q_bar, 2) * Z_t * (1 - Math.exp(-137*window.beta/Math.pow(Z_t, 0.67))) * Math.pow(Math.pow(window.gamma, 2) - 1, -2) * Math.pow(ion_charge/q_bar, b) * rho * window.beta * c; // [s^-1]
			lambda_ec += lambda_ec_temp;
		}

	}else{
		var default_state = true;
		default_state *= (Number(document.getElementById('ion_energy_MeVu').value)>=40 && Number(document.getElementById('ion_energy_MeVu').value)<=1000);
		if (ion_charge <= 18){
			var K = 1;
		}else{
			var K = 1.20 - 0.01 * ion_Z;
		}
	
		for (var i=0; i<6; i++) {
			let Z_t = Number(document.getElementById(`Gas${i}_Z_t`).value);
			let quantity = Number(document.getElementById(`Gas${i}_quantity_t`).value);
			var rho = gas_rho * quantity;
			var lambda_ec_temp =  K * 1e-33 * Z_t * Math.pow(ion_Z, 5) * Math.pow(window.beta, -9) * 1e-4 * rho * speed_c;
			lambda_ec += lambda_ec_temp;
			default_state *= (Z_t < 18);
		}
		if (default_state){
	    	document.getElementById('EC_method').innerHTML = "Dmitriev Method";
		}else{
	    	document.getElementById('EC_method').innerHTML = "(Dmitriev Method)";
		}
	}
	var tau_ec = 1/lambda_ec;
	var result_ec = document.getElementById('result_lifetime_ec');
	result_ec.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_ec,3) + '</span>';
}

// binding Energy of ions from NIST Atomic Spectra Database Ionization Energies Data
// format: 'element': {Z-charge: ionization energy [eV]}
var bindEnergy = {'H':{0: 13.598434599702}, 'He': {0: 24.587389011, 1: 54.417765486}, 'Li':{3: 5.391714996, 2: 75.6400970, 1: 122.45435913}, 'Be': {3: 18.21115, 2: 153.896205, 1: 217.71858459}, 'B': {3: 37.93059, 2: 259.3715, 1: 340.2260225}, 'C': {3: 64.49352, 2: 392.090518, 1: 489.99320779}, 'N': {3: 97.8901, 2: 552.06733, 1: 667.0461377}, 'O': {3: 138.1189, 2: 739.32683, 1: 871.409913}, 'F': {3: 185.1868, 2: 953.89805, 1: 1103.1175302}, 'Ne': {3: 239.0970, 2: 1195.80784, 1: 1362.199256}, 'Na': {3: 299.856, 2: 1465.134502, 1: 1648.702285}, 'Mg': {3: 367.489, 2: 1761.80488, 1: 1962.663889}, 'Al': {3: 442.005, 2: 2085.97702, 1: 2304.140359}, 'Si': {3: 523.415, 2: 2437.65815, 1: 2673.177958}, 'P': {3: 611.741, 2: 2816.90879, 1: 3069.842145}, 'S': {3: 706.994, 2: 3223.7807, 1: 3494.188518}, 'Cl': {3: 809.198, 2: 3658.3438, 1: 3946.29179}, 'Ar': {3: 918.375, 2: 4120.6657, 1: 4426.22407}, 'K': {3: 1034.542, 2: 4610.87018, 1: 4934.04979}, 'Ca': {3: 1157.726, 2: 5128.8578, 1: 5469.86358}, 'Sc': {3: 1287.957, 2: 5674.9037, 1: 6033.75643}, 'Ti': {3: 1425.257, 2: 6249.0226, 1: 6625.81023}, 'V': {3: 1569.656, 2: 6851.3109, 1: 7246.12624}, 'Cr': {3: 1721.183, 2: 7481.8628, 1: 7894.80289}, 'Mn': {3: 1879.873, 2: 8140.7872, 1: 8571.95438}, 'Fe': {3: 2045.759, 2: 8828.1879, 1: 9277.6886}, 'Co': {3: 2218.876, 2: 9544.1833, 1: 10012.1297}, 'Ni': {3: 2399.259, 2: 10288.8862, 1: 10775.3948}, 'Cu': {3: 2586.954, 2: 11062.4313, 1: 11567.6237}, 'Zn': {3: 2781.996, 2: 11864.9399, 1: 12388.9427}, 'Ga': {3: 2984.426, 2: 12696.5575, 1: 13239.5029}, 'Ge': {3: 3194.293, 2: 13557.4208, 1: 14119.4457}, 'As': {3: 3411.643, 2: 14447.678, 1: 15028.9251}, 'Se': {3: 3636.526, 2: 15367.491, 1: 15968.1075}, 'Br': {3: 3868.986, 2: 16317.011, 1: 16937.1497}, 'Kr': {3: 4109.083, 2: 17296.424, 1: 17936.2405}, 'Rb': {3: 4356.865, 2: 18305.884, 1: 18965.5484}, 'Sr': {3: 4612.397, 2: 19345.588, 1: 20025.2673}, 'Y': {3: 4875.731, 2: 20415.717, 1: 21115.588}, 'Zr': {3: 5146.935, 2: 21516.469, 1: 22236.712}, 'Nb': {3: 5426.066, 2: 22648.046, 1: 23388.850}, 'Mo': {3: 5713.194, 2: 23810.654, 1: 24572.213}, 'Tc': {3: 6008.391, 2: 25004.533, 1: 25787.047}, 'Ru': {3: 6311.721, 2: 26229.895, 1: 27033.564}, 'Rh': {3: 6623.262, 2: 27486.983, 1: 28312.031}, 'Pd': {3: 6943.097, 2: 28776.034, 1: 29622.678}, 'Ag': {3: 7271.298, 2: 30097.318, 1: 30965.780}, 'Cd': {3: 7607.95, 2: 31451.062, 1: 32341.587}, 'In': {3: 7953.14, 2: 32837.592, 1: 33750.404}, 'Sn': {3: 8306.95, 2: 34257.143, 1: 35192.501}, 'Sb': {3: 8669.48, 2: 35710.028, 1: 36668.183}, 'Te': {3: 9040.83, 2: 37196.522, 1: 38177.740}, 'I': {3: 9421.10, 2: 38716.996, 1: 39721.549}, 'Xe': {3: 9810.37, 2: 40271.724, 1: 41299.892}, 'Cs': {3: 10208.78, 2: 41861.075, 1: 42913.144}, 'Ba': {3: 10616.42, 2: 43485.366, 1: 44561.633}, 'La': {3: 11033.40, 2: 45144.996, 1: 46245.77}, 'Ce': {3: 11459.85, 2: 46840.306, 1: 47965.89}, 'Pr': {3: 11895.89, 2: 48571.71, 1: 49722.44}, 'Nd': {3: 12341.66, 2: 50339.59, 1: 51515.78}, 'Pm': {3: 12797.26, 2: 52144.29, 1: 53346.31}, 'Sm': {3: 13262.85, 2: 53986.12, 1: 55214.30}, 'Eu': {3: 13738.58, 2: 55865.92, 1: 57120.64}, 'Gd': {3: 14224.57, 2: 57783.90, 1: 59065.54}, 'Tb': {3: 14721.02, 2: 59739.3, 1: 61050.1}, 'Dy': {3: 15228.06, 2: 61736.56, 1: 63073.23}, 'Ho': {3: 15745.77, 2: 63772.43, 1: 65137.13}, 'Er': {3: 16274.56, 2: 65848.24, 1: 67241.48}, 'Tm': {3: 16814.34, 2: 67965.26, 1: 69387.45}, 'Yb': {3: 17365.44, 2: 70123.04, 1: 71574.63}, 'Lu': {3: 17928.05, 2: 72322.91, 1: 73804.35}, 'Hf': {3: 18502.32, 2: 74565.93, 1: 76077.70}, 'Ta': {3: 19088.51, 2: 76852.03, 1: 78394.63}, 'W': {3: 19686.74, 2: 79181.94, 1: 80755.91}, 'Re': {3: 20297.40, 2: 81556.90, 1: 83162.41}, 'Os': {3: 20920.60, 2: 83976.21, 1: 85614.42}, 'Ir': {3: 21556.60, 2: 86438.9, 1: 88113.6}, 'Pt': {3: 22205.7, 2: 88955.18, 1: 90659.84}, 'Au': {3: 22868.1, 2: 91515.82, 1: 93254.62}, 'Hg': {3: 23544.1, 2: 94124.70, 1: 95898.19}, 'Tl': {3: 24234.1, 2: 96783.21, 1: 98592.12}, 'Pb': {3: 24938.2, 2: 99491.85, 1: 101336.7}, 'Bi': {3: 25656.9, 2: 102251.76, 1: 104133.4}, 'Po': {3: 26390.4, 2: 105064.3, 1: 106983.4}, 'At': {3: 27139.0, 2: 107923.4, 1: 109887.2}, 'Rn': {3: 27903.1, 2: 110842.0, 1: 112842.2}, 'Fr': {3: 28683.4, 2: 113817.2, 1: 115857.5}, 'Ra': {3: 29479.8, 2: 116848.7, 1: 118929.5}, 'Ac': {3: 30293.1, 2: 119938.6, 1: 122063.1}, 'Th': {3: 31122.8, 2: 123086.4, 1: 125250.3}, 'Pa': {3: 31971.6, 2: 126296.6, 1: 128507}, 'U': {3: 32836.5, 2: 129570.3, 1: 131816.2}, 'Np': {3: 33722.2, 2: 132901.8, 1: 135202}, 'Pu': {3: 34625.8, 2: 136299.2, 1: 138640.2}, 'Am': {3: 35549.4, 2: 139769.5, 1: 142153.5}, 'Cm': {3: 36493.0, 2: 143299.6, 1: 145740.1}, 'Bk': {3: 37457.6, 2: 146904.7, 1: 149398}, 'Cf': {3: 38443.5, 2: 150579.3, 1: 153124}, 'Es': {3: 39451.4, 2: 154328.1, 1: 156927}, 'Fm': {3: 40482.2, 2: 158152.5, 1: 160808}, 'Md': {3: 41548, 1: 164764}, 'No': {3: 42632, 1: 168804}, 'Lr': {3: 43759, 1: 172928}, 'Rf': {1: 177142}, 'Db': {1: 181445}, 'Sg': {1: 185835}, 'Bh': {1: 190329}, 'Hs': {1: 194911}, 'Mt': {1: 199605}, 'Ds': {1: 204394}};

function LifeTimeEL(){
	/*
		estimate 1/e lifetime of the ion from Electron Stripping (EL)
		EL is one of the main mechanism of the beam loss in the ring.
		EL is forbidden for the bare ion.
		with assumption: ion_beta / alpha > ion_Z

		default: ion velocity beta > Z/137
		EL@ D. Habs, et al. First experiments with the Heidelberg test storage ring TSR. Nucl. Instr. Meth. Phys. Res. B. 43. (1989): 390-410
		option: this article gives a good result when ion_Z >= 36; gamma = (1 + T/931.5) <= 1.1
		single-EL@ B. Franzke. Vacuum requirements for heavy ion synchrotrons. IEEE Trans. Nucl. Sci. NS-28 (1981): 2116-2118

	*/
	const k_b = 1.380649e-23; // [J/K] Boltzmann's constant
	const alpha_0 = 5.29177210903e-11; // [m], Bohr radius
	const alpha = 7.2973525693e-3; // fine-structure constant
	
	let gas_pressure = Number(document.getElementById('gas_pressure').value); // [mbar]
	let gas_temperture = Number(document.getElementById('gas_temperture').value) + 273.15; //[K]
	var gas_rho = gas_pressure * 1e2 / k_b / gas_temperture; // [m^-3]

	let ion_element = document.getElementById('ion_element').value;
	let ion_charge = Number(document.getElementById('ion_charge').value);
	let ion_Z = Number(document.getElementById('ion_Z').value);
	var ion_num = ion_Z - ion_charge;
	
	var lambda_el = 0;
	if (ion_Z >= 36 && window.gamma <= 1.1){
		document.getElementById('EL_method').innerHTML = "Franzke Method";
		var q_bar = ion_Z * (1 - Math.exp(-137*window.beta/Math.pow(ion_Z, 0.67)));
		if (ion_charge < q_bar){
			var b = -2.3;
		}else{
			var b = -4;
		}
		for (var i=0; i<6; i++){
			let Z_t = Number(document.getElementById(`Gas${i}_Z_t`).value);
			let quantity = Number(document.getElementById(`Gas${i}_quantity_t`).value);
			var rho = gas_rho * quantity; // [m^-3]
			var lambda_el_temp = 3.5 * Math.pow(10, -18+Math.pow(0.71 * Math.log(ion_charge), 1.5)) * 1e-4 * Math.pow(q_bar, -2) * Z_t * (1 - Math.exp(-137*window.beta/Math.pow(Z_t, 0.67))) * Math.pow(Math.pow(window.gamma, 2) - 1, -1/2) * Math.pow(ion_charge/q_bar, b) * rho * window.beta * speed_c; // [s^-1]
			lambda_el += lambda_el_temp;
		}	
	}else{
		document.getElementById('EL_method').innerHTML = "Habs Method";
		for (var i=0; i<6; i++){
			let Z_t = Number(document.getElementById(`Gas${i}_Z_t`).value);
			let quantity = Number(document.getElementById(`Gas${i}_quantity_t`).value);
			var rho = gas_rho * quantity;
			if (ion_Z >= Z_t){
				var I_ratio = 0;
				for (var j=1; j<ion_num+1; j++){
					I_ratio += bindEnergy['H'][0] / bindEnergy[ion_element][j];
				}
				var lambda_el_temp = 4 * Math.PI * Math.pow(alpha_0, 2) * Math.pow(alpha / window.beta, 2) * (Z_t + 1) * Z_t * I_ratio * rho * window.beta * speed_c; // [s^-1]
			}else if (ion_Z <= Math.pow(Z_t, 1/3)) {
				var I_ration = 0;
				for (var j=1; j<ion_num+1; j++){
					I_ratio += Math.pow(bindEnergy['H'][0] / bindEnergy[ion_element][j], 1/2);
				}
				var lambda_el_temp =  Math.PI * Math.pow(alpha_0, 2) * alpha / window.beta * Math.pow(Z_t, 2/3) * I_ratio * rho * window.beta * speed_c;
			}else{
				var lambda_el_temp = 0;
			}
			lambda_el += lambda_el_temp;
			console.log(lambda_el);
		}
	}
	var tau_el = 1/lambda_el;
	var result_el = document.getElementById('result_lifetime_el');
	result_el.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_el,3) + '</span>';
}

function LifeTimeNR(){
	/*
		estimate 1/e lifetime of the ion from Nuclear Reaction (NR)
		
		NR @ X.Y. Zhang. Study of Lifetime of Light Ion Beam stored in CSRm for Internal Target Experiment. (2005) MSc thesis 
	*/
	const k_b = 1.380649e-23; // [J/K] Boltzmann's constant

	let gas_pressure = Number(document.getElementById('gas_pressure').value); // [mbar]
	let gas_temperture = Number(document.getElementById('gas_temperture').value) + 273.15; // [K]
	var gas_rho = gas_pressure * 1e2 / k_b / gas_temperture; // [m^-3]

	var ion_A = Math.round(window.ion_mass);

	var lambda_nr = 0;
	for (var i=0; i<6; i++){
		let Z_t = Number(document.getElementById(`Gas${i}_Z_t`).value);
		let quantity = Number(document.getElementById(`Gas${i}_quantity_t`).value);
		var rho = gas_rho * quantity; // [m^-3]
		var lambda_nr_temp = Math.PI * Math.pow(Math.pow(ion_A, 1/3) + Math.pow(Z_t * 2, 1/3), 2) * 1e-30 / window.beta * rho * speed_c;
		lambda_nr += lambda_nr_temp;
	}
	var tau_nr = 1/lambda_nr;
	var result_nr = document.getElementById('result_lifetime_nr');
	result_nr.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_nr,3) + '</span>';
}

