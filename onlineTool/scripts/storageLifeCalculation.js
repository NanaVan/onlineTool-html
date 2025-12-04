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
			A_nu = Number(document.getElementById('interTG_befDec_A_nu').value); // [πmm mrad]
			average_beta = Number(document.getElementById('interTG_befDec_beta').value); // [m]
			break;
		case 1:
			A_nu = Number(document.getElementById('interTG_AfDec_A_nu').value); // [πmm mrad]
			average_beta = Number(document.getElementById('interTG_AfDec_beta').value); // [m]
			break;
		case 2:
			A_nu = Number(document.getElementById('isochronous_1.43_A_nu').value);
			average_beta = Number(document.getElementById('isochronous_1.43_beta').value);
			break;
		case 3:
			A_nu = Number(document.getElementById('isochronous_1.67_A_nu').value);
			average_beta = Number(document.getElementById('isochronous_1.67_beta').value);
			break;
		case 4:
			A_nu = Number(document.getElementById('normal_large_A_nu').value);
			average_beta = Number(document.getElementById('normal_large_beta').value);
			break;
		case 5:
			A_nu = Number(document.getElementById('normal_small_A_nu').value);
			average_beta = Number(document.getElementById('normal_small_beta').value);
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
			var lambda_ec_temp = 2e-24 * 1e-4 * Math.pow(ion_Z, 1/2) * Math.pow(q_bar, 2) * Z_t * (1 - Math.exp(-137*window.beta/Math.pow(Z_t, 0.67))) * Math.pow(Math.pow(window.gamma, 2) - 1, -2) * Math.pow(ion_charge/q_bar, b) * rho * window.beta * speed_c; // [s^-1]
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
					I_ratio += bindEnergy['H'][1] / bindEnergy[ion_element][j];
				}
				var lambda_el_temp = 4 * Math.PI * Math.pow(alpha_0, 2) * Math.pow(alpha / window.beta, 2) * (Z_t + 1) * Z_t * I_ratio * rho * window.beta * speed_c; // [s^-1]
			}else if (ion_Z <= Math.pow(Z_t, 1/3)) {
				var I_ration = 0;
				for (var j=1; j<ion_num+1; j++){
					I_ratio += Math.pow(bindEnergy['H'][1] / bindEnergy[ion_element][j], 1/2);
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




// controller

function LifeTimeTotal(status){
	if (status){
		var lambda_nr = 0;
		var lambda_REC = 0;
		var lambda_ES = 0;
		var lambda_EC = 0;
		var lambda_EL = 0;
		var lambda_NR = 0;
		if (document.getElementById('lifetimeREC_check').checked == true){
			lambda_REC = 1/Number(document.getElementById('result_lifetime_rec').textContent);
		}
		if (document.getElementById('lifetimeScattering_check').checked == true){
			lambda_ES = 1/Number(document.getElementById('result_lifetime_es').textContent);
		}
		if (document.getElementById('lifetimeElectronCapture_check').checked == true){
			lambda_EC = 1/Number(document.getElementById('result_lifetime_ec').textContent);
		}
		if (document.getElementById('lifetimeElectronStripping_check').checked == true){
			lambda_EL = 1/Number(document.getElementById('result_lifetime_el').textContent);
		}
		if (document.getElementById('lifetimeNuclearReaction_check').checked == true){
			lambda_NR = 1/Number(document.getElementById('result_lifetime_nr').textContent);
		}
		lambda_total = lambda_REC + lambda_ES + lambda_EC + lambda_EL + lambda_NR;
		var tau_total = 1/lambda_total;
		var result_total_s = document.getElementById('result_lifetime_total_s');
		result_total_s.innerHTML = '<span style="color:red; font-family:Roboto;">' + toSci(tau_total,3) + '</span>';
		var result_total_min = document.getElementById('result_lifetime_total_min');
		result_total_min.innerHTML = '<span style="color:red; font-family:Roboto;">' + toSci(tau_total/60,3) + '</span>';

		var result_total_hour = document.getElementById('result_lifetime_total_hour');
		result_total_hour.innerHTML = '<span style="color:red; font-family:Roboto;">' + toSci(tau_total/3600,3) + '</span>';

		document.getElementById('ratio_lifetime_REC').innerHTML = '<span style="color:red; font-family:Roboto;">' + (lambda_REC/lambda_total*100).toFixed(3) + '</span>';
		document.getElementById('ratio_lifetime_ES').innerHTML = '<span style="color:red; font-family:Roboto;">' + (lambda_ES/lambda_total*100).toFixed(3) + '</span>';
		document.getElementById('ratio_lifetime_EC').innerHTML = '<span style="color:red; font-family:Roboto;">' + (lambda_EC/lambda_total*100).toFixed(3) + '</span>';
		document.getElementById('ratio_lifetime_EL').innerHTML = '<span style="color:red; font-family:Roboto;">' + (lambda_EL/lambda_total*100).toFixed(3) + '</span>';
		document.getElementById('ratio_lifetime_NR').innerHTML = '<span style="color:red; font-family:Roboto;">' + (lambda_NR/lambda_total*100).toFixed(3) + '</span>';

	}else{
		var result_total_s = document.getElementById('result_lifetime_total_s');
		result_total_s.innerHTML = "";
		var result_total_min = document.getElementById('result_lifetime_total_min');
		result_total_min.innerHTML = "";
		var result_total_hour = document.getElementById('result_lifetime_total_hour');
		result_total_hour.innerHTML = "";
		document.getElementById('ratio_lifetime_REC').innerHTML = "";
		document.getElementById('ratio_lifetime_ES').innerHTML = "";
		document.getElementById('ratio_lifetime_EC').innerHTML = "";
		document.getElementById('ratio_lifetime_EL').innerHTML = "";
		document.getElementById('ratio_lifetime_NR').innerHTML = "";	
	}
}

function InformationLock(){
	if (document.getElementById('lifetimeREC_check').checked == true || document.getElementById('lifetimeScattering_check').checked == true || document.getElementById('lifetimeElectronCapture_check').checked == true || document.getElementById('lifetimeElectronStripping_check').checked == true || document.getElementById('lifetimeNuclearReaction_check').checked == true) {
		SettingInformationStatus(true);					
		RingModeInformationStatus(true);
		GasInformationStatus(true);
		LifeTimeTotal(true);
	}else{
		SettingInformationStatus(false);					
		RingModeInformationStatus(false);
		GasInformationStatus(false);
		LifeTimeTotal(false);
	}
}

function SettingInformationStatus(status){
	if (status){
		document.getElementById('ion_element').disabled = true;
		document.getElementById('ion_Z').disabled = true;
		document.getElementById('ion_A').disabled = true;
		document.getElementById('ion_charge').disabled = true;
		document.getElementById('ion_energy_MeVu').disabled = true;
		document.getElementById('ion_energy_AMeV').disabled = true;
		document.getElementById('ion_totalKineticEnergy').disabled = true;
		document.getElementById('ion_Brho').disabled = true;
		document.getElementById('ion_gamma').disabled = true;
		document.getElementById('ion_beta').disabled = true;
		document.getElementById('ion_velocity').disabled = true;
	}else{
		document.getElementById('ion_element').disabled = false;
		document.getElementById('ion_Z').disabled = false;
		document.getElementById('ion_A').disabled = false;
		document.getElementById('ion_charge').disabled = false;
		document.getElementById('ion_energy_MeVu').disabled = false;
		document.getElementById('ion_energy_AMeV').disabled = false;
		document.getElementById('ion_totalKineticEnergy').disabled = false;
		document.getElementById('ion_Brho').disabled = false;
		document.getElementById('ion_gamma').disabled = false;
		document.getElementById('ion_beta').disabled = false;
		document.getElementById('ion_velocity').disabled = false;
	}
}

function RingModeInformationStatus(status){
	if (status){
		document.getElementById('interTG_A_nu').disabled = true;
		document.getElementById('interTG_beta').disabled = true;
		document.getElementById('isochronous_A_nu').disabled = true;
		document.getElementById('isochronous_beta').disabled = true;
		document.getElementById('normal_A_nu').disabled = true;
		document.getElementById('normal_beta').disabled = true;
	}else{
		document.getElementById('interTG_A_nu').disabled = false;
		document.getElementById('interTG_beta').disabled = false;
		document.getElementById('isochronous_A_nu').disabled = false;
		document.getElementById('isochronous_beta').disabled = false;
		document.getElementById('normal_A_nu').disabled = false;
		document.getElementById('normal_beta').disabled = false;
	}
}

function GasInformationStatus(status){
	if (status){
		document.getElementById('gas_pressure').disabled = true;
		document.getElementById('gas_temperture').disabled = true;
		for (var i=0; i<6; i++){
			document.getElementById(`Gas${i}_Z_t`).disabled = true;
			document.getElementById(`Gas${i}_quantity_t`).disabled = true;
		}
	}else{
		document.getElementById('gas_pressure').disabled = false;
		document.getElementById('gas_temperture').disabled = false;
		for (var i=0; i<6; i++){
			document.getElementById(`Gas${i}_Z_t`).disabled = false;
			document.getElementById(`Gas${i}_quantity_t`).disabled = false;
		}
	}
}

function selectREC(e){
	if (e.checked == true){
		console.log('flag REC 1');
		document.getElementById('T_v').disabled = true;
		document.getElementById('I_e').disabled = true;
		document.getElementById('r_0').disabled = true;
		document.getElementById('len_C').disabled = true;
		document.getElementById('len_cool').disabled = true;
		LifeTimeREC();
	}else{
		console.log('flag REC 0');
		var result = document.getElementById('result_lifetime_rec');
		result.innerHTML = "";
		document.getElementById('T_v').disabled = false;
		document.getElementById('I_e').disabled = false;
		document.getElementById('r_0').disabled = false;
		document.getElementById('len_C').disabled = false;
		document.getElementById('len_cool').disabled = false;
	}
	InformationLock();
}

function selectES(e){
	if (e.checked == true){
		console.log('flag ES 1');
		document.getElementById('selectMode').disabled = true;
		document.getElementById('emittance').disabled = true;
		LifeTimeES();
	}else{
		console.log('flag ES 0');
		var result_ss = document.getElementById('result_lifetime_ss');
		result_lifetime_ss.innerHTML = "";
		var result_ms = document.getElementById('result_lifetime_ms');
		result_lifetime_ms.innerHTML = "";
		var result_es = document.getElementById('result_lifetime_es');
		result_lifetime_es.innerHTML = "";
		document.getElementById('selectMode').disabled = false;
		document.getElementById('emittance').disabled = false;
	}
	InformationLock();
}

function selectEC(e){
	if (e.checked == true){
		console.log('flag EC 1');
		LifeTimeEC();
	}else{
		console.log('flag EC 0');
		var result_ec = document.getElementById('result_lifetime_ec');
		result_ec.innerHTML = "";
		var method_ec = document.getElementById('EC_method');
		method_ec.innerHTML = "";
	}
	InformationLock();
}

function selectEL(e){
	if (e.checked == true){
		console.log('flag EL 1');
		LifeTimeEL();
	}else{
		console.log('flag EL 0');
		var result_el = document.getElementById('result_lifetime_el');
		result_el.innerHTML = "";
		var method_ec = document.getElementById('EL_method');
		method_ec.innerHTML = "";
	}
	InformationLock();
}

function selectNR(e){
	if (e.checked == true){
		console.log('flag NR 1');
		LifeTimeNR();
	}else{
		console.log('flag NR 0');
		var result_el = document.getElementById('result_lifetime_nr');
		result_el.innerHTML = "";
	}
	InformationLock();
}

