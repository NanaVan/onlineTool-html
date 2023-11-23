async function checkIon(){
	const [SQL, buf] = await Promise.all([
		initSqlJs({locateFile: file=> `/dist/${file}`}),
		fetch("/dist/ionic_data.sqlite").then(res => res.arrayBuffer())
		// ionic_data is from NUBASE2020
	]);		
	const db = new SQL.Database(new Uint8Array(buf));
	let stmt = db.prepare("SELECT min(A), max(A) FROM IONICDATA WHERE Z=$zval AND ISOMERIC=0");
	let result = stmt.get({'$zval': document.getElementById('ion_Z').value});
	//console.log(result);
	if (document.getElementById('ion_A').value < result[0] || document.getElementById('ion_A').value > result[1]){
		console.log('out of range');
		document.getElementById('ion_A').value = result[0];
	}
	stmt.free();
	stmt = db.prepare("SELECT MASS, TYPE, HALFLIFE FROM IONICDATA WHERE Z=$zval AND A=$aval AND Q=$qval AND ISOMERIC=0");
	result = stmt.get({'$zval': document.getElementById('ion_Z').value, '$aval': document.getElementById('ion_A').value, '$qval': document.getElementById('ion_charge').value});
	//console.log(result);
	document.getElementById('ion_info_static').rows[4].cells[1].innerHTML = Number(result[0]).toFixed(5);
    document.getElementById('ion_info_static').rows[5].cells[1].innerHTML = result[1];
    document.getElementById('ion_info_static').rows[6].cells[1].innerHTML = result[2];
	window.ion_mass = result[0];
	db.close();
}

function limiter_ion_element(input){
	var index = elements.indexOf(input.value);
	console.log(index);
	if (index == -1 || index < 3){
		index = Number(document.getElementById('ion_Z').value)-1;
		input.value = elements[index];
	} else {
		document.getElementById('ion_Z').value = String(index + 1);
	}
	document.getElementById('ion_charge').value = document.getElementById('ion_Z').value;
	document.getElementById('lifetimeElectronStripping_check').disabled = true;
	initialIon_info_dynamic();
}

function limiter_ion_Z(input){
	if (input.value > 110 || input.value < 4) input.value = '26';
	document.getElementById('ion_element').value = elements[Number(input.value)-1];
	document.getElementById('ion_charge').value = input.value;
	document.getElementById('lifetimeElectronStripping_check').disabled = true;
	initialIon_info_dynamic();
}

function limiter_ion_A(input){
	document.getElementById('lifetimeElectronStripping_check').disabled = true;
	initialIon_info_dynamic();
}

function limiter_ion_charge(input){
	if (input.value >= document.getElementById('ion_Z').value){
		input.value = document.getElementById('ion_Z').value;
		document.getElementById('lifetimeElectronStripping_check').disabled = true;
	}else if (input.value < document.getElementById('ion_Z').value-3){
		input.value = document.getElementById('ion_Z').value-3;
		document.getElementById('lifetimeElectronStripping_check').disabled = false;
	}else{
		document.getElementById('lifetimeElectronStripping_check').disabled = false;
	}
	initialIon_info_dynamic();
}


async function initialIon_info_dynamic(){
	await checkIon();
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var Brho = 1; // [Tm]
	var gamma_beta = Brho / window.ion_mass * Q / speed_c / u2kg * elementary_charge;
	window.beta = gamma_beta / Math.sqrt(1 + Math.pow(gamma_beta, 2));
	window.velocity = window.beta * speed_c * 1e-7;
	window.gamma = 1 / Math.sqrt(1 - Math.pow(window.beta, 2));
	var energy_MeVu = (window.gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho;
	document.getElementById('ion_gamma').value = window.gamma.toFixed(5);
	document.getElementById('ion_beta').value = window.beta.toFixed(5);
	document.getElementById('ion_velocity').value = window.velocity.toFixed(5);
	harmCalc();
	calc_U_e();
} 


function change_ion_energy_MeVu(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var energy_MeVu = Number(input.value);
	var TKE = window.ion_mass * energy_MeVu;
	var energy_AMeV = TKE / A;
	window.gamma = 1 + energy_MeVu * MeV2u;
	window.beta = Math.sqrt(1 - 1/Math.pow(window.gamma,2));
	window.velocity = window.beta * speed_c * 1e-7;
	var Brho = window.gamma * window.beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = window.gamma.toFixed(5);
	document.getElementById('ion_beta').value = window.beta.toFixed(5);
	document.getElementById('ion_velocity').value = window.velocity.toFixed(5);
	harmCalc();
	calc_U_e();
}

function change_ion_energy_AMeV(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var energy_AMeV = Number(input.value);
	var TKE = A * energy_AMeV;
	var energy_MeVu = TKE / window.ion_mass;
	window.gamma = 1 + energy_MeVu * MeV2u;
	window.beta = Math.sqrt(1 - 1/Math.pow(window.gamma,2));
	window.velocity = window.beta * speed_c * 1e-7;
	var Brho = window.gamma * window.beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = window.gamma.toFixed(5);
	document.getElementById('ion_beta').value = window.beta.toFixed(5);
	document.getElementById('ion_velocity').value = window.velocity.toFixed(5);
	harmCalc();
	calc_U_e();
}

function change_ion_totalKineticEnergy(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var TKE = Number(input.value);
	var energy_MeVu = TKE / window.ion_mass;
	var energy_AMeV = TKE / A;
	window.gamma = 1 + energy_MeVu * MeV2u;
	window.beta = Math.sqrt(1 - 1/Math.pow(window.gamma,2));
	window.velocity = window.beta * speed_c * 1e-7;
	var Brho = window.gamma * window.beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = window.gamma.toFixed(5);
	document.getElementById('ion_beta').value = window.beta.toFixed(5);
	document.getElementById('ion_velocity').value = window.velocity.toFixed(5);
	harmCalc();
	calc_U_e();
}

function change_ion_Brho(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var Brho = Number(input.value);
	var gamma_beta = Brho / window.ion_mass * Q / speed_c / u2kg * elementary_charge;
	window.beta = gamma_beta / Math.sqrt(1 + Math.pow(gamma_beta, 2));
	window.velocity = window.beta * speed_c * 1e-7;
	window.gamma = 1 / Math.sqrt(1 - Math.pow(window.beta, 2));
	var energy_MeVu = (window.gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_gamma').value = window.gamma.toFixed(5);
	document.getElementById('ion_beta').value = window.beta.toFixed(5);
	document.getElementById('ion_velocity').value = window.velocity.toFixed(5);
	harmCalc();
	calc_U_e();
}

function change_ion_gamma(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	window.gamma = Number(input.value);
	window.beta = Math.sqrt(1 - 1/Math.pow(window.gamma,2));
	window.velocity = window.beta * speed_c * 1e-7;
	var Brho = window.gamma * window.beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	var energy_MeVu = (window.gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_beta').value = window.beta.toFixed(5);
	document.getElementById('ion_velocity').value = window.velocity.toFixed(5);
	harmCalc();
	calc_U_e();
}

function change_ion_beta(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	window.beta = Number(input.value);
	window.gamma = 1 / Math.sqrt(1 - Math.pow(window.beta, 2));
	window.velocity = window.beta * speed_c * 1e-7;
	var Brho = window.gamma * window.beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	var energy_MeVu = (window.gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = window.gamma.toFixed(5);
	document.getElementById('ion_velocity').value = window.velocity.toFixed(5);
	harmCalc();
	calc_U_e();
}

function change_ion_velocity(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	window.velocity = Number(input.value);
	window.beta = window.velocity * 1e7 / speed_c;
	window.gamma = 1 / Math.sqrt(1 - Math.pow(window.beta, 2));
	var Brho = window.gamma * window.beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	var energy_MeVu = (window.gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = window.gamma.toFixed(5);
	document.getElementById('ion_beta').value = window.beta.toFixed(5);
	harmCalc();
	calc_U_e();
}
