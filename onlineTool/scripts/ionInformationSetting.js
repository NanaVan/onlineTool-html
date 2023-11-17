function toSci(xx, digi){
  //return xx.toExponential(digi).replace(/e\+?/, ' x 10^');
  return xx.toExponential(digi);
}

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

const elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds'];

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
	initialIon_info_dynamic();
	//document.getElementById('lifetimeElectronStripping_check').disabled = true;
}

function limiter_ion_Z(input){
	if (input.value > 110 || input.value < 4) input.value = '26';
	document.getElementById('ion_element').value = elements[Number(input.value)-1];
	document.getElementById('ion_charge').value = input.value;
	initialIon_info_dynamic();
	//document.getElementById('lifetimeElectronStripping_check').disabled = true;
}

function limiter_ion_A(input){
	initialIon_info_dynamic();
	//document.getElementById('lifetimeElectronStripping_check').disabled = true;
}

function limiter_ion_charge(input){
	if (input.value >= document.getElementById('ion_Z').value){
		input.value = document.getElementById('ion_Z').value;
		//document.getElementById('lifetimeElectronStripping_check').disabled = true;
	}else if (input.value < document.getElementById('ion_Z').value-3){
		input.value = document.getElementById('ion_Z').value-3;
		//document.getElementById('lifetimeElectronStripping_check').disabled = false;
	}else{
		//document.getElementById('lifetimeElectronStripping_check').disabled = false;
	}
	initialIon_info_dynamic();
}

// parameter from CODATA 2018
const speed_c = 299792458 // [m/s] speed of light in vacuum
const elementary_charge = 1.602176634e-19 // [C] elementary charge
const me = 5.48579909065e-4 // [u] electron mass in u
const u2kg = 1.66053906660e-27 // [kg/u] amount of kg per 1 u
const MeV2u = 1.07354410233e-3 // [u/MeV] amount of u per 1 MeV

async function initialIon_info_dynamic(){
	await checkIon();
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var Brho = 1; // [Tm]
	var gamma_beta = Brho / window.ion_mass * Q / speed_c / u2kg * elementary_charge;
	var beta = gamma_beta / Math.sqrt(1 + Math.pow(gamma_beta, 2));
	var velocity = beta * speed_c * 1e-7;
	var gamma = 1 / Math.sqrt(1 - Math.pow(beta, 2));
	var energy_MeVu = (gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho;
	document.getElementById('ion_gamma').value = gamma.toFixed(5);
	document.getElementById('ion_beta').value = beta.toFixed(5);
	document.getElementById('ion_velocity').value = velocity.toFixed(5);
} 


function change_ion_energy_MeVu(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var energy_MeVu = Number(input.value);
	var TKE = window.ion_mass * energy_MeVu;
	var energy_AMeV = TKE / A;
	var gamma = 1 + energy_MeVu * MeV2u;
	var beta = Math.sqrt(1 - 1/Math.pow(gamma,2));
	var velocity = beta * speed_c * 1e-7;
	var Brho = gamma * beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = gamma.toFixed(5);
	document.getElementById('ion_beta').value = beta.toFixed(5);
	document.getElementById('ion_velocity').value = velocity.toFixed(5);
}

function change_ion_energy_AMeV(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var energy_AMeV = Number(input.value);
	var TKE = A * energy_AMeV;
	var energy_MeVu = TKE / window.ion_mass;
	var gamma = 1 + energy_MeVu * MeV2u;
	var beta = Math.sqrt(1 - 1/Math.pow(gamma,2));
	var velocity = beta * speed_c * 1e-7;
	var Brho = gamma * beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = gamma.toFixed(5);
	document.getElementById('ion_beta').value = beta.toFixed(5);
	document.getElementById('ion_velocity').value = velocity.toFixed(5);
}

function change_ion_totalKineticEnergy(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var TKE = Number(input.value);
	var energy_MeVu = TKE / window.ion_mass;
	var energy_AMeV = TKE / A;
	var gamma = 1 + energy_MeVu * MeV2u;
	var beta = Math.sqrt(1 - 1/Math.pow(gamma,2));
	var velocity = beta * speed_c * 1e-7;
	var Brho = gamma * beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = gamma.toFixed(5);
	document.getElementById('ion_beta').value = beta.toFixed(5);
	document.getElementById('ion_velocity').value = velocity.toFixed(5);
}

function change_ion_Brho(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var Brho = Number(input.value);
	var gamma_beta = Brho / window.ion_mass * Q / speed_c / u2kg * elementary_charge;
	var beta = gamma_beta / Math.sqrt(1 + Math.pow(gamma_beta, 2));
	var velocity = beta * speed_c * 1e-7;
	var gamma = 1 / Math.sqrt(1 - Math.pow(beta, 2));
	var energy_MeVu = (gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_gamma').value = gamma.toFixed(5);
	document.getElementById('ion_beta').value = beta.toFixed(5);
	document.getElementById('ion_velocity').value = velocity.toFixed(5);
}

function change_ion_gamma(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var gamma = Number(input.value);
	var beta = Math.sqrt(1 - 1/Math.pow(gamma,2));
	var velocity = beta * speed_c * 1e-7;
	var Brho = gamma * beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	var energy_MeVu = (gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_beta').value = beta.toFixed(5);
	document.getElementById('ion_velocity').value = velocity.toFixed(5);
}

function change_ion_beta(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var beta = Number(input.value);
	var gamma = 1 / Math.sqrt(1 - Math.pow(beta, 2));
	var velocity = beta * speed_c * 1e-7;
	var Brho = gamma * beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	var energy_MeVu = (gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = gamma.toFixed(5);
	document.getElementById('ion_velocity').value = velocity.toFixed(5);
}

function change_ion_velocity(input){
	var Q = Number(document.getElementById('ion_charge').value);
	var A = Number(document.getElementById('ion_A').value);
	var velocity = Number(input.value);
	var beta = velocity * 1e7 / speed_c;
	var gamma = 1 / Math.sqrt(1 - Math.pow(beta, 2));
	var Brho = gamma * beta * window.ion_mass / Q * speed_c * u2kg / elementary_charge;
	var energy_MeVu = (gamma - 1) / MeV2u;
	var TKE = window.ion_mass * energy_MeVu; 
	var energy_AMeV = TKE / A;
	document.getElementById('ion_energy_MeVu').value = energy_MeVu.toFixed(5);
	document.getElementById('ion_energy_AMeV').value = energy_AMeV.toFixed(5);
	document.getElementById('ion_totalKineticEnergy').value = TKE.toFixed(5);
	document.getElementById('ion_Brho').value = Brho.toFixed(5);
	document.getElementById('ion_gamma').value = gamma.toFixed(5);
	document.getElementById('ion_beta').value = beta.toFixed(5);
}
