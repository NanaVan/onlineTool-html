function toSci(xx, digi){
  //return xx.toExponential(digi).replace(/e\+?/, ' x 10^');
  return xx.toExponential(digi);
}

checkIon();

async function checkIon(){
	const [SQL, buf] = await Promise.all([
		initSqlJs({locateFile: file=> `/dist/${file}`}),
		fetch("/dist/ionic_data.sqlite").then(res => res.arrayBuffer())
	]);		
	const db = new SQL.Database(new Uint8Array(buf));
	let stmt = db.prepare("SELECT min(A), max(A) FROM IONICDATA WHERE Z=$zval AND ISOMERIC=0");
	let result = stmt.get({'$zval': document.getElementById('ion_Z').value});
	console.log(result);
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
	checkIon();
	//document.getElementById('lifetimeElectronStripping_check').disabled = true;
}

function limiter_ion_Z(input){
	if (input.value > 110 || input.value < 4) input.value = '26';
	document.getElementById('ion_element').value = elements[Number(input.value)-1];
	document.getElementById('ion_charge').value = input.value;
	checkIon();
	//document.getElementById('lifetimeElectronStripping_check').disabled = true;
}

function limiter_ion_A(input){
	console.log('A check');
	checkIon();
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
	checkIon();
}

// 
