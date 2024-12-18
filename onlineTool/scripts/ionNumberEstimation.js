async function beamCalc(case_undefault){
	const [SQL, buf] =  await Promise.all([
		initSqlJs({locateFile: file=> `/dist/${file}`}),
		fetch("/dist/ionic_data.sqlite").then(res => res.arrayBuffer())
	]);	
	const db = new SQL.Database(new Uint8Array(buf));
	let stmt = db.prepare("SELECT MASS FROM IONICDATA WHERE Z=$zval AND Q=$qval AND A=$aval AND ISOMERIC=0");
	let primaryBeam_result = stmt.get({'$zval': elements.indexOf(document.getElementById('primaryBeam_element').value)+1, '$qval': document.getElementById('primaryBeam_Q').value, '$aval': document.getElementById('primaryBeam_A').value});
	let secondaryBeam_result = stmt.get({'$zval': elements.indexOf(document.getElementById('secondaryBeam_element').value)+1, '$qval': document.getElementById('secondaryBeam_Q').value, '$aval': document.getElementById('secondaryBeam_A').value});
	stmt.free();
	db.close();
	if (!primaryBeam_result.length){
  		document.getElementById('primaryBeam').rows[7].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + 'No such ion!' + '</span>';
	}else{
		let primaryBeam_Q = document.getElementById('primaryBeam_Q').value;
  		let DCCTm = document.getElementById('DCCTm').value; // [μA]
  		let LengthCSRm = document.getElementById('CSRmL').value; // [m]
		let primaryBeam_mass = primaryBeam_result[0]; // [MeV]
		if (case_undefault){
			let primaryBeam_Brho = document.getElementById('primaryBeam_Brho').value; // [Tm]
			let primaryBeam_gamma_beta = primaryBeam_Brho / primaryBeam_mass * primaryBeam_Q / speed_c / u2kg * elementary_charge;
			var primaryBeam_beta = primaryBeam_gamma_beta / Math.sqrt(1 + Math.pow(primaryBeam_gamma_beta, 2));
			let primaryBeam_gamma = 1 / Math.sqrt(1 - Math.pow(primaryBeam_beta, 2));
			var primaryBeam_energy = (primaryBeam_gamma - 1) / MeV2u;
			document.getElementById('primaryBeam_energy').value = primaryBeam_energy.toFixed(5);
		}else{
			let primaryBeam_energy = document.getElementById('primaryBeam_energy').value; // [MeV/u]
			let primaryBeam_gamma = 1 + primaryBeam_energy * MeV2u;
			var primaryBeam_beta = Math.sqrt(1 - 1/Math.pow(primaryBeam_gamma,2));
			var primaryBeam_Brho = primaryBeam_gamma * primaryBeam_beta * primaryBeam_mass / primaryBeam_Q * speed_c * u2kg / elementary_charge; // [Tm]
			document.getElementById('primaryBeam_Brho').value = primaryBeam_Brho.toFixed(5);
		}
		var PBppp = DCCTm*1e-6 / primaryBeam_Q / (primaryBeam_beta*speed_c) * LengthCSRm / elementary_charge;
  		document.getElementById('primaryBeam').rows[7].cells[1].innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(PBppp,3) + '</span>';
	}	
	if (!secondaryBeam_result.length){
  		document.getElementById('secondaryBeam').rows[7].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + 'No such ion!' + '</span>';
	}else{
		let secondaryBeam_Q = document.getElementById('secondaryBeam_Q').value;
  		let DCCTe = document.getElementById('DCCTe').value;  // [μA]
  		let LengthCSRe = document.getElementById('CSReL').value; // [m]
		let secondaryBeam_mass = secondaryBeam_result[0]; // [MeV]
		if (case_undefault){
			let secondaryBeam_energy = document.getElementById('secondaryBeam_energy').value; // [MeV/u]
			let secondaryBeam_gamma = 1 + secondaryBeam_energy * MeV2u;
			var secondaryBeam_beta = Math.sqrt(1 - 1/Math.pow(secondaryBeam_gamma,2));
			var secondaryBeam_Brho = secondaryBeam_gamma * secondaryBeam_beta * secondaryBeam_mass / secondaryBeam_Q * speed_c * u2kg / elementary_charge; // [Tm]
			document.getElementById('secondaryBeam_Brho').value = secondaryBeam_Brho.toFixed(5);
		}else{
			let secondaryBeam_Brho = document.getElementById('secondaryBeam_Brho').value; // [Tm]
			let secondaryBeam_gamma_beta = secondaryBeam_Brho / secondaryBeam_mass * secondaryBeam_Q / speed_c / u2kg * elementary_charge;
			var secondaryBeam_beta = secondaryBeam_gamma_beta / Math.sqrt(1 + Math.pow(secondaryBeam_gamma_beta, 2));
			let secondaryBeam_gamma = 1 / Math.sqrt(1 - Math.pow(secondaryBeam_beta, 2));
			var secondaryBeam_energy = (secondaryBeam_gamma - 1) / MeV2u;
			document.getElementById('secondaryBeam_energy').value = secondaryBeam_energy.toFixed(5);
		}
		var SBppp = DCCTe*1e-6 / secondaryBeam_Q / (secondaryBeam_beta*speed_c) * LengthCSRe / elementary_charge;
  		document.getElementById('secondaryBeam').rows[7].cells[1].innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(SBppp,3) + '</span>';
	}
	if (primaryBeam_result.length && secondaryBeam_result.length){
  		var TE = Math.round(SBppp/PBppp * 10000) / 100.0;
  		document.getElementById('convertor').rows[1].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + TE + "%" + '</span>';
		
	}else{
  		document.getElementById('convertor').rows[1].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + "Error" + '</span>';
	}
}


function limiter_primaryBeam_element(input){
	var index = elements.indexOf(input.value);
	console.log(index);
	if (index == -1 || index < 3){
		index = Number('Fe')-1;
		input.value = elements[index];
	}
}

function limiter_secondaryBeam_element(input){
	var index = elements.indexOf(input.value);
	console.log(index);
	if (index == -1 || index < 3){
		index = Number('Fe')-1;
		input.value = elements[index];
	}
}

document.getElementById('primaryBeam_element').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13 ){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('primaryBeam_A').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13 ){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('primaryBeam_Q').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13 ){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('DCCTm').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('CSRmL').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('primaryBeam_energy').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('primaryBeam_Brho').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(true);      
    }
  }, false
);

document.getElementById('secondaryBeam_element').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('secondaryBeam_A').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('secondaryBeam_Q').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('DCCTe').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('CSReL').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);
document.getElementById('secondaryBeam_energy').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(true);      
    }
  }, false
);
document.getElementById('secondaryBeam_Brho').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      beamCalc(false);      
    }
  }, false
);

