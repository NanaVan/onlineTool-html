async function initialHarm_calc(){
	await initialIon_info_dynamic();
	harmCalc();
}

function harmCalc(){
	let velocity = document.getElementById('ion_velocity').value; // [cm/ns]
	let lenOrbit = document.getElementById('length_orbit').value; // [m]
	let cenFreq = document.getElementById('center_frequency').value; // [MHz]
	let span = document.getElementById('span').value; // [kHz]

	var frequency = velocity / lenOrbit * 1e4; // [kHz]
	document.getElementById('result_rev_freq').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + (frequency).toFixed(5) + '</span>';
	document.getElementById('result_rev_time').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + (1e6/frequency).toFixed(6) + '</span>';
	var harm_start = Math.ceil((cenFreq*1e3-span/2)/frequency);
	var harm_until = Math.floor((cenFreq*1e3+span/2)/frequency);
	
	const cols = 3;
	var harmTable = document.getElementById('harmResult');
	if (harm_until >= harm_start) {
		var rows = harm_until - harm_start + 1;
		var harmTBody = document.createElement('tbody');
		for (var i=0; i<rows; i++){
			var row = document.createElement('tr');
			var harmLabel = new Array(harm_start+i, ((harm_start+i)*frequency-cenFreq*1e3).toFixed(6), ((harm_start+i)*frequency*1e-3).toFixed(9));
			for (let j=0; j<cols; j++) {
				var cell = document.createElement('td');
				var cellText = document.createTextNode(harmLabel[j]);
				cell.appendChild(cellText);
				row.appendChild(cell);
			}
			harmTBody.appendChild(row);
		}
	}else{
		var harmTBody = document.createElement('tbody');
		var row = document.createElement('tr')
		var harmLabel = new Array("", "", "");
		for (let j=0; j<cols; j++) {
			var cell = document.createElement('td');
			var cellText = document.createTextNode(harmLabel[j]);
			cell.appendChild(cellText);
			row.appendChild(cell);
		}
		harmTBody.appendChild(row);
	}
	var old_tbody = document.getElementById('old_tbody');
	old_tbody.parentNode.replaceChild(harmTBody, old_tbody);
	harmTBody.setAttribute('id', 'old_tbody');
}

document.getElementById('ion_velocity').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      harmCalc();      
    }
  }, false
);

document.getElementById('length_orbit').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      harmCalc();      
    }
  }, false
);

document.getElementById('center_frequency').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      harmCalc();      
    }
  }, false
);

document.getElementById('span').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      harmCalc();      
    }
  }, false
);

