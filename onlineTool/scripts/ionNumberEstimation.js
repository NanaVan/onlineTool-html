function toSci(xx, digi){
  //return xx.toExponential(digi).replace(/e\+?/, ' x 10^');
  return xx.toExponential(digi);
}

function BeamCalc(){
  let eee = 624151;
  let chargeStateP = document.getElementById('ChargeStateP').value; 
  let DCCTm = document.getElementById('DCCTm').value; 
  let LengthCSRm = document.getElementById('CSRmL').value;
  let PBv = document.getElementById('PBv').value;

  var PBppp = DCCTm/chargeStateP/PBv*LengthCSRm*eee;

  let chargeStateS = document.getElementById('ChargeStateS').value; 
  let DCCTe = document.getElementById('DCCTe').value; 
  let LengthCSRe = document.getElementById('CSReL').value;
  let SBv = document.getElementById('SBv').value;

  var SBppp = DCCTe/chargeStateS/SBv*LengthCSRe*eee;

  document.getElementById('primaryBeam').rows[5].cells[1].innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(PBppp,3) + '</span>';
  document.getElementById('secondaryBeam').rows[5].cells[1].innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(SBppp,3) + '</span>';
  var TE = Math.round(SBppp/PBppp * 10000) / 100.0;
  document.getElementById('convertor').rows[1].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + TE + "%" + '</span>';
}

BeamCalc();

document.getElementById('ChargeStateP').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13 ){
      BeamCalc();      
    }
  }, false
);
document.getElementById('DCCTm').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      BeamCalc();      
    }
  }, false
);
document.getElementById('CSRmL').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      BeamCalc();      
    }
  }, false
);
document.getElementById('PBv').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      BeamCalc();      
    }
  }, false
);
document.getElementById('ChargeStateS').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      BeamCalc();      
    }
  }, false
);
document.getElementById('DCCTe').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      BeamCalc();      
    }
  }, false
);
document.getElementById('CSReL').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      BeamCalc();      
    }
  }, false
);
document.getElementById('SBv').addEventListener('keypress', 
  function(e){
    if(e.keyCode == 13){
      BeamCalc();      
    }
  }, false
);

