// openTab.js
const all_tab_ids = ['SRing', 'CSRe']; // 确保这里和所有模块中的 all_tab_ids 一致

function openRing(evt, ringName) {
  var i, tabcontent, tablinks;

  tabcontent = document.getElementsByClassName("tabcontent");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  tablinks = document.getElementsByClassName("tablinks");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  document.getElementById(ringName).style.display = "block";

  if (evt && evt.currentTarget) {
    evt.currentTarget.className += " active";
  } else {
      const tabButton = document.querySelector(`.tablinks[onclick*="openRing(event, '${ringName}')"]`);
      if (tabButton) {
          tabButton.className += " active";
      }
  }

  // --- 协调各模块的计算触发 ---
  // ionInformationSetting 优先，因为它产生其他模块所需的核心离子数据
  if (typeof initialIon_info_dynamic === 'function') {
      // initialIon_info_dynamic 内部会触发 harmCalc 和 calc_U_e
      initialIon_info_dynamic(ringName).catch(e => console.error(`Error in initialIon_info_dynamic for tab ${ringName}:`, e));
  } else {
      console.warn("initialIon_info_dynamic function not found. Ion info will not update on tab switch.");
  }

  // ionNumberEstimation 是独立的，可以在任何时候触发
  if (typeof beamCalc === 'function') {
      beamCalc(ringName, false).catch(e => console.error(`Error in beamCalc for tab ${ringName}:`, e));
  } else {
      console.warn("beamCalc function not found. Beam estimation will not update on tab switch.");
  }
}

// helper function to get an element by its base ID and the current tab ID
function getElementForTab(tabID, baseID) {
  const fullID = `${baseID}_${tabID}`;
  const element = document.getElementById(fullID);
  return element;
}
