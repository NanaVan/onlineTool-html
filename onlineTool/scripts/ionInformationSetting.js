// ionInformationSetting.js

// 假设以下全局函数和常量在加载此文件时已可用：
// - getElementForTab(tabId, baseId) (来自 openTab.js)
// - getTabState(tabId) (来自 index.html 的内联 script)
// - harmCalc(tabId) (来自 harmonicCalculation.js)
// - calc_U_e(tabId) (来自 storageLifeCalculation.js)
// - 所有 constantValues.js 中的常量 (如 elements, speed_c, u2kg, elementary_charge, MeV2u)


/**
 * 连接数据库，查询离子数据，并更新相应 Tab 的离子静态信息。
 * 也会更新 tabStates[tabId].ion_mass。
 * @param {string} tabId 当前活动的 Tab ID。
 */
async function checkIon(tabId) {
    // --- 每次执行时重新加载 SQL.js ---
    const [SQL, buf] = await Promise.all([
        initSqlJs({ locateFile: file => `/dist/${file}` }),
        fetch("/dist/ionic_data.sqlite").then(res => res.arrayBuffer())
    ]);
    const db = new SQL.Database(new Uint8Array(buf));
    let stmt = db.prepare("SELECT MASS FROM IONICDATA WHERE Z=$zval AND Q=$qval AND A=$aval AND ISOMERIC=0");

    const tabState = getTabState(tabId); // 获取当前 Tab 的状态

    let ionZInput = getElementForTab(tabId, 'ion_Z');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionInfoStaticTable = getElementForTab(tabId, 'ion_info_static');

    if (!ionZInput || !ionAInput || !ionChargeInput || !ionInfoStaticTable) {
        console.error(`Missing ion input/display elements for tab: ${tabId}. Cannot perform checkIon.`);
        db.close();
        return;
    }

    const currentZ = parseInt(ionZInput.value);
    const currentA = parseInt(ionAInput.value);
    const currentQ = parseInt(ionChargeInput.value);

    if (isNaN(currentZ) || isNaN(currentA) || isNaN(currentQ)) {
        console.warn(`Invalid Z, A, or Q value for tab ${tabId}. Cannot query database.`);
        if (ionInfoStaticTable.rows[4] && ionInfoStaticTable.rows[4].cells[1]) ionInfoStaticTable.rows[4].cells[1].innerHTML = '';
        if (ionInfoStaticTable.rows[5] && ionInfoStaticTable.rows[5].cells[1]) ionInfoStaticTable.rows[5].cells[1].innerHTML = '';
        if (ionInfoStaticTable.rows[6] && ionInfoStaticTable.rows[6].cells[1]) ionInfoStaticTable.rows[6].cells[1].innerHTML = '';
        tabState.ion_mass = 0; // 重置 Tab 状态中的质量
        db.close();
        return;
    }

    // 检查 A 范围并修正
    let stmtA = db.prepare("SELECT min(A), max(A) FROM IONICDATA WHERE Z=$zval AND ISOMERIC=0");
    let resultA = stmtA.get({'$zval': currentZ});
    stmtA.free();

    let correctedA = currentA;
    if (resultA && resultA[0] && resultA[1]) {
        if (currentA < resultA[0] || currentA > resultA[1]){
            console.log(`Ion A value (${currentA}) out of range [${resultA[0]}, ${resultA[1]}] for Z=${currentZ} in tab ${tabId}. Correcting to ${resultA[0]}.`);
            ionAInput.value = resultA[0]; // 修正 A 值
            correctedA = resultA[0];
        }
    } else {
        console.warn(`Could not determine A range for Z=${currentZ} in tab ${tabId}.`);
    }

    // 查询离子详细信息
    let stmtMass = db.prepare("SELECT MASS, TYPE, HALFLIFE FROM IONICDATA WHERE Z=$zval AND A=$aval AND Q=$qval AND ISOMERIC=0");
    let resultMass = stmtMass.get({'$zval': currentZ, '$aval': correctedA, '$qval': currentQ});
    stmtMass.free();
    db.close();

    if (!resultMass) {
        console.warn(`No ion data found for Z=${currentZ}, A=${correctedA}, Q=${currentQ} in tab ${tabId}.`);
        if (ionInfoStaticTable.rows[4] && ionInfoStaticTable.rows[4].cells[1]) ionInfoStaticTable.rows[4].cells[1].innerHTML = '';
        if (ionInfoStaticTable.rows[5] && ionInfoStaticTable.rows[5].cells[1]) ionInfoStaticTable.rows[5].cells[1].innerHTML = '';
        if (ionInfoStaticTable.rows[6] && ionInfoStaticTable.rows[6].cells[1]) ionInfoStaticTable.rows[6].cells[1].innerHTML = '';
        tabState.ion_mass = 0; // 重置 Tab 状态中的质量
        return;
    }

    // 更新静态信息表格
    // 注意：HTML 示例中的表格索引可能需要调整，这里假设 mass 在第一行 (index 0)
    if (ionInfoStaticTable.rows[4] && ionInfoStaticTable.rows[4].cells[1]) {
        ionInfoStaticTable.rows[4].cells[1].innerHTML = Number(resultMass[0]).toFixed(5);
    }
    if (ionInfoStaticTable.rows[5] && ionInfoStaticTable.rows[5].cells[1]) {
        ionInfoStaticTable.rows[5].cells[1].innerHTML = resultMass[1];
    }
    if (ionInfoStaticTable.rows[6] && ionInfoStaticTable.rows[6].cells[1]) {
        ionInfoStaticTable.rows[6].cells[1].innerHTML = resultMass[2];
    }

    tabState.ion_mass = resultMass[0]; // 更新 Tab 状态中的离子质量
}


/**
 * 限制和处理离子元素名称输入。
 * @param {HTMLInputElement} input 输入元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function limiter_ion_element(input, tabId) {
    let ionZInput = getElementForTab(tabId, 'ion_Z');
    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let lifetimeElectronStrippingCheck = getElementForTab(tabId, 'lifetimeElectronStripping_check');

    if (!ionZInput || !ionChargeInput || !lifetimeElectronStrippingCheck) {
        console.error(`Missing limiter_ion_element related inputs for tab: ${tabId}`);
        return;
    }

    let index = elements.indexOf(input.value);
    if (index === -1) {
        let currentZ = parseInt(ionZInput.value);
        if (!isNaN(currentZ) && currentZ >= 1 && currentZ <= elements.length) {
            input.value = elements[currentZ - 1];
            ionZInput.value = currentZ;
        } else {
            input.value = 'Fe';
            ionZInput.value = elements.indexOf('Fe') + 1;
        }
    } else {
        ionZInput.value = String(index + 1);
    }
    ionChargeInput.value = ionZInput.value;
    lifetimeElectronStrippingCheck.disabled = true;
    initialIon_info_dynamic(tabId);
}


/**
 * 限制和处理离子 Z 值输入。
 * @param {HTMLInputElement} input 输入元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function limiter_ion_Z(input, tabId) {
    let ionElementInput = getElementForTab(tabId, 'ion_element');
    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let lifetimeElectronStrippingCheck = getElementForTab(tabId, 'lifetimeElectronStripping_check');

    if (!ionElementInput || !ionChargeInput || !lifetimeElectronStrippingCheck) {
        console.error(`Missing limiter_ion_Z related inputs for tab: ${tabId}`);
        return;
    }

    let val = parseInt(input.value);
    if (isNaN(val) || val < 1 || val > elements.length) {
        input.value = '26';
        val = 26;
    }
    ionElementInput.value = elements[val - 1];
    ionChargeInput.value = val;
    lifetimeElectronStrippingCheck.disabled = true;
    initialIon_info_dynamic(tabId);
}


/**
 * 限制和处理离子 A 值输入。
 * @param {HTMLInputElement} input 输入元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function limiter_ion_A(input, tabId) {
    let lifetimeElectronStrippingCheck = getElementForTab(tabId, 'lifetimeElectronStripping_check');
    if (lifetimeElectronStrippingCheck) {
        lifetimeElectronStrippingCheck.disabled = true;
    }
    initialIon_info_dynamic(tabId);
}


/**
 * 限制和处理离子电荷输入。
 * @param {HTMLInputElement} input 输入元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function limiter_ion_charge(input, tabId) {
    let ionZInput = getElementForTab(tabId, 'ion_Z');
    let lifetimeElectronStrippingCheck = getElementForTab(tabId, 'lifetimeElectronStripping_check');

    if (!ionZInput || !lifetimeElectronStrippingCheck) {
        console.error(`Missing limiter_ion_charge related inputs for tab: ${tabId}`);
        return;
    }

    const currentIonZ = parseInt(ionZInput.value);
    let currentIonCharge = parseInt(input.value);

    if (isNaN(currentIonZ) || isNaN(currentIonCharge)) {
        console.warn(`Invalid Z or Q value for tab ${tabId}.`);
        return;
    }

    if (currentIonCharge < 0) {
        currentIonCharge = 0;
        input.value = currentIonCharge;
    } else if (currentIonCharge > currentIonZ) {
        input.value = currentIonZ;
        currentIonCharge = currentIonZ;
    } else if (currentIonCharge < currentIonZ - 3) {
        input.value = currentIonZ - 3;
        currentIonCharge = currentIonZ - 3;
    }

    if (currentIonCharge < currentIonZ) {
        lifetimeElectronStrippingCheck.disabled = false;
    } else {
        lifetimeElectronStrippingCheck.disabled = true;
    }
    initialIon_info_dynamic(tabId);
}


/**
 * 初始化并计算离子动态信息（beta, gamma, velocity, energy等）。
 * 每次离子参数变化时都会调用。
 * 更新 tabStates[tabId] 中的 beta, gamma, velocity 等，并触发 harmCalc 和 calc_U_e。
 * @param {string} tabId 当前活动的 Tab ID。
 */
async function initialIon_info_dynamic(tabId) {
    await checkIon(tabId); // 首先更新离子质量

    const tabState = getTabState(tabId); // 获取当前 Tab 的状态

    // 检查 ion_mass 是否有效，如果无效则停止计算
    if (tabState.ion_mass === 0 || isNaN(tabState.ion_mass)) {
        console.warn(`ion_mass is 0 or invalid for tab ${tabId}. Cannot calculate dynamic ion info.`);
        // 清空所有动态信息显示
        clearDynamicIonInfo(tabId);
        // 也清空 tabState 中的动态值
        tabState.beta = 0;
        tabState.gamma = 0;
        tabState.velocity = 0;
        return;
    }

    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionBrhoInput = getElementForTab(tabId, 'ion_Brho');
    let ionEnergyMeVuInput = getElementForTab(tabId, 'ion_energy_MeVu');
    let ionEnergyAMeVInput = getElementForTab(tabId, 'ion_energy_AMeV');
    let ionTotalKineticEnergyInput = getElementForTab(tabId, 'ion_totalKineticEnergy');
    let ionGammaInput = getElementForTab(tabId, 'ion_gamma');
    let ionBetaInput = getElementForTab(tabId, 'ion_beta');
    let ionVelocityInput = getElementForTab(tabId, 'ion_velocity');

    if (!ionChargeInput || !ionAInput || !ionBrhoInput || !ionEnergyMeVuInput || !ionVelocityInput ||
        !ionEnergyAMeVInput || !ionTotalKineticEnergyInput || !ionGammaInput || !ionBetaInput) {
        console.error(`Missing required dynamic ion info input/display elements for tab: ${tabId}`);
        clearDynamicIonInfo(tabId);
        tabState.beta = 0;
        tabState.gamma = 0;
        tabState.velocity = 0;
        return;
    }

    var Q = parseFloat(ionChargeInput.value);
    var A = parseFloat(ionAInput.value);
    var Brho = 1; // 你的代码默认为 1 [Tm]

    if (isNaN(Q) || isNaN(A) || Q === 0 || A === 0) {
        console.warn(`Invalid Q or A for tab ${tabId}. Cannot calculate dynamic ion info.`);
        clearDynamicIonInfo(tabId);
        tabState.beta = 0;
        tabState.gamma = 0;
        tabState.velocity = 0;
        return;
    }

    // 物理量计算，并更新 Tab 状态
    var gamma_beta = Brho / tabState.ion_mass * Q / speed_c / u2kg * elementary_charge;
    tabState.beta = gamma_beta / Math.sqrt(1 + Math.pow(gamma_beta, 2));
    tabState.velocity = tabState.beta * speed_c * 1e-7; // [cm/s]
    tabState.gamma = 1 / Math.sqrt(1 - Math.pow(tabState.beta, 2));

    var energy_MeVu = (tabState.gamma - 1) / MeV2u;
    var TKE = tabState.ion_mass * energy_MeVu;
    var energy_AMeV = TKE / A;

    // 更新 UI 元素
    ionEnergyMeVuInput.value = energy_MeVu.toFixed(5);
    ionEnergyAMeVInput.value = energy_AMeV.toFixed(5);
    ionTotalKineticEnergyInput.value = TKE.toFixed(5);
    ionBrhoInput.value = Brho.toFixed(5);
    ionGammaInput.value = tabState.gamma.toFixed(5);
    ionBetaInput.value = tabState.beta.toFixed(5);
    ionVelocityInput.value = tabState.velocity.toFixed(5);

    // 触发其他模块的计算
    if (typeof harmCalc === 'function') {
        harmCalc(tabId);
    }
    if (typeof calc_U_e === 'function') {
        calc_U_e(tabId);
    }
}

/**
 * 清空动态离子信息显示。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function clearDynamicIonInfo(tabId) {
    const idsToClear = ['ion_energy_MeVu', 'ion_energy_AMeV', 'ion_totalKineticEnergy', 'ion_Brho', 'ion_gamma', 'ion_beta', 'ion_velocity'];
    idsToClear.forEach(baseId => {
        const el = getElementForTab(tabId, baseId);
        if (el) el.value = '';
    });
}


function change_ion_energy_MeVu(input, tabId) {
    const tabState = getTabState(tabId);
    if (tabState.ion_mass === 0 || isNaN(tabState.ion_mass)) { console.warn(`Ion mass is 0 or invalid for tab ${tabId}. Cannot calculate.`); clearDynamicIonInfo(tabId); return; }

    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionEnergyAMeVInput = getElementForTab(tabId, 'ion_energy_AMeV');
    let ionTotalKineticEnergyInput = getElementForTab(tabId, 'ion_totalKineticEnergy');
    let ionBrhoInput = getElementForTab(tabId, 'ion_Brho');
    let ionGammaInput = getElementForTab(tabId, 'ion_gamma');
    let ionBetaInput = getElementForTab(tabId, 'ion_beta');
    let ionVelocityInput = getElementForTab(tabId, 'ion_velocity');

    if (!ionChargeInput || !ionAInput || !ionEnergyAMeVInput || !ionTotalKineticEnergyInput ||
        !ionBrhoInput || !ionGammaInput || !ionBetaInput || !ionVelocityInput) {
        console.error(`Missing required display elements for change_ion_energy_MeVu in tab ${tabId}.`);
        return;
    }

    let Q = parseFloat(ionChargeInput.value);
    let A = parseFloat(ionAInput.value);
    let energy_MeVu = parseFloat(input.value);

    if (isNaN(Q) || isNaN(A) || isNaN(energy_MeVu) || Q === 0 || A === 0 || energy_MeVu < 0) {
        console.warn(`Invalid input for change_ion_energy_MeVu in tab ${tabId}.`);
        clearDynamicIonInfo(tabId);
        tabState.beta = 0;
        tabState.gamma = 0;
        tabState.velocity = 0;
        return;
    }

    let TKE = tabState.ion_mass * energy_MeVu;
    let energy_AMeV = TKE / A;

    tabState.gamma = 1 + energy_MeVu * MeV2u;
    tabState.beta = Math.sqrt(1 - 1 / Math.pow(tabState.gamma, 2));
    tabState.velocity = tabState.beta * speed_c * 1e-7;

    let Brho = tabState.gamma * tabState.beta * tabState.ion_mass / Q * speed_c * u2kg / elementary_charge;

    ionEnergyAMeVInput.value = energy_AMeV.toFixed(5);
    ionTotalKineticEnergyInput.value = TKE.toFixed(5);
    ionBrhoInput.value = Brho.toFixed(5);
    ionGammaInput.value = tabState.gamma.toFixed(5);
    ionBetaInput.value = tabState.beta.toFixed(5);
    ionVelocityInput.value = tabState.velocity.toFixed(5);

    if (typeof harmCalc === 'function') { harmCalc(tabId); }
    if (typeof calc_U_e === 'function') { calc_U_e(tabId); }
}

function change_ion_energy_AMeV(input, tabId) {
    const tabState = getTabState(tabId);
    if (tabState.ion_mass === 0 || isNaN(tabState.ion_mass)) { console.warn(`Ion mass is 0 or invalid for tab ${tabId}. Cannot calculate.`); clearDynamicIonInfo(tabId); return; }

    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionEnergyMeVuInput = getElementForTab(tabId, 'ion_energy_MeVu');
    let ionTotalKineticEnergyInput = getElementForTab(tabId, 'ion_totalKineticEnergy');
    let ionBrhoInput = getElementForTab(tabId, 'ion_Brho');
    let ionGammaInput = getElementForTab(tabId, 'ion_gamma');
    let ionBetaInput = getElementForTab(tabId, 'ion_beta');
    let ionVelocityInput = getElementForTab(tabId, 'ion_velocity');

    if (!ionChargeInput || !ionAInput || !ionEnergyMeVuInput || !ionTotalKineticEnergyInput ||
        !ionBrhoInput || !ionGammaInput || !ionBetaInput || !ionVelocityInput) {
        console.error(`Missing required display elements for change_ion_energy_AMeV in tab ${tabId}.`);
        return;
    }

    let Q = parseFloat(ionChargeInput.value);
    let A = parseFloat(ionAInput.value);
    let energy_AMeV = parseFloat(input.value);

    if (isNaN(Q) || isNaN(A) || isNaN(energy_AMeV) || Q === 0 || A === 0 || energy_AMeV < 0) {
        console.warn(`Invalid input for change_ion_energy_AMeV in tab ${tabId}.`);
        clearDynamicIonInfo(tabId);
        tabState.beta = 0;
        tabState.gamma = 0;
        tabState.velocity = 0;
        return;
    }

    let TKE = A * energy_AMeV;
    let energy_MeVu = TKE / tabState.ion_mass;

    tabState.gamma = 1 + energy_MeVu * MeV2u;
    tabState.beta = Math.sqrt(1 - 1 / Math.pow(tabState.gamma, 2));
    tabState.velocity = tabState.beta * speed_c * 1e-7;

    let Brho = tabState.gamma * tabState.beta * tabState.ion_mass / Q * speed_c * u2kg / elementary_charge;

    ionEnergyMeVuInput.value = energy_MeVu.toFixed(5);
    ionTotalKineticEnergyInput.value = TKE.toFixed(5);
    ionBrhoInput.value = Brho.toFixed(5);
    ionGammaInput.value = tabState.gamma.toFixed(5);
    ionBetaInput.value = tabState.beta.toFixed(5);
    ionVelocityInput.value = tabState.velocity.toFixed(5);

    if (typeof harmCalc === 'function') { harmCalc(tabId); }
    if (typeof calc_U_e === 'function') { calc_U_e(tabId); }
}

function change_ion_totalKineticEnergy(input, tabId) {
    const tabState = getTabState(tabId);
    if (tabState.ion_mass === 0 || isNaN(tabState.ion_mass)) { console.warn(`Ion mass is 0 or invalid for tab ${tabId}. Cannot calculate.`); clearDynamicIonInfo(tabId); return; }

    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionEnergyMeVuInput = getElementForTab(tabId, 'ion_energy_MeVu');
    let ionEnergyAMeVInput = getElementForTab(tabId, 'ion_energy_AMeV');
    let ionBrhoInput = getElementForTab(tabId, 'ion_Brho');
    let ionGammaInput = getElementForTab(tabId, 'ion_gamma');
    let ionBetaInput = getElementForTab(tabId, 'ion_beta');
    let ionVelocityInput = getElementForTab(tabId, 'ion_velocity');

    if (!ionChargeInput || !ionAInput || !ionEnergyMeVuInput || !ionEnergyAMeVInput ||
        !ionBrhoInput || !ionGammaInput || !ionBetaInput || !ionVelocityInput) {
        console.error(`Missing required display elements for change_ion_totalKineticEnergy in tab ${tabId}.`);
        return;
    }

    let Q = parseFloat(ionChargeInput.value);
    let A = parseFloat(ionAInput.value);
    let TKE = parseFloat(input.value);

    if (isNaN(Q) || isNaN(A) || isNaN(TKE) || Q === 0 || A === 0 || TKE < 0) {
        console.warn(`Invalid input for change_ion_totalKineticEnergy in tab ${tabId}.`);
        clearDynamicIonInfo(tabId);
        tabState.beta = 0;
        tabState.gamma = 0;
        tabState.velocity = 0;
        return;
    }

    let energy_MeVu = TKE / tabState.ion_mass;
    let energy_AMeV = TKE / A;

    tabState.gamma = 1 + energy_MeVu * MeV2u;
    tabState.beta = Math.sqrt(1 - 1 / Math.pow(tabState.gamma, 2));
    tabState.velocity = tabState.beta * speed_c * 1e-7;

    let Brho = tabState.gamma * tabState.beta * tabState.ion_mass / Q * speed_c * u2kg / elementary_charge;

    ionEnergyMeVuInput.value = energy_MeVu.toFixed(5);
    ionEnergyAMeVInput.value = energy_AMeV.toFixed(5);
    ionBrhoInput.value = Brho.toFixed(5);
    ionGammaInput.value = tabState.gamma.toFixed(5);
    ionBetaInput.value = tabState.beta.toFixed(5);
    ionVelocityInput.value = tabState.velocity.toFixed(5);

    if (typeof harmCalc === 'function') { harmCalc(tabId); }
    if (typeof calc_U_e === 'function') { calc_U_e(tabId); }
}

function change_ion_Brho(input, tabId) {
    const tabState = getTabState(tabId);
    if (tabState.ion_mass === 0 || isNaN(tabState.ion_mass)) { console.warn(`Ion mass is 0 or invalid for tab ${tabId}. Cannot calculate.`); clearDynamicIonInfo(tabId); return; }

    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionEnergyMeVuInput = getElementForTab(tabId, 'ion_energy_MeVu');
    let ionEnergyAMeVInput = getElementForTab(tabId, 'ion_energy_AMeV');
    let ionTotalKineticEnergyInput = getElementForTab(tabId, 'ion_totalKineticEnergy');
    let ionGammaInput = getElementForTab(tabId, 'ion_gamma');
    let ionBetaInput = getElementForTab(tabId, 'ion_beta');
    let ionVelocityInput = getElementForTab(tabId, 'ion_velocity');

    if (!ionChargeInput || !ionAInput || !ionEnergyMeVuInput || !ionEnergyAMeVInput ||
        !ionTotalKineticEnergyInput || !ionGammaInput || !ionBetaInput || !ionVelocityInput) {
        console.error(`Missing required display elements for change_ion_Brho in tab ${tabId}.`);
        return;
    }

    let Q = parseFloat(ionChargeInput.value);
    let A = parseFloat(ionAInput.value);
    let Brho = parseFloat(input.value);

    if (isNaN(Q) || isNaN(A) || isNaN(Brho) || Q === 0 || A === 0 || Brho < 0) {
        console.warn(`Invalid input for change_ion_Brho in tab ${tabId}.`);
        clearDynamicIonInfo(tabId);
        tabState.beta = 0;
        tabState.gamma = 0;
        tabState.velocity = 0;
        return;
    }

    let gamma_beta = Brho / tabState.ion_mass * Q / speed_c / u2kg * elementary_charge;
    tabState.beta = gamma_beta / Math.sqrt(1 + Math.pow(gamma_beta, 2));
    tabState.velocity = tabState.beta * speed_c * 1e-7;
    tabState.gamma = 1 / Math.sqrt(1 - Math.pow(tabState.beta, 2));

    let energy_MeVu = (tabState.gamma - 1) / MeV2u;
    let TKE = tabState.ion_mass * energy_MeVu;
    let energy_AMeV = TKE / A;

    ionEnergyMeVuInput.value = energy_MeVu.toFixed(5);
    ionEnergyAMeVInput.value = energy_AMeV.toFixed(5);
    ionTotalKineticEnergyInput.value = TKE.toFixed(5);
    ionGammaInput.value = tabState.gamma.toFixed(5);
    ionBetaInput.value = tabState.beta.toFixed(5);
    ionVelocityInput.value = tabState.velocity.toFixed(5);

    if (typeof harmCalc === 'function') { harmCalc(tabId); }
    if (typeof calc_U_e === 'function') { calc_U_e(tabId); }
}

function change_ion_gamma(input, tabId) {
    const tabState = getTabState(tabId);
    if (tabState.ion_mass === 0 || isNaN(tabState.ion_mass)) { console.warn(`Ion mass is 0 or invalid for tab ${tabId}. Cannot calculate.`); clearDynamicIonInfo(tabId); return; }

    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionEnergyMeVuInput = getElementForTab(tabId, 'ion_energy_MeVu');
    let ionEnergyAMeVInput = getElementForTab(tabId, 'ion_energy_AMeV');
    let ionTotalKineticEnergyInput = getElementForTab(tabId, 'ion_totalKineticEnergy');
    let ionBrhoInput = getElementForTab(tabId, 'ion_Brho');
    let ionBetaInput = getElementForTab(tabId, 'ion_beta');
    let ionVelocityInput = getElementForTab(tabId, 'ion_velocity');

    if (!ionChargeInput || !ionAInput || !ionEnergyMeVuInput || !ionEnergyAMeVInput ||
        !ionTotalKineticEnergyInput || !ionBrhoInput || !ionBetaInput || !ionVelocityInput) {
        console.error(`Missing required display elements for change_ion_gamma in tab ${tabId}.`);
        return;
    }

    let Q = parseFloat(ionChargeInput.value);
    let A = parseFloat(ionAInput.value);
    tabState.gamma = parseFloat(input.value);

    if (isNaN(Q) || isNaN(A) || isNaN(tabState.gamma) || Q === 0 || A === 0 || tabState.gamma < 1) {
        console.warn(`Invalid input for change_ion_gamma in tab ${tabId}.`);
        clearDynamicIonInfo(tabId);
        tabState.beta = 0;
        tabState.velocity = 0;
        return;
    }

    tabState.beta = Math.sqrt(1 - 1 / Math.pow(tabState.gamma, 2));
    tabState.velocity = tabState.beta * speed_c * 1e-7;
    let Brho = tabState.gamma * tabState.beta * tabState.ion_mass / Q * speed_c * u2kg / elementary_charge;

    let energy_MeVu = (tabState.gamma - 1) / MeV2u;
    let TKE = tabState.ion_mass * energy_MeVu;
    let energy_AMeV = TKE / A;

    ionEnergyMeVuInput.value = energy_MeVu.toFixed(5);
    ionEnergyAMeVInput.value = energy_AMeV.toFixed(5);
    ionTotalKineticEnergyInput.value = TKE.toFixed(5);
    ionBrhoInput.value = Brho.toFixed(5);
    ionBetaInput.value = tabState.beta.toFixed(5);
    ionVelocityInput.value = tabState.velocity.toFixed(5);

    if (typeof harmCalc === 'function') { harmCalc(tabId); }
    if (typeof calc_U_e === 'function') { calc_U_e(tabId); }
}

function change_ion_beta(input, tabId) {
    const tabState = getTabState(tabId);
    if (tabState.ion_mass === 0 || isNaN(tabState.ion_mass)) { console.warn(`Ion mass is 0 or invalid for tab ${tabId}. Cannot calculate.`); clearDynamicIonInfo(tabId); return; }

    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionEnergyMeVuInput = getElementForTab(tabId, 'ion_energy_MeVu');
    let ionEnergyAMeVInput = getElementForTab(tabId, 'ion_energy_AMeV');
    let ionTotalKineticEnergyInput = getElementForTab(tabId, 'ion_totalKineticEnergy');
    let ionBrhoInput = getElementForTab(tabId, 'ion_Brho');
    let ionGammaInput = getElementForTab(tabId, 'ion_gamma');
    let ionVelocityInput = getElementForTab(tabId, 'ion_velocity');

    if (!ionChargeInput || !ionAInput || !ionEnergyMeVuInput || !ionEnergyAMeVInput ||
        !ionTotalKineticEnergyInput || !ionBrhoInput || !ionGammaInput || !ionVelocityInput) {
        console.error(`Missing required display elements for change_ion_beta in tab ${tabId}.`);
        return;
    }

    let Q = parseFloat(ionChargeInput.value);
    let A = parseFloat(ionAInput.value);
    tabState.beta = parseFloat(input.value);

    if (isNaN(Q) || isNaN(A) || isNaN(tabState.beta) || Q === 0 || A === 0 || tabState.beta < 0 || tabState.beta >= 1) {
        console.warn(`Invalid input for change_ion_beta in tab ${tabId}.`);
        clearDynamicIonInfo(tabId);
        tabState.gamma = 0;
        tabState.velocity = 0;
        return;
    }

    tabState.gamma = 1 / Math.sqrt(1 - Math.pow(tabState.beta, 2));
    tabState.velocity = tabState.beta * speed_c * 1e-7;
    let Brho = tabState.gamma * tabState.beta * tabState.ion_mass / Q * speed_c * u2kg / elementary_charge;

    let energy_MeVu = (tabState.gamma - 1) / MeV2u;
    let TKE = tabState.ion_mass * energy_MeVu;
    let energy_AMeV = TKE / A;

    ionEnergyMeVuInput.value = energy_MeVu.toFixed(5);
    ionEnergyAMeVInput.value = energy_AMeV.toFixed(5);
    ionTotalKineticEnergyInput.value = TKE.toFixed(5);
    ionBrhoInput.value = Brho.toFixed(5);
    ionGammaInput.value = tabState.gamma.toFixed(5);
    ionVelocityInput.value = tabState.velocity.toFixed(5);

    if (typeof harmCalc === 'function') { harmCalc(tabId); }
    if (typeof calc_U_e === 'function') { calc_U_e(tabId); }
}

function change_ion_velocity(input, tabId) {
    const tabState = getTabState(tabId);
    if (tabState.ion_mass === 0 || isNaN(tabState.ion_mass)) { console.warn(`Ion mass is 0 or invalid for tab ${tabId}. Cannot calculate.`); clearDynamicIonInfo(tabId); return; }

    let ionChargeInput = getElementForTab(tabId, 'ion_charge');
    let ionAInput = getElementForTab(tabId, 'ion_A');
    let ionEnergyMeVuInput = getElementForTab(tabId, 'ion_energy_MeVu');
    let ionEnergyAMeVInput = getElementForTab(tabId, 'ion_energy_AMeV');
    let ionTotalKineticEnergyInput = getElementForTab(tabId, 'ion_totalKineticEnergy');
    let ionBrhoInput = getElementForTab(tabId, 'ion_Brho');
    let ionGammaInput = getElementForTab(tabId, 'ion_gamma');
    let ionBetaInput = getElementForTab(tabId, 'ion_beta');

    if (!ionChargeInput || !ionAInput || !ionEnergyMeVuInput || !ionEnergyAMeVInput ||
        !ionTotalKineticEnergyInput || !ionBrhoInput || !ionGammaInput || !ionBetaInput) {
        console.error(`Missing required display elements for change_ion_velocity in tab ${tabId}.`);
        return;
    }

    let Q = parseFloat(ionChargeInput.value);
    let A = parseFloat(ionAInput.value);
    tabState.velocity = parseFloat(input.value);

    if (isNaN(Q) || isNaN(A) || isNaN(tabState.velocity) || Q === 0 || A === 0 || tabState.velocity < 0) {
        console.warn(`Invalid input for change_ion_velocity in tab ${tabId}.`);
        clearDynamicIonInfo(tabId);
        tabState.beta = 0;
        tabState.gamma = 0;
        return;
    }

    tabState.beta = tabState.velocity * 1e7 / speed_c;
    if (tabState.beta >= 1) {
        console.warn(`Calculated beta >= 1 for tab ${tabId}. Velocity too high. Setting to 0.999999999.`);
        tabState.beta = 0.999999999;
    }

    tabState.gamma = 1 / Math.sqrt(1 - Math.pow(tabState.beta, 2));
    let Brho = tabState.gamma * tabState.beta * tabState.ion_mass / Q * speed_c * u2kg / elementary_charge;

    let energy_MeVu = (tabState.gamma - 1) / MeV2u;
    let TKE = tabState.ion_mass * energy_MeVu;
    let energy_AMeV = TKE / A;

    ionEnergyMeVuInput.value = energy_MeVu.toFixed(5);
    ionEnergyAMeVInput.value = energy_AMeV.toFixed(5);
    ionTotalKineticEnergyInput.value = TKE.toFixed(5);
    ionBrhoInput.value = Brho.toFixed(5);
    ionGammaInput.value = tabState.gamma.toFixed(5);
    ionBetaInput.value = tabState.beta.toFixed(5);
    // ion_velocity input is already updated by user

    if (typeof harmCalc === 'function') { harmCalc(tabId); }
    if (typeof calc_U_e === 'function') { calc_U_e(tabId); }
}

// --- DOMContentLoaded 事件监听器 ---
document.addEventListener('DOMContentLoaded', async function() {
    // --- 定义一个辅助函数来绑定事件，减少重复代码 ---
    function bindInputEvents(baseId, handlerFunction, eventType = 'change') {
        all_tab_ids.forEach(tabId => {
            const inputElement = getElementForTab(tabId, baseId);
            if (inputElement) {
                inputElement.addEventListener(eventType, function() {
                    handlerFunction(inputElement, tabId);
                });
                if (eventType === 'change' || eventType === 'input') {
                    inputElement.addEventListener('keypress', function(e) {
                        if (e.keyCode === 13) {
                            handlerFunction(inputElement, tabId);
                        }
                    });
                }
            }
        });
    }

    // --- 绑定 ionInformationSetting 相关的输入框事件 ---
    bindInputEvents('ion_element', limiter_ion_element);
    bindInputEvents('ion_Z', limiter_ion_Z);
    bindInputEvents('ion_A', limiter_ion_A);
    bindInputEvents('ion_charge', limiter_ion_charge);
    bindInputEvents('ion_energy_MeVu', change_ion_energy_MeVu);
    bindInputEvents('ion_energy_AMeV', change_ion_energy_AMeV);
    bindInputEvents('ion_totalKineticEnergy', change_ion_totalKineticEnergy);
    bindInputEvents('ion_Brho', change_ion_Brho);
    bindInputEvents('ion_gamma', change_ion_gamma);
    bindInputEvents('ion_beta', change_ion_beta);
    bindInputEvents('ion_velocity', change_ion_velocity);


    // --- 初始离子信息设置 ---
    // 这个函数也会触发 harmCalc 和 calc_U_e
    all_tab_ids.forEach(async tabId => {
        try {
            console.log(`Initializing ion dynamic information for tab: ${tabId}`);
            await initialIon_info_dynamic(tabId);
        } catch (error) {
            console.error(`Error during initial initialIon_info_dynamic for tab ${tabId}:`, error);
        }
    });

    // 禁用电子剥离寿命检查 (如果你的应用逻辑需要在初始化时禁用)
    all_tab_ids.forEach(tabId => {
        const lifetimeElectronStrippingCheck = getElementForTab(tabId, 'lifetimeElectronStripping_check');
        if (lifetimeElectronStrippingCheck) {
            lifetimeElectronStrippingCheck.disabled = true;
        } else {
            console.warn(`Element with ID 'lifetimeElectronStripping_check_${tabId}' not found. Skipping disable.`);
        }
    });
});
