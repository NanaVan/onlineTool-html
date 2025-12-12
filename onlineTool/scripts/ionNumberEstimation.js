async function beamCalc(tabId, case_undefault) {
    if (!all_tab_ids.includes(tabId)) {
        console.error("Invalid tabID provided to beamCalc:", tabId);
        console.trace();
        return;
    }

    // --- 每次执行时重新加载 SQL.js ---
    const [SQL, buf] = await Promise.all([
        initSqlJs({ locateFile: file => `/dist/${file}` }),
        fetch("/dist/ionic_data.sqlite").then(res => res.arrayBuffer())
    ]);
    const db = new SQL.Database(new Uint8Array(buf));
    let stmt = db.prepare("SELECT MASS FROM IONICDATA WHERE Z=$zval AND Q=$qval AND A=$aval AND ISOMERIC=0");

    // collecting Primary Beam data
    const primaryBeamElementInput = getElementForTab(tabId, 'primaryBeam_element');
    const primaryBeamQInput = getElementForTab(tabId, 'primaryBeam_Q');
    const primaryBeamAInput = getElementForTab(tabId, 'primaryBeam_A');
    const DCCTmInput = getElementForTab(tabId, 'DCCTm');
    const CSRmLInput = getElementForTab(tabId, 'CSRmL');
    const primaryBeamBrhoInput = getElementForTab(tabId, 'primaryBeam_Brho');
    const primaryBeamEnergyInput = getElementForTab(tabId, 'primaryBeam_energy');
    const primaryBeamTable = getElementForTab(tabId, 'primaryBeam'); // table elements

    // collecting Secondary Beam data
    const secondaryBeamElementInput = getElementForTab(tabId, 'secondaryBeam_element');
    const secondaryBeamQInput = getElementForTab(tabId, 'secondaryBeam_Q');
    const secondaryBeamAInput = getElementForTab(tabId, 'secondaryBeam_A');
    const DCCTeInput = getElementForTab(tabId, 'DCCTe');
    const CSReLInput = getElementForTab(tabId, 'CSReL');
    const secondaryBeamBrhoInput = getElementForTab(tabId, 'secondaryBeam_Brho');
    const secondaryBeamEnergyInput = getElementForTab(tabId, 'secondaryBeam_energy');
    const secondaryBeamTable = getElementForTab(tabId, 'secondaryBeam'); // table elements

    // 获取 Convertor 表格
    const convertorTable = getElementForTab(tabId, 'convertor');

    // 检查所有必要的元素是否都已找到
    if (!primaryBeamElementInput || !primaryBeamQInput || !primaryBeamAInput ||
        !DCCTmInput || !CSRmLInput || !primaryBeamBrhoInput || !primaryBeamEnergyInput || !primaryBeamTable ||
        !secondaryBeamElementInput || !secondaryBeamQInput || !secondaryBeamAInput ||
        !DCCTeInput || !CSReLInput || !secondaryBeamBrhoInput || !secondaryBeamEnergyInput || !secondaryBeamTable ||
        !convertorTable) {
        console.error("Missing one or more required elements for beamCalc in tab:", tabId);
        db.close();
        return;
    }


    let primaryBeam_result = stmt.get({
        '$zval': elements.indexOf(primaryBeamElementInput.value) + 1,
        '$qval': primaryBeamQInput.value,
        '$aval': primaryBeamAInput.value
    });
    let secondaryBeam_result = stmt.get({
        '$zval': elements.indexOf(secondaryBeamElementInput.value) + 1,
        '$qval': secondaryBeamQInput.value,
        '$aval': secondaryBeamAInput.value
    });
    stmt.free();
    db.close();

    let PBppp = 0;
    let SBppp = 0;

    // === Primary Beam Calculation ===
    if (!primaryBeam_result || primaryBeam_result.length === 0) {
        if (primaryBeamTable && primaryBeamTable.rows[7] && primaryBeamTable.rows[7].cells[1]) { // 第7行是结果行
            primaryBeamTable.rows[7].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + 'No such ion!' + '</span>';
        }
    } else {
        let primaryBeam_Q = parseFloat(primaryBeamQInput.value);
        let DCCTm = parseFloat(DCCTmInput.value); // [μA]
        let LengthCSRm = parseFloat(CSRmLInput.value); // [m]
        let primaryBeam_mass = parseFloat(primaryBeam_result[0]); // [MeV]

        let primaryBeam_beta;
        let primaryBeam_gamma;

        if (case_undefault) { // 用户输入 Brho
            let primaryBeam_Brho = parseFloat(primaryBeamBrhoInput.value); // [Tm]
            let primaryBeam_gamma_beta = primaryBeam_Brho / primaryBeam_mass * primaryBeam_Q / speed_c / u2kg * elementary_charge;
            primaryBeam_beta = primaryBeam_gamma_beta / Math.sqrt(1 + Math.pow(primaryBeam_gamma_beta, 2));
            primaryBeam_gamma = 1 / Math.sqrt(1 - Math.pow(primaryBeam_beta, 2));
            let primaryBeam_energy = (primaryBeam_gamma - 1) / MeV2u;
            primaryBeamEnergyInput.value = primaryBeam_energy.toFixed(5);
        } else { // 用户输入 energy
            let primaryBeam_energy = parseFloat(primaryBeamEnergyInput.value); // [MeV/u]
            primaryBeam_gamma = 1 + primaryBeam_energy * MeV2u;
            primaryBeam_beta = Math.sqrt(1 - 1 / Math.pow(primaryBeam_gamma, 2));
            let primaryBeam_Brho = primaryBeam_gamma * primaryBeam_beta * primaryBeam_mass / primaryBeam_Q * speed_c * u2kg / elementary_charge; // [Tm]
            primaryBeamBrhoInput.value = primaryBeam_Brho.toFixed(5);
        }
        PBppp = DCCTm * 1e-6 / primaryBeam_Q / (primaryBeam_beta * speed_c) * LengthCSRm / elementary_charge;
        if (primaryBeamTable && primaryBeamTable.rows[7] && primaryBeamTable.rows[7].cells[1]) {
            primaryBeamTable.rows[7].cells[1].innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(PBppp, 3) + '</span>';
        }
    }

    // === Secondary Beam Calculation ===
    if (!secondaryBeam_result || secondaryBeam_result.length === 0) {
        if (secondaryBeamTable && secondaryBeamTable.rows[7] && secondaryBeamTable.rows[7].cells[1]) {
            secondaryBeamTable.rows[7].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + 'No such ion!' + '</span>';
        }
    } else {
        let secondaryBeam_Q = parseFloat(secondaryBeamQInput.value);
        let DCCTe = parseFloat(DCCTeInput.value); // [μA]
        let LengthCSRe = parseFloat(CSReLInput.value); // [m]
        let secondaryBeam_mass = parseFloat(secondaryBeam_result[0]); // [MeV]

        let secondaryBeam_beta;
        let secondaryBeam_gamma;

        if (case_undefault) { // 用户输入 energy
            let secondaryBeam_energy = parseFloat(secondaryBeamEnergyInput.value); // [MeV/u]
            secondaryBeam_gamma = 1 + secondaryBeam_energy * MeV2u;
            secondaryBeam_beta = Math.sqrt(1 - 1 / Math.pow(secondaryBeam_gamma, 2));
            let secondaryBeam_Brho = secondaryBeam_gamma * secondaryBeam_beta * secondaryBeam_mass / secondaryBeam_Q * speed_c * u2kg / elementary_charge; // [Tm]
            secondaryBeamBrhoInput.value = secondaryBeam_Brho.toFixed(5);
        } else { // 用户输入 Brho
            let secondaryBeam_Brho = parseFloat(secondaryBeamBrhoInput.value); // [Tm]
            let secondaryBeam_gamma_beta = secondaryBeam_Brho / secondaryBeam_mass * secondaryBeam_Q / speed_c / u2kg * elementary_charge;
            secondaryBeam_beta = secondaryBeam_gamma_beta / Math.sqrt(1 + Math.pow(secondaryBeam_gamma_beta, 2));
            secondaryBeam_gamma = 1 / Math.sqrt(1 - Math.pow(secondaryBeam_beta, 2));
            let secondaryBeam_energy = (secondaryBeam_gamma - 1) / MeV2u;
            secondaryBeamEnergyInput.value = secondaryBeam_energy.toFixed(5);
        }
        SBppp = DCCTe * 1e-6 / secondaryBeam_Q / (secondaryBeam_beta * speed_c) * LengthCSRe / elementary_charge;
        if (secondaryBeamTable && secondaryBeamTable.rows[7] && secondaryBeamTable.rows[0].cells[1]) {
            secondaryBeamTable.rows[7].cells[1].innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(SBppp, 3) + '</span>';
        }
    }

    // === Convertor Table Update ===
    if (primaryBeam_result && secondaryBeam_result && PBppp > 0) {
        var TE = Math.round(SBppp / PBppp * 10000) / 100.0;
        if (convertorTable && convertorTable.rows[1] && convertorTable.rows[1].cells[1]) {
            convertorTable.rows[1].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + TE + "%" + '</span>';
        }
    } else {
        if (convertorTable && convertorTable.rows[1] && convertorTable.rows[1].cells[1]) {
            convertorTable.rows[1].cells[1].innerHTML = '<span style="color:red; font-family:Roboto;">' + "Error" + '</span>';
        }
    }
}


/**
 * 限制和处理 Primary Beam 元素名称输入。
 * @param {HTMLInputElement} input 输入元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function limiter_primaryBeam_element(input, tabId) {
    var index = elements.indexOf(input.value);
    if (index == -1 || index < 3) {
        index = elements.indexOf('Fe');
        if (index === -1) {
            index = 0;
        }
        input.value = elements[index];
    }
    beamCalc(tabId, false); // 触发 beamCalc 重新计算
}

/**
 * 限制和处理 Secondary Beam 元素名称输入。
 * @param {HTMLInputElement} input 输入元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function limiter_secondaryBeam_element(input, tabId) {
    var index = elements.indexOf(input.value);
    if (index == -1 || index < 3) {
        index = elements.indexOf('Fe');
        if (index === -1) {
            index = 0;
        }
        input.value = elements[index];
    }
    beamCalc(tabId, false); // 触发 beamCalc 重新计算
}


// --- DOMContentLoaded 事件监听器 ---
document.addEventListener('DOMContentLoaded', async function() {
    // --- 定义一个辅助函数来绑定事件，减少重复代码 ---
    function bindInputEvents(baseId, handlerFunction, eventType = 'change') {
        all_tab_ids.forEach(tabId => {
            const inputElement = getElementForTab(tabId, baseId);
            if (inputElement) {
                inputElement.addEventListener(eventType, function() {
                    handlerFunction(inputElement, tabId); // 传递 inputElement 和 tabId
                });
                if (eventType === 'change' || eventType === 'input') {
                    inputElement.addEventListener('keypress', function(e) {
                        if (e.keyCode === 13) {
                            handlerFunction(inputElement, tabId); // 传递 inputElement 和 tabId
                        }
                    });
                }
            }
        });
    }

    // --- 绑定 beamCalc 相关的输入框事件 ---
    const beamCalcTriggerIds = [
        'primaryBeam_A', 'primaryBeam_Q', 'DCCTm', 'CSRmL', 'primaryBeam_energy', 'primaryBeam_Brho',
        'secondaryBeam_A', 'secondaryBeam_Q', 'DCCTe', 'CSReL', 'secondaryBeam_energy', 'secondaryBeam_Brho'
    ];

    all_tab_ids.forEach(tabId => {
        beamCalcTriggerIds.forEach(baseId => {
            const inputElement = getElementForTab(tabId, baseId);
            if (inputElement) {
                inputElement.addEventListener('keypress', function(e) {
                    if (e.keyCode === 13) {
                        let case_undefault = false;
                        if (baseId === 'primaryBeam_Brho' || baseId === 'secondaryBeam_energy') {
                            case_undefault = true;
                        }
                        beamCalc(tabId, case_undefault);
                    }
                }, false);
                inputElement.addEventListener('change', function() {
                    let case_undefault = false;
                    if (baseId === 'primaryBeam_Brho' || baseId === 'secondaryBeam_energy') {
                        case_undefault = true;
                    }
                    beamCalc(tabId, case_undefault);
                });
            }
        });
    });

    // 绑定元素名称的 limiter 函数
    bindInputEvents('primaryBeam_element', limiter_primaryBeam_element);
    bindInputEvents('secondaryBeam_element', limiter_secondaryBeam_element);

    // --- 初始调用 beamCalc ---
    // openRing 内部也会调用 beamCalc，但为了确保所有 Tab 在页面加载时都有一次初始计算，这里可以再次触发。
    // 如果 openRing 的调用顺序能保证 beamCalc 在 ionInformationSetting 之后，则这里的循环可以移除。
    // 为了稳妥，我们可以让它独立执行，因为它是独立的。
    for (const tabId of all_tab_ids) {
        try {
            await beamCalc(tabId, false);
        } catch (error) {
            console.error(`Error during initial beamCalc for tab ${tabId}:`, error);
        }
    }
});

