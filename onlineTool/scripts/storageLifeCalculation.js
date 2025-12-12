function calc_U_e(tabId) {
    const tabState = getTabState(tabId); // 获取当前 Tab 的状态

    // 检查所需的数据是否在 tabState 中
    if (typeof tabState.gamma === 'undefined' || tabState.gamma === null) {
        console.warn(`Invalid in tabState for tab ${tabId}. Cannot calculate U_e.`);
        let U_e_output = getElementForTab(tabId, 'U_e');
        if (U_e_output) U_e_output.innerHTML = '';
        return;
    }

    // me_kg 假设是电子质量的常数
    // elementary_charge 假设是基本电荷的常数
    // speed_c 假设是光速的常数

    // 计算 U_e
    var U_e = (tabState.gamma - 1) * me_kg * Math.pow(speed_c, 2) / elementary_charge * 1e-3; // [kV]

    // 更新 UI 元素
    let U_e_output = getElementForTab(tabId, 'U_e');
    if (U_e_output) {
        U_e_output.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + U_e.toFixed(3) + '</span>';
    } else {
        console.warn(`Element with ID 'U_e_${tabId}' not found for tab ${tabId}. Cannot display U_e.`);
    }
}


function LifeTimeREC(tabId) {
	/*
		estimate 1/e lifetime of the ion from Radioactive Electron Capture (REC)
		REC is the main mechanism of the beam loss under the ECooling

		REC@ H.Poth. Electron cooling: Theory, experiment, application. Physics Reports. 196 (1990): 135-297
    */

	const lifetimeREC_check = getElementForTab(tabId, 'lifetimeREC_check');
    if (!lifetimeREC_check || !lifetimeREC_check.checked) {
        // 如果复选框不存在或未被选中，则清空结果并退出
        getElementForTab(tabId, 'result_lifetime_rec').innerHTML = '';
        return;
    }
    const tabState = getTabState(tabId);

    // 从 Tab 状态中获取依赖的物理量
    const Q = parseFloat(getElementForTab(tabId, 'ion_charge').value); // ion_charge 也在 tabState 之外
    const e_gamma = tabState.gamma;
    const beta = tabState.beta; // 需要 beta，如果 tabState 中没有，需要从 gamma 推导或从 DOM 获取

    // 从 DOM 中获取其他输入
    const T_v = parseFloat(getElementForTab(tabId, 'T_v').value);
    const I_e = parseFloat(getElementForTab(tabId, 'I_e').value);
    const r_0 = parseFloat(getElementForTab(tabId, 'r_0').value);
    const len_C = parseFloat(getElementForTab(tabId, 'len_C').value);
    const len_cool = parseFloat(getElementForTab(tabId, 'len_cool').value);

    // 验证输入
    if (isNaN(Q) || isNaN(e_gamma) || isNaN(beta) || isNaN(T_v) || isNaN(I_e) || isNaN(r_0) || isNaN(len_C) || isNaN(len_cool) ||
        Q === 0 || e_gamma <= 1 || beta <= 0 || T_v === 0 || I_e === 0 || r_0 === 0 || len_C === 0 || len_cool === 0) {
        console.warn(`Invalid input for LifeTimeREC in tab ${tabId}.`);
        let result = getElementForTab(tabId, 'result_lifetime_rec');
        if (result) result.innerHTML = '';
        return;
    }

    // REC 计算
    var alpha_rec = 3.02e-13 * Math.pow(Q, 2) / T_v * (Math.log(11.32 * Q / Math.pow(T_v, 1 / 2)) + 0.14 * Math.pow((T_v / Math.pow(Q, 2)), 1 / 3)); // [cm^3/s]
    var rho = I_e / Math.PI / Math.pow(r_0, 2) / elementary_charge / speed_c / Math.pow(1 - 1 / Math.pow(e_gamma, 2), 1 / 2) * 1e-6; // [cm^-3]
    if (isNaN(rho) || rho === 0) { // 防止除以零或 rho 无效
        console.warn(`Calculated rho for REC is invalid in tab ${tabId}.`);
        let result = getElementForTab(tabId, 'result_lifetime_rec');
        if (result) result.innerHTML = '';
        return;
    }

    var tau_rec = Math.pow(e_gamma, 2) * len_C / len_cool / alpha_rec / rho; // [s]

    // 更新 UI
    let result = getElementForTab(tabId, 'result_lifetime_rec');
    if (result) {
        result.innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_rec, 3) + '</span>';
    }
}

function LifeTimeES(tabId) {
	/*
		estimate 1/e lifetime of the ion from Elastic Scattering (ES)
	   	ES is one of the main mechanism of the beam loss with the reaction with the residual gas in the storage ring

	   	with assumption: round vacuum chamber, beta_x = beta_z = <beta_y>, varepsilon_x = varepsilon_z = varepsilon
        single-scattering@ H.Poth. Electron cooling: Theory, experiment, application. Physics Reports. 196 (1990): 135-297
        multiple-scattering@ B.Franzke. Interaction of stored ion beams with the residual gas. CAS-CERN Accelerator School: 4th advanced accelerator physics course. (1992): 100-119
	*/
	const lifetimeREC_check = getElementForTab(tabId, 'lifetimeScattering_check');
    if (!lifetimeREC_check || !lifetimeREC_check.checked) {
        // 如果复选框不存在或未被选中，则清空结果并退出
        getElementForTab(tabId, 'result_lifetime_ss').innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_ms').innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_es').innerHTML = '';
        return;
    }
    const tabState = getTabState(tabId);

    // 从 Tab 状态中获取依赖的物理量
    const beta = tabState.beta;
    const gamma = tabState.gamma;
    const ion_mass_u = tabState.ion_mass; // 离子质量，单位 u
    const ion_Z_val = parseFloat(getElementForTab(tabId, 'ion_Z').value);

    // 常量
    const k_b = 1.380649e-23; // [J/K] Boltzmann's constant
    const r_e = 2.8179403262e-15; // [m] classical electron radius
    const u = 1.66053906660e-27; // [kg] atomic mass unit 

    // 从 DOM 中获取其他输入
    const gas_pressure_input = getElementForTab(tabId, 'gas_pressure');
    const gas_temperture_input = getElementForTab(tabId, 'gas_temperture');
    const selectMode = getElementForTab(tabId, 'selectMode');
    const emittance_input = getElementForTab(tabId, 'emittance');

    if (!gas_pressure_input || !gas_temperture_input || !selectMode || !emittance_input ||
        isNaN(beta) || beta <= 0 || isNaN(gamma) || gamma <= 1 || isNaN(ion_mass_u) || ion_mass_u <= 0 || isNaN(ion_Z_val) || ion_Z_val <= 0) {
        console.warn(`Missing or invalid core parameters for LifeTimeES in tab ${tabId}.`);
        getElementForTab(tabId, 'result_lifetime_ss').innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_ms').innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_es').innerHTML = '';
        return;
    }

    const gas_pressure = parseFloat(gas_pressure_input.value); // [mbar]
    const gas_temperture = parseFloat(gas_temperture_input.value) + 273.15; // [K]
    const gas_rho = gas_pressure * 1e2 / k_b / gas_temperture; // [m^-3]

    // 假设 ion_mass 变量是 MeV/c^2，如果 ion_mass_u 是 MeV/c^2，则需要转换成 u
    // 这里根据原始代码，`window.ion_mass / u` 意味着 `window.ion_mass` 单位是 MeV，所以需要除以 MeV2u 得到 amu (u)，再乘以 u 得到 kg
    // 但在 ionInformationSetting 中 ion_mass 存的是 MeV/c^2, 所以需要除以 MeV2u 得到 u
    var r_i = me_kg * r_e / (ion_mass_u / MeV2u * u); // [m] 修正离子质量单位

    // 获取模式相关参数
    var mode_index = selectMode.selectedIndex;
    var A_nu = 0;
    var average_beta_val = 0;
	if (tabId === 'SRing'){
    	switch (mode_index) {
    	    case 0:
    	        A_nu = parseFloat(getElementForTab(tabId, 'interTG_befDec_A_nu').value); // [πmm mrad]
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'interTG_befDec_beta').value); // [m]
    	        break;
    	    case 1:
    	        A_nu = parseFloat(getElementForTab(tabId, 'interTG_AfDec_A_nu').value); // [πmm mrad]
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'interTG_AfDec_beta').value); // [m]
    	        break;
    	    case 2:
    	        A_nu = parseFloat(getElementForTab(tabId, 'isochronous_1.43_A_nu').value);
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'isochronous_1.43_beta').value);
    	        break;
    	    case 3:
    	        A_nu = parseFloat(getElementForTab(tabId, 'isochronous_1.67_A_nu').value);
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'isochronous_1.67_beta').value);
    	        break;
    	    case 4:
    	        A_nu = parseFloat(getElementForTab(tabId, 'normal_large_A_nu').value);
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'normal_large_beta').value);
    	        break;
    	    case 5:
    	        A_nu = parseFloat(getElementForTab(tabId, 'normal_small_A_nu').value);
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'normal_small_beta').value);
    	        break;
    	    default:
    	        console.warn(`Invalid mode selected for LifeTimeES in tab ${tabId}.`);
    	        getElementForTab(tabId, 'result_lifetime_ss').innerHTML = '';
    	        getElementForTab(tabId, 'result_lifetime_ms').innerHTML = '';
    	        getElementForTab(tabId, 'result_lifetime_es').innerHTML = '';
    	        return;
		}
    }else if(tabId === 'CSRe'){
		switch (mode_index) {
    	    case 0:
    	        A_nu = parseFloat(getElementForTab(tabId, 'interTG_A_nu').value); // [πmm mrad]
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'interTG_beta').value); // [m]
    	        break;
    	    case 1:
    	        A_nu = parseFloat(getElementForTab(tabId, 'isochronous_A_nu').value);
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'isochronous_beta').value);
    	        break;
			case 2:
    	        A_nu = parseFloat(getElementForTab(tabId, 'normal_A_nu').value);
    	        average_beta_val = parseFloat(getElementForTab(tabId, 'normal_beta').value);
    	        break;
			default:
    	        console.warn(`Invalid mode selected for LifeTimeES in tab ${tabId}.`);
    	        getElementForTab(tabId, 'result_lifetime_ss').innerHTML = '';
    	        getElementForTab(tabId, 'result_lifetime_ms').innerHTML = '';
    	        getElementForTab(tabId, 'result_lifetime_es').innerHTML = '';
    	        return;
		}
	}else{
			console.warn(`Invalid mode selected for LifeTimeES in tab ${tabId}.`);
    	    getElementForTab(tabId, 'result_lifetime_ss').innerHTML = '';
    	    getElementForTab(tabId, 'result_lifetime_ms').innerHTML = '';
    	    getElementForTab(tabId, 'result_lifetime_es').innerHTML = '';
    	    return;
	}

    const emittance = parseFloat(emittance_input.value); // [m]
    var sqr_angle_acceptance = A_nu * 1e-6 / average_beta_val * Math.pow(1 - emittance / A_nu, 2); // only value, no unit

    var lambda_ss_part = 0;
    var lambda_ms_part = 0;
    for (var i = 0; i < 6; i++) {
        let Z_t = parseFloat(getElementForTab(tabId, `Gas${i}_Z_t`).value);
        let quantity = parseFloat(getElementForTab(tabId, `Gas${i}_quantity_t`).value);
        if (isNaN(Z_t) || isNaN(quantity)) continue;
        var lambda_ss_temp = Math.pow(Z_t, 2) * quantity;
        var lambda_ms_temp = Math.pow(Z_t, 2) * quantity * Math.log(204 / Math.pow(Z_t, 1 / 3));
        lambda_ss_part += lambda_ss_temp;
        lambda_ms_part += lambda_ms_temp;
    }

    if (isNaN(sqr_angle_acceptance) || sqr_angle_acceptance <= 0 || isNaN(gas_rho) || gas_rho <= 0 ||
        isNaN(lambda_ss_part) || lambda_ss_part === 0 || isNaN(lambda_ms_part) || lambda_ms_part === 0 ||
        isNaN(r_i) || r_i === 0) {
        console.warn(`Calculated intermediate parameters for LifeTimeES are invalid in tab ${tabId}.`);
        getElementForTab(tabId, 'result_lifetime_ss').innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_ms').innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_es').innerHTML = '';
        return;
    }

    var tau_ss = Math.pow(beta, 3) * Math.pow(gamma, 2) * sqr_angle_acceptance / 4 / Math.PI / Math.pow(ion_Z_val, 2) / Math.pow(r_i, 2) / speed_c / gas_rho / lambda_ss_part;
    var tau_ms = Math.pow(beta, 3) * Math.pow(gamma, 2) * (A_nu - emittance) * 1e-6 / average_beta_val / 32 / Math.pow(ion_Z_val, 2) / Math.pow(r_i, 2) / speed_c / gas_rho / lambda_ms_part;
    var tau_es = 1 / (1 / tau_ss + 1 / tau_ms);

    // 更新 UI
    getElementForTab(tabId, 'result_lifetime_ss').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_ss, 3) + '</span>';
    getElementForTab(tabId, 'result_lifetime_ms').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_ms, 3) + '</span>';
    getElementForTab(tabId, 'result_lifetime_es').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_es, 3) + '</span>';
}


function LifeTimeEC(tabId) {
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
    const tabState = getTabState(tabId);

    // 从 Tab 状态中获取依赖的物理量
    const beta = tabState.beta;
    const gamma = tabState.gamma;

    // 常量
    const k_b = 1.380649e-23; // [J/K] Boltzmann's constant

    // 从 DOM 中获取其他输入
    const gas_pressure_input = getElementForTab(tabId, 'gas_pressure');
    const gas_temperture_input = getElementForTab(tabId, 'gas_temperture');
    const ion_Z_input = getElementForTab(tabId, 'ion_Z');
    const ion_charge_input = getElementForTab(tabId, 'ion_charge');
    const ion_energy_MeVu_input = getElementForTab(tabId, 'ion_energy_MeVu');
    const ec_method_output = getElementForTab(tabId, 'EC_method');

    if (!gas_pressure_input || !gas_temperture_input || !ion_Z_input || !ion_charge_input || !ion_energy_MeVu_input ||
        isNaN(beta) || beta <= 0 || isNaN(gamma)) {
        console.warn(`Missing or invalid core parameters for LifeTimeEC in tab ${tabId}.`);
        if (ec_method_output) ec_method_output.innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_ec').innerHTML = '';
        return;
    }

    const gas_pressure = parseFloat(gas_pressure_input.value); // [mbar]
    const gas_temperture = parseFloat(gas_temperture_input.value) + 273.15; // [K]
    const gas_rho = gas_pressure * 1e2 / k_b / gas_temperture; // [m^-3]

    const ion_Z = parseFloat(ion_Z_input.value);
    const ion_charge = parseFloat(ion_charge_input.value);
    const ion_energy_MeVu = parseFloat(ion_energy_MeVu_input.value);

    var lambda_ec = 0;
    var current_method = "";

    // 检查并设置方法
    if ((gamma - 1) / MeV2u * 1e3 / Math.pow(ion_charge, 0.7) / Math.pow(18, 1.25) > 10 &&
        (gamma - 1) / MeV2u * 1e3 / Math.pow(ion_charge, 0.7) < 1000 && ion_charge >= 3) {
        current_method = "Schlacher Method";
        for (var i = 0; i < 6; i++) {
            let Z_t = parseFloat(getElementForTab(tabId, `Gas${i}_Z_t`).value);
            let quantity = parseFloat(getElementForTab(tabId, `Gas${i}_quantity_t`).value);
            if (isNaN(Z_t) || isNaN(quantity)) continue;

            var rho_partial = gas_rho * quantity; // [m^-3]
            var E_temp = (gamma - 1) / MeV2u * 1e3 / Math.pow(ion_Z, 0.7) / Math.pow(Z_t, 1.25); // [keV]
            var lambda_ec_temp = 1.1e-8 / Math.pow(E_temp, 4.8) * (1 - Math.exp(-0.037 * Math.pow(E_temp, 2.2))) * (1 - Math.exp(-2.44e-5 * Math.pow(E_temp, 2.6))) * 1e-4 * rho_partial * beta * speed_c * Math.pow(ion_charge, 0.5) / Math.pow(Z_t, 1.8); // [s^-1]
            lambda_ec += lambda_ec_temp;
        }
    } else if (ion_Z >= 36 && gamma <= 1.1) {
        current_method = "Franzke Method";
        var q_bar = ion_Z * (1 - Math.exp(-137 * beta / Math.pow(ion_Z, 0.67)));
        var b = (ion_charge < q_bar) ? 4 : 2;

        for (var i = 0; i < 6; i++) {
            let Z_t = parseFloat(getElementForTab(tabId, `Gas${i}_Z_t`).value);
            let quantity = parseFloat(getElementForTab(tabId, `Gas${i}_quantity_t`).value);
            if (isNaN(Z_t) || isNaN(quantity)) continue;

            var rho_partial = gas_rho * quantity; // [m^-3]
            var lambda_ec_temp = 2e-24 * 1e-4 * Math.pow(ion_Z, 1 / 2) * Math.pow(q_bar, 2) * Z_t * (1 - Math.exp(-137 * beta / Math.pow(Z_t, 0.67))) * Math.pow(Math.pow(gamma, 2) - 1, -2) * Math.pow(ion_charge / q_bar, b) * rho_partial * beta * speed_c; // [s^-1]
            lambda_ec += lambda_ec_temp;
        }
    } else {
        var default_state = true;
        default_state *= (ion_energy_MeVu >= 40 && ion_energy_MeVu <= 1000);
        var K = (ion_charge <= 18) ? 1 : (1.20 - 0.01 * ion_Z);

        for (var i = 0; i < 6; i++) {
            let Z_t = parseFloat(getElementForTab(tabId, `Gas${i}_Z_t`).value);
            let quantity = parseFloat(getElementForTab(tabId, `Gas${i}_quantity_t`).value);
            if (isNaN(Z_t) || isNaN(quantity)) continue;

            var rho_partial = gas_rho * quantity;
            var lambda_ec_temp = K * 1e-33 * Z_t * Math.pow(ion_Z, 5) * Math.pow(beta, -9) * 1e-4 * rho_partial * speed_c;
            lambda_ec += lambda_ec_temp;
            default_state *= (Z_t < 18);
        }
        current_method = default_state ? "Dmitriev Method" : "(Dmitriev Method)";
    }

    if (ec_method_output) {
        ec_method_output.innerHTML = current_method;
    }

    var tau_ec = (lambda_ec === 0 || isNaN(lambda_ec)) ? Infinity : 1 / lambda_ec;
    getElementForTab(tabId, 'result_lifetime_ec').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_ec, 3) + '</span>';
}


function LifeTimeEL(tabId) {
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

    const tabState = getTabState(tabId);

    // 从 Tab 状态中获取依赖的物理量
    const beta = tabState.beta;
    const gamma = tabState.gamma;

    // 常量
    const k_b = 1.380649e-23; // [J/K] Boltzmann's constant
    const alpha_0 = 5.29177210903e-11; // [m], Bohr radius
    const alpha = 7.2973525693e-3; // fine-structure constant

    // 从 DOM 中获取其他输入
    const gas_pressure_input = getElementForTab(tabId, 'gas_pressure');
    const gas_temperture_input = getElementForTab(tabId, 'gas_temperture');
    const ion_element_input = getElementForTab(tabId, 'ion_element');
    const ion_charge_input = getElementForTab(tabId, 'ion_charge');
    const ion_Z_input = getElementForTab(tabId, 'ion_Z');
    const el_method_output = getElementForTab(tabId, 'EL_method');
    const lifetimeElectronStrippingCheck = getElementForTab(tabId, 'lifetimeElectronStripping_check');

    if (!gas_pressure_input || !gas_temperture_input || !ion_element_input || !ion_charge_input || !ion_Z_input ||
        isNaN(beta) || beta <= 0 || isNaN(gamma)) {
        console.warn(`Missing or invalid core parameters for LifeTimeEL in tab ${tabId}.`);
        if (el_method_output) el_method_output.innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_el').innerHTML = '';
        return;
    }

    // 只有当 lifetimeElectronStrippingCheck 启用时才进行计算
    if (!lifetimeElectronStrippingCheck.checked) {
        if (el_method_output) el_method_output.innerHTML = '';
        getElementForTab(tabId, 'result_lifetime_el').innerHTML = '';
        return;
    }

    const gas_pressure = parseFloat(gas_pressure_input.value); // [mbar]
    const gas_temperture = parseFloat(gas_temperture_input.value) + 273.15; // [K]
    const gas_rho = gas_pressure * 1e2 / k_b / gas_temperture; // [m^-3]

    const ion_element = ion_element_input.value;
    const ion_charge = parseFloat(ion_charge_input.value);
    const ion_Z = parseFloat(ion_Z_input.value);
    const ion_num = ion_Z - ion_charge; // 束流中的电子数

    var lambda_el = 0;
    var current_method = "";

    // 检查并设置方法
    if (ion_Z >= 36 && gamma <= 1.1) {
        current_method = "Franzke Method";
        var q_bar = ion_Z * (1 - Math.exp(-137 * beta / Math.pow(ion_Z, 0.67)));
        var b = (ion_charge < q_bar) ? -2.3 : -4;

        for (var i = 0; i < 6; i++) {
            let Z_t = parseFloat(getElementForTab(tabId, `Gas${i}_Z_t`).value);
            let quantity = parseFloat(getElementForTab(tabId, `Gas${i}_quantity_t`).value);
            if (isNaN(Z_t) || isNaN(quantity)) continue;

            var rho_partial = gas_rho * quantity; // [m^-3]
            var lambda_el_temp = 3.5 * Math.pow(10, -18 + Math.pow(0.71 * Math.log(ion_charge), 1.5)) * 1e-4 * Math.pow(q_bar, -2) * Z_t * (1 - Math.exp(-137 * beta / Math.pow(Z_t, 0.67))) * Math.pow(Math.pow(gamma, 2) - 1, -1 / 2) * Math.pow(ion_charge / q_bar, b) * rho_partial * beta * speed_c; // [s^-1]
            lambda_el += lambda_el_temp;
        }
    } else {
        current_method = "Habs Method";
        for (var i = 0; i < 6; i++) {
            let Z_t = parseFloat(getElementForTab(tabId, `Gas${i}_Z_t`).value);
            let quantity = parseFloat(getElementForTab(tabId, `Gas${i}_quantity_t`).value);
            if (isNaN(Z_t) || isNaN(quantity)) continue;

            var rho_partial = gas_rho * quantity;
            var lambda_el_temp = 0;

            if (ion_num > 0) { // 只有非裸离子才可能发生电子剥离
                if (ion_Z >= Z_t) {
                    var I_ratio = 0;
                    for (var j = 1; j < ion_num + 1; j++) {
                        // 确保 bindEnergy[ion_element][j] 存在且不为零
                        if (bindEnergy && bindEnergy['H'] && bindEnergy['H'][1] && bindEnergy[ion_element] && bindEnergy[ion_element][j]) {
                            I_ratio += bindEnergy['H'][1] / bindEnergy[ion_element][j];
                        } else {
                            console.warn(`Missing bindEnergy data for ${ion_element} electron shell ${j} in tab ${tabId}.`);
                        }
                    }
                    lambda_el_temp = 4 * Math.PI * Math.pow(alpha_0, 2) * Math.pow(alpha / beta, 2) * (Z_t + 1) * Z_t * I_ratio * rho_partial * beta * speed_c; // [s^-1]
                } else if (ion_Z <= Math.pow(Z_t, 1 / 3)) {
                    var I_ratio = 0;
                    for (var j = 1; j < ion_num + 1; j++) {
                        if (bindEnergy && bindEnergy['H'] && bindEnergy['H'][1] && bindEnergy[ion_element] && bindEnergy[ion_element][j]) {
                            I_ratio += Math.pow(bindEnergy['H'][1] / bindEnergy[ion_element][j], 1 / 2);
                        } else {
                            console.warn(`Missing bindEnergy data for ${ion_element} electron shell ${j} in tab ${tabId}.`);
                        }
                    }
                    lambda_el_temp = Math.PI * Math.pow(alpha_0, 2) * alpha / beta * Math.pow(Z_t, 2 / 3) * I_ratio * rho_partial * beta * speed_c;
                }
            }
            lambda_el += lambda_el_temp;
        }
    }

    if (el_method_output) {
        el_method_output.innerHTML = current_method;
    }

    var tau_el = (lambda_el === 0 || isNaN(lambda_el)) ? Infinity : 1 / lambda_el;
    getElementForTab(tabId, 'result_lifetime_el').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_el, 3) + '</span>';
}



function LifeTimeNR(tabId) {
	/*
		estimate 1/e lifetime of the ion from Nuclear Reaction (NR)
		
		NR @ X.Y. Zhang. Study of Lifetime of Light Ion Beam stored in CSRm for Internal Target Experiment. (2005) MSc thesis 
	*/

    const tabState = getTabState(tabId);

    // 从 Tab 状态中获取依赖的物理量
    const beta = tabState.beta;
    const ion_mass_u = tabState.ion_mass; // 离子质量 MeV/c^2

    // 常量
    const k_b = 1.380649e-23; // [J/K] Boltzmann's constant

    // 从 DOM 中获取其他输入
    const gas_pressure_input = getElementForTab(tabId, 'gas_pressure');
    const gas_temperture_input = getElementForTab(tabId, 'gas_temperture');

    if (!gas_pressure_input || !gas_temperture_input ||
        isNaN(beta) || beta <= 0 || isNaN(ion_mass_u) || ion_mass_u <= 0) {
        console.warn(`Missing or invalid core parameters for LifeTimeNR in tab ${tabId}.`);
        getElementForTab(tabId, 'result_lifetime_nr').innerHTML = '';
        return;
    }

    const gas_pressure = parseFloat(gas_pressure_input.value); // [mbar]
    const gas_temperture = parseFloat(gas_temperture_input.value) + 273.15; // [K]
    const gas_rho = gas_pressure * 1e2 / k_b / gas_temperture; // [m^-3]

    // 假设 ion_mass 是 MeV/c^2，所以除以 MeV2u 得到 u
    var ion_A = Math.round(ion_mass_u / MeV2u); // 转换为原子质量单位并四舍五入

    var lambda_nr = 0;
    for (var i = 0; i < 6; i++) {
        let Z_t = parseFloat(getElementForTab(tabId, `Gas${i}_Z_t`).value);
        let quantity = parseFloat(getElementForTab(tabId, `Gas${i}_quantity_t`).value);
        if (isNaN(Z_t) || isNaN(quantity)) continue;

        var rho_partial = gas_rho * quantity; // [m^-3]
        // 修正：原始公式中 Math.pow(Z_t * 2, 1/3) 似乎不常见。通常是 Math.pow(A_t, 1/3) 或 Math.pow(Z_t, 1/3)。
        // 暂时保持原始公式逻辑，但请验证其物理意义。
        var lambda_nr_temp = Math.PI * Math.pow(Math.pow(ion_A, 1 / 3) + Math.pow(Z_t * 2, 1 / 3), 2) * 1e-30 / beta * rho_partial * speed_c;
        lambda_nr += lambda_nr_temp;
    }

    var tau_nr = (lambda_nr === 0 || isNaN(lambda_nr)) ? Infinity : 1 / lambda_nr;
    getElementForTab(tabId, 'result_lifetime_nr').innerHTML = '<span style="color:darkorange; font-family:Roboto;">' + toSci(tau_nr, 3) + '</span>';
}


// controller
function LifeTimeTotal(tabId, status) {
    if (status) {
        var lambda_REC = 0;
        var lambda_ES = 0;
        var lambda_EC = 0;
        var lambda_EL = 0;
        var lambda_NR = 0;

        // 获取各项寿命的倒数 (衰变常数 lambda)
        // 注意：这里需要检查 innerHTML 是否为空，并且值是否有效
        function getLambda(elementId) {
            const el = getElementForTab(tabId, elementId);
            if (el && el.textContent) {
                const val = parseFloat(toSciInverse(el.textContent)); // 使用 toSciInverse 确保能解析科学计数法
                if (!isNaN(val) && val > 0) return 1 / val;
            }
            return 0;
        }

        if (getElementForTab(tabId, 'lifetimeREC_check').checked) {
            lambda_REC = getLambda('result_lifetime_rec');
        }
        if (getElementForTab(tabId, 'lifetimeScattering_check').checked) {
            lambda_ES = getLambda('result_lifetime_es');
        }
        if (getElementForTab(tabId, 'lifetimeElectronCapture_check').checked) {
            lambda_EC = getLambda('result_lifetime_ec');
        }
        if (getElementForTab(tabId, 'lifetimeElectronStripping_check').checked) {
            lambda_EL = getLambda('result_lifetime_el');
        }
        if (getElementForTab(tabId, 'lifetimeNuclearReaction_check').checked) {
            lambda_NR = getLambda('result_lifetime_nr');
        }

        var lambda_total = lambda_REC + lambda_ES + lambda_EC + lambda_EL + lambda_NR;
        var tau_total = (lambda_total === 0 || isNaN(lambda_total)) ? Infinity : 1 / lambda_total;

        // 更新总寿命显示
        getElementForTab(tabId, 'result_lifetime_total_s').innerHTML = '<span style="color:red; font-family:Roboto;">' + toSci(tau_total, 3) + '</span>';
        getElementForTab(tabId, 'result_lifetime_total_min').innerHTML = '<span style="color:red; font-family:Roboto;">' + toSci(tau_total / 60, 3) + '</span>';
        getElementForTab(tabId, 'result_lifetime_total_hour').innerHTML = '<span style="color:red; font-family:Roboto;">' + toSci(tau_total / 3600, 3) + '</span>';

        // 更新各项寿命占比
        getElementForTab(tabId, 'ratio_lifetime_REC').innerHTML = '<span style="color:red; font-family:Roboto;">' + ((lambda_total !== 0) ? (lambda_REC / lambda_total * 100).toFixed(3) : 0) + '</span>';
        getElementForTab(tabId, 'ratio_lifetime_ES').innerHTML = '<span style="color:red; font-family:Roboto;">' + ((lambda_total !== 0) ? (lambda_ES / lambda_total * 100).toFixed(3) : 0) + '</span>';
        getElementForTab(tabId, 'ratio_lifetime_EC').innerHTML = '<span style="color:red; font-family:Roboto;">' + ((lambda_total !== 0) ? (lambda_EC / lambda_total * 100).toFixed(3) : 0) + '</span>';
        getElementForTab(tabId, 'ratio_lifetime_EL').innerHTML = '<span style="color:red; font-family:Roboto;">' + ((lambda_total !== 0) ? (lambda_EL / lambda_total * 100).toFixed(3) : 0) + '</span>';
        getElementForTab(tabId, 'ratio_lifetime_NR').innerHTML = '<span style="color:red; font-family:Roboto;">' + ((lambda_total !== 0) ? (lambda_NR / lambda_total * 100).toFixed(3) : 0) + '</span>';

    } else {
        // 清空结果
        getElementForTab(tabId, 'result_lifetime_total_s').innerHTML = "";
        getElementForTab(tabId, 'result_lifetime_total_min').innerHTML = "";
        getElementForTab(tabId, 'result_lifetime_total_hour').innerHTML = "";
        getElementForTab(tabId, 'ratio_lifetime_REC').innerHTML = "";
        getElementForTab(tabId, 'ratio_lifetime_ES').innerHTML = "";
        getElementForTab(tabId, 'ratio_lifetime_EC').innerHTML = "";
        getElementForTab(tabId, 'ratio_lifetime_EL').innerHTML = "";
        getElementForTab(tabId, 'ratio_lifetime_NR').innerHTML = "";
    }
}

// 辅助函数：从 toSci() 格式的字符串中解析数字
function toSciInverse(sciString) {
    if (!sciString) return NaN;
    // 移除 HTML span 标签
    const cleanedString = sciString.replace(/<span[^>]*>|<\/span>/g, '').trim();
    return parseFloat(cleanedString);
}


/**
 * 根据复选框状态锁定/解锁输入并触发总寿命计算。
 * @param {string} tabId 当前活动的 Tab ID。
 */
//function InformationLock(tabId) {
//    const isAnyChecked = getElementForTab(tabId, 'lifetimeREC_check').checked ||
//                         getElementForTab(tabId, 'lifetimeScattering_check').checked ||
//                         getElementForTab(tabId, 'lifetimeElectronCapture_check').checked ||
//                         getElementForTab(tabId, 'lifetimeElectronStripping_check').checked ||
//                         getElementForTab(tabId, 'lifetimeNuclearReaction_check').checked;
//
//    SettingInformationStatus(tabId, isAnyChecked);
//    RingModeInformationStatus(tabId, isAnyChecked);
//    GasInformationStatus(tabId, isAnyChecked);
//    LifeTimeTotal(tabId, isAnyChecked);
//}
//
function InformationLock(tabId) {
    const lifetimeREC_check_el = getElementForTab(tabId, 'lifetimeREC_check');
    const lifetimeScattering_check_el = getElementForTab(tabId, 'lifetimeScattering_check');
    const lifetimeElectronCapture_check_el = getElementForTab(tabId, 'lifetimeElectronCapture_check');
    const lifetimeElectronStripping_check_el = getElementForTab(tabId, 'lifetimeElectronStripping_check');
    const lifetimeNuclearReaction_check_el = getElementForTab(tabId, 'lifetimeNuclearReaction_check');

    // 检查所有复选框是否存在，如果不存在则视为未选中，避免报错
    const isREC_checked = lifetimeREC_check_el ? lifetimeREC_check_el.checked : false;
    const isES_checked = lifetimeScattering_check_el ? lifetimeScattering_check_el.checked : false;
    const isEC_checked = lifetimeElectronCapture_check_el ? lifetimeElectronCapture_check_el.checked : false;
    const isEL_checked = lifetimeElectronStripping_check_el ? lifetimeElectronStripping_check_el.checked : false;
    const isNR_checked = lifetimeNuclearReaction_check_el ? lifetimeNuclearReaction_check_el.checked : false;

    // 只要有一个复选框被选中，就视为需要“锁定”相关输入
    const isAnyLifetimeChecked = isREC_checked || isES_checked || isEC_checked || isEL_checked || isNR_checked;

    // 根据 isAnyLifetimeChecked 的状态来控制各个区域的禁用状态
    // 如果 isAnyLifetimeChecked 为 true，则 status 为 true（禁用）；
    // 如果 isAnyLifetimeChecked 为 false，则 status 为 false（启用/解锁）。
    SettingInformationStatus(tabId, isAnyLifetimeChecked);
    RingModeInformationStatus(tabId, isAnyLifetimeChecked);
    GasInformationStatus(tabId, isAnyLifetimeChecked);

    // 无论是否锁定，都尝试计算总寿命（如果 allChecked 为 false，LifeTimeTotal 内部会清空显示）
    LifeTimeTotal(tabId, isAnyLifetimeChecked);
}


/**
 * 锁定/解锁设置信息区域的输入。
 * @param {string} tabId 当前活动的 Tab ID。
 * @param {boolean} status true 为锁定，false 为解锁。
 */
function SettingInformationStatus(tabId, status) {
    const elementsToLock = [
        'ion_element', 'ion_Z', 'ion_A', 'ion_charge', 'ion_energy_MeVu',
        'ion_energy_AMeV', 'ion_totalKineticEnergy', 'ion_Brho', 'ion_gamma',
        'ion_beta', 'ion_velocity'
    ];
    elementsToLock.forEach(baseId => {
        const el = getElementForTab(tabId, baseId);
        if (el) el.disabled = status;
    });
}

/**
 * 锁定/解锁环模式信息区域的输入。
 * @param {string} tabId 当前活动的 Tab ID。
 * @param {boolean} status true 为锁定，false 为解锁。
 */
function RingModeInformationStatus(tabId, status) {
    const elementsToLock = [
        'selectMode', // select 元素
        'emittance',
        'interTG_befDec_A_nu', 'interTG_befDec_beta',
        'interTG_AfDec_A_nu', 'interTG_AfDec_beta',
        'isochronous_1.43_A_nu', 'isochronous_1.43_beta',
        'isochronous_1.67_A_nu', 'isochronous_1.67_beta',
        'normal_large_A_nu', 'normal_large_beta',
        'normal_small_A_nu', 'normal_small_beta'
    ];
    elementsToLock.forEach(baseId => {
        const el = getElementForTab(tabId, baseId);
        if (el) el.disabled = status;
    });
}

/**
 * 锁定/解锁气体信息区域的输入。
 * @param {string} tabId 当前活动的 Tab ID。
 * @param {boolean} status true 为锁定，false 为解锁。
 */
function GasInformationStatus(tabId, status) {
    const elementsToLock = ['gas_pressure', 'gas_temperture'];
    elementsToLock.forEach(baseId => {
        const el = getElementForTab(tabId, baseId);
        if (el) el.disabled = status;
    });
    for (var i = 0; i < 6; i++) {
        const Z_t_el = getElementForTab(tabId, `Gas${i}_Z_t`);
        const quantity_t_el = getElementForTab(tabId, `Gas${i}_quantity_t`);
        if (Z_t_el) Z_t_el.disabled = status;
        if (quantity_t_el) quantity_t_el.disabled = status;
    }
}


/**
 * 处理 REC 复选框状态变化。
 * @param {HTMLInputElement} checkboxElement 触发事件的复选框元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function selectREC(checkboxElement, tabId) {
    const elementsToToggle = ['T_v', 'I_e', 'r_0', 'len_C', 'len_cool'];
    elementsToToggle.forEach(baseId => {
        const el = getElementForTab(tabId, baseId);
        if (el) el.disabled = checkboxElement.checked;
    });

    if (checkboxElement.checked) {
        LifeTimeREC(tabId);
    } else {
        getElementForTab(tabId, 'result_lifetime_rec').innerHTML = "";
    }
    InformationLock(tabId);
}

/**
 * 处理 ES 复选框状态变化。
 * @param {HTMLInputElement} checkboxElement 触发事件的复选框元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function selectES(checkboxElement, tabId) {
    const elementsToToggle = ['selectMode', 'emittance'];
    elementsToToggle.forEach(baseId => {
        const el = getElementForTab(tabId, baseId);
        if (el) el.disabled = checkboxElement.checked;
    });

    // 模式选择器下的输入框也需要根据主复选框来启用/禁用
    const modeSelect = getElementForTab(tabId, 'selectMode');
    if(modeSelect) modeSelect.disabled = checkboxElement.checked;

    const modeSpecificInputs = [
        'interTG_befDec_A_nu', 'interTG_befDec_beta',
        'interTG_AfDec_A_nu', 'interTG_AfDec_beta',
        'isochronous_1.43_A_nu', 'isochronous_1.43_beta',
        'isochronous_1.67_A_nu', 'isochronous_1.67_beta',
        'normal_large_A_nu', 'normal_large_beta',
        'normal_small_A_nu', 'normal_small_beta'
    ];
    modeSpecificInputs.forEach(baseId => {
        const el = getElementForTab(tabId, baseId);
        if (el) el.disabled = checkboxElement.checked;
    });


    if (checkboxElement.checked) {
        LifeTimeES(tabId);
    } else {
        getElementForTab(tabId, 'result_lifetime_ss').innerHTML = "";
        getElementForTab(tabId, 'result_lifetime_ms').innerHTML = "";
        getElementForTab(tabId, 'result_lifetime_es').innerHTML = "";
    }
    InformationLock(tabId);
}

/**
 * 处理 EC 复选框状态变化。
 * @param {HTMLInputElement} checkboxElement 触发事件的复选框元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function selectEC(checkboxElement, tabId) {
    if (checkboxElement.checked) {
        LifeTimeEC(tabId);
    } else {
        getElementForTab(tabId, 'result_lifetime_ec').innerHTML = "";
        getElementForTab(tabId, 'EC_method').innerHTML = "";
    }
    InformationLock(tabId);
}

/**
 * 处理 EL 复选框状态变化。
 * @param {HTMLInputElement} checkboxElement 触发事件的复选框元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function selectEL(checkboxElement, tabId) {
    if (checkboxElement.checked) {
        LifeTimeEL(tabId);
    } else {
        getElementForTab(tabId, 'result_lifetime_el').innerHTML = "";
        getElementForTab(tabId, 'EL_method').innerHTML = "";
    }
    InformationLock(tabId);
}

/**
 * 处理 NR 复选框状态变化。
 * @param {HTMLInputElement} checkboxElement 触发事件的复选框元素。
 * @param {string} tabId 当前活动的 Tab ID。
 */
function selectNR(checkboxElement, tabId) {
    if (checkboxElement.checked) {
        LifeTimeNR(tabId);
    } else {
        getElementForTab(tabId, 'result_lifetime_nr').innerHTML = "";
    }
    InformationLock(tabId);
}


// --- DOMContentLoaded 事件监听器 ---
document.addEventListener('DOMContentLoaded', async function() {
    // 辅助函数，用于绑定输入事件到指定处理函数
    function bindInputToCalc(baseId, handlerFunction) {
        all_tab_ids.forEach(tabId => {
            const inputElement = getElementForTab(tabId, baseId);
            if (inputElement) {
                inputElement.addEventListener('change', () => handlerFunction(tabId));
                inputElement.addEventListener('keypress', (e) => {
                    if (e.keyCode === 13) handlerFunction(tabId);
                });
            }
        });
    }

    // 辅助函数，用于绑定复选框事件到指定处理函数
    function bindCheckboxToCalc(baseId, handlerFunction) {
        all_tab_ids.forEach(tabId => {
            const checkboxElement = getElementForTab(tabId, baseId);
            if (checkboxElement) {
                checkboxElement.addEventListener('change', () => handlerFunction(checkboxElement, tabId));
            }
        });
    }

    // --- 绑定各项寿命计算相关的输入事件 ---
    // REC inputs
    bindInputToCalc('T_v', LifeTimeREC);
    bindInputToCalc('I_e', LifeTimeREC);
    bindInputToCalc('r_0', LifeTimeREC);
    bindInputToCalc('len_C', LifeTimeREC);
    bindInputToCalc('len_cool', LifeTimeREC);

    // ES inputs (以及气体参数)
    bindInputToCalc('gas_pressure', () => { LifeTimeES(tabId); LifeTimeEC(tabId); LifeTimeEL(tabId); LifeTimeNR(tabId); }); // 气体压力会影响所有散射和捕获机制
    bindInputToCalc('gas_temperture', () => { LifeTimeES(tabId); LifeTimeEC(tabId); LifeTimeEL(tabId); LifeTimeNR(tabId); });
    for (let i = 0; i < 6; i++) {
        bindInputToCalc(`Gas${i}_Z_t`, () => { LifeTimeES(tabId); LifeTimeEC(tabId); LifeTimeEL(tabId); LifeTimeNR(tabId); });
        bindInputToCalc(`Gas${i}_quantity_t`, () => { LifeTimeES(tabId); LifeTimeEC(tabId); LifeTimeEL(tabId); LifeTimeNR(tabId); });
    }
    bindInputToCalc('selectMode', LifeTimeES);
    bindInputToCalc('emittance', LifeTimeES);

    // EC inputs (没有额外输入，依赖 ionInformationSetting 的 Z, Q, energy)
    // EL inputs (没有额外输入，依赖 ionInformationSetting 的 Z, Q, element)
    // NR inputs (没有额外输入，依赖 ionInformationSetting 的 mass)


    // --- 绑定复选框事件 ---
    bindCheckboxToCalc('lifetimeREC_check', selectREC);
    bindCheckboxToCalc('lifetimeScattering_check', selectES);
    bindCheckboxToCalc('lifetimeElectronCapture_check', selectEC);
    bindCheckboxToCalc('lifetimeElectronStripping_check', selectEL);
    bindCheckboxToCalc('lifetimeNuclearReaction_check', selectNR);


    // --- 初始计算 ---
    // 页面加载时，为所有 Tab 调用一次主要的计算函数
    // 注意：这些函数可能依赖于 tabStates 中已有的 beta, gamma, ion_mass
    // openTab.js 中的 openRing 应该会触发 initialIon_info_dynamic，从而更新 tabStates。
    // 因此，这里的初始调用可以放在 openRing 之后，或者像现在这样，等待 tabState 有值再计算。

    all_tab_ids.forEach(tabId => {
        try {
            // 首先计算 U_e，它是一个独立的量
            calc_U_e(tabId);

            // 初始时根据 checkbox 状态来判断是否计算各项寿命
            // 如果 checkbox 默认是未勾选的，那么这里不会触发 LifeTime* 计算，而是 InformationLock 会清空显示
            InformationLock(tabId);

            // 如果默认有勾选，则在这里会触发相应的 LifeTime* 计算
            if (getElementForTab(tabId, 'lifetimeREC_check').checked) LifeTimeREC(tabId);
            if (getElementForTab(tabId, 'lifetimeScattering_check').checked) LifeTimeES(tabId);
            if (getElementForTab(tabId, 'lifetimeElectronCapture_check').checked) LifeTimeEC(tabId);
            if (getElementForTab(tabId, 'lifetimeElectronStripping_check').checked) LifeTimeEL(tabId);
            if (getElementForTab(tabId, 'lifetimeNuclearReaction_check').checked) LifeTimeNR(tabId);
            LifeTimeTotal(tabId, true); // 确保总寿命也计算一次
        } catch (error) {
            console.error(`Error during initial storage life calculations for tab ${tabId}:`, error);
        }
    });

    // 为 selectMode 元素添加 change 事件监听器
    // 当模式改变时，需要重新计算 LifeTimeES
    all_tab_ids.forEach(tabId => {
        const selectModeElement = getElementForTab(tabId, 'selectMode');
        if (selectModeElement) {
            selectModeElement.addEventListener('change', () => LifeTimeES(tabId));
        }
    });

});
