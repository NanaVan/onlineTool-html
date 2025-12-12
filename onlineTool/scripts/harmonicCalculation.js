function updateHarmonicTable(tabID, data, cols) {
    let harmTableContainer = getElementForTab(tabID, 'harmResult'); // 假设表格容器的 base ID 是 'harmResult'
    if (!harmTableContainer) {
        console.error(`Table container 'harmResult_${tabID}' not found for tab ${tabID}.`);
        return;
    }

    let newTBody = document.createElement('tbody');
    // 如果没有数据，或者数据为空数组，创建一个空行
    if (!data || data.length === 0) {
        let row = document.createElement('tr');
        for (let j = 0; j < cols; j++) {
            let cell = document.createElement('td');
            cell.appendChild(document.createTextNode("")); // 空格单元格
            row.appendChild(cell);
        }
        newTBody.appendChild(row);
    } else {
        // 遍历数据并创建表格行
        data.forEach(rowData => {
            let row = document.createElement('tr');
            for (let j = 0; j < cols; j++) {
                let cell = document.createElement('td');
                let cellText = document.createTextNode(rowData[j] || ""); // 防止数据为 undefined/null
                cell.appendChild(cellText);
                row.appendChild(cell);
            }
            newTBody.appendChild(row);
        });
    }

    // 查找并替换旧的 tbody 元素
    let oldTBody = harmTableContainer.querySelector('tbody');
    if (oldTBody) {
        harmTableContainer.replaceChild(newTBody, oldTBody);
    } else {
        // 如果容器中还没有 tbody，就直接添加
        harmTableContainer.appendChild(newTBody);
    }
}

/**
 * 辅助函数：更新革命频率和革命时间的结果显示。
 * @param {string} tabID 当前活动的 Tab ID。
 * @param {number|null} revFreqValue 革命频率值，如果为 null 则清空显示。
 * @param {number|null} revTimeValue 革命时间值，如果为 null 则清空显示。
 * @param {Array<Array<string>>|null} tableData 表格数据，如果为 null 则清空表格。
 */
function updateHarmonicResultsDisplay(tabID, revFreqValue, revTimeValue, tableData) {
    let resultRevFreqElement = getElementForTab(tabID, 'result_rev_freq');
    let resultRevTimeElement = getElementForTab(tabID, 'result_rev_time');

    if (resultRevFreqElement) {
        resultRevFreqElement.innerHTML = revFreqValue !== null ?
            `<span style="color:darkorange; font-family:Roboto;">${revFreqValue.toFixed(5)}</span>` : '';
    }
    if (resultRevTimeElement) {
        resultRevTimeElement.innerHTML = revTimeValue !== null ?
            `<span style="color:darkorange; font-family:Roboto;">${revTimeValue.toFixed(6)}</span>` : '';
    }
    // 同时也更新表格
    updateHarmonicTable(tabID, tableData, 3); // 假设表格始终是 3 列
}


/**
 * 执行谐波计算并更新相应 Tab 的 UI。
 * @param {string} tabID 当前活动的 Tab ID。
 */
function harmCalc(tabID) {
    const tabState = getTabState(tabID); // 获取当前 Tab 的独立状态

    // 1. 检查必要的前置数据 (velocity)
    if (typeof tabState.velocity === 'undefined' || tabState.velocity === null || tabState.velocity <= 0) {
        console.warn(`Velocity not set or invalid for tab ${tabID}. Cannot perform harmonic calculation.`);
        updateHarmonicResultsDisplay(tabID, null, null, []); // 清空结果
        return;
    }

    // 2. 获取 Tab 专属的输入元素值
    let lenOrbitInput = getElementForTab(tabID, 'length_orbit');
    let cenFreqInput = getElementForTab(tabID, 'center_frequency');
    let spanInput = getElementForTab(tabID, 'span');

    if (!lenOrbitInput || !cenFreqInput || !spanInput) {
        console.error(`Missing one or more input elements for harmonic calculation in tab: ${tabID}`);
        updateHarmonicResultsDisplay(tabID, null, null, []); // 清空结果
        return;
    }

    let lenOrbit = parseFloat(lenOrbitInput.value); // [m]
    let cenFreq = parseFloat(cenFreqInput.value); // [MHz]
    let span = parseFloat(spanInput.value); // [kHz]

    // 3. 验证输入值
    if (isNaN(lenOrbit) || isNaN(cenFreq) || isNaN(span) || lenOrbit === 0) {
        console.warn(`Invalid input for harmonic calculation in tab ${tabID}. Check length_orbit, center_frequency, span.`);
        updateHarmonicResultsDisplay(tabID, null, null, []); // 清空结果
        return;
    }

    // 4. 执行核心计算
    // 革命频率 [kHz]
    var frequency = tabState.velocity / lenOrbit * 1e4; // tabState.velocity is assumed to be in cm/s

    // 计算谐波范围
    var harm_start = Math.ceil((cenFreq * 1e3 - span / 2) / frequency);
    var harm_until = Math.floor((cenFreq * 1e3 + span / 2) / frequency);

    // 5. 准备表格数据
    const cols = 3; // 固定为 3 列
    var harmTableBodyData = [];

    if (harm_until >= harm_start) {
        for (var i = 0; i < (harm_until - harm_start + 1); i++) {
            var harmonicNumber = harm_start + i;
            var frequencyDifference = (harmonicNumber * frequency - cenFreq * 1e3).toFixed(6);
            var actualFrequency = (harmonicNumber * frequency * 1e-3).toFixed(9); // 转换为 MHz

            harmTableBodyData.push([
                harmonicNumber.toString(),
                frequencyDifference,
                actualFrequency
            ]);
        }
    } else {
        // 如果没有有效的谐波范围，添加一个空行
        harmTableBodyData.push(["", "", ""]);
    }

    // 6. 更新 UI
    updateHarmonicResultsDisplay(tabID, frequency, 1e6 / frequency, harmTableBodyData);
}

const harmonicInputBaseIds = ['length_orbit', 'center_frequency', 'span'];

all_tab_ids.forEach(tabID => {
    harmonicInputBaseIds.forEach(baseID => {
        const inputElement = getElementForTab(tabID, baseID);
        if (inputElement) {
            // 监听 'keypress' 事件，特别是回车键
            inputElement.addEventListener('keypress', function(e) {
                if (e.keyCode === 13) { // 检查是否是回车键
                    harmCalc(tabID); // 调用 harmCalc 并传入当前 Tab ID
                }
            });
            // 也可以监听 'change' 事件，当输入框失去焦点或值改变时触发
            inputElement.addEventListener('change', function() {
                harmCalc(tabID); // 调用 harmCalc 并传入当前 Tab ID
            });
        }
    });
});


