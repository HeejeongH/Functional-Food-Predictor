// API Base URL
const API_BASE_URL = '';

// 전역 변수
let collectedProteins = [];
let lastTrainedModel = null;

// 데이터 수집 (Step 1 사전 실행)
async function collectData() {
    const targetGenes = document.getElementById('target-genes').value.trim();
    
    if (!targetGenes) {
        alert('수집할 단백질 리스트를 입력해주세요.');
        return;
    }
    
    const btn = document.getElementById('collect-btn');
    btn.disabled = true;
    btn.innerHTML = '<i class="fas fa-spinner fa-spin"></i> 수집 중...';
    
    try {
        const targetList = targetGenes.split(',').map(s => s.trim()).filter(s => s);
        
        const response = await axios.post(`${API_BASE_URL}/api/data/collect`, {
            target_list: targetList,
            standard_type: 'IC50'
        });
        
        // 수집된 단백질 저장
        collectedProteins = targetList;
        
        // UI 업데이트
        document.getElementById('collected-proteins').classList.remove('hidden');
        document.getElementById('collected-proteins-list').innerHTML = 
            `수집된 단백질: <strong>${collectedProteins.join(', ')}</strong> (${response.data.total_compounds}개 화합물)`;
        
        // 단백질 선택 드롭다운 활성화
        const proteinSelect = document.getElementById('protein-name');
        proteinSelect.disabled = false;
        proteinSelect.innerHTML = collectedProteins.map(p => 
            `<option value="${p}">${p}</option>`
        ).join('');
        
        // Step 1 체크박스 자동 선택
        document.getElementById('step1').checked = true;
        updateStepCard(1);
        updateStepStatus('step1', 'success');
        
        alert('데이터 수집 완료! 이제 단백질을 선택하고 다음 단계를 진행하세요.');
        
    } catch (error) {
        alert('데이터 수집 실패: ' + (error.response?.data?.detail || error.message));
    } finally {
        btn.disabled = false;
        btn.innerHTML = '<i class="fas fa-download"></i> 데이터 수집';
    }
}

// 단계 선택 토글
function toggleStep(stepNum) {
    const checkbox = document.getElementById(`step${stepNum}`);
    checkbox.checked = !checkbox.checked;
    updateStepCard(stepNum);
    updateStepConfig(stepNum);
}

// 단계 직접 선택
function selectStep(stepNum) {
    const checkbox = document.getElementById(`step${stepNum}`);
    checkbox.checked = true;
    updateStepCard(stepNum);
    updateStepConfig(stepNum);
}

// 단계별 설정 폼 표시/숨김
function updateStepConfig(stepNum) {
    // Step 1: 데이터 수집 폼은 항상 표시
    const step1Config = document.getElementById('step1-config');
    if (step1Config) {
        step1Config.style.display = 'block'; // 항상 표시
    }
    
    // Step 5와 Step 6 체크 상태 확인
    const step5Checkbox = document.getElementById('step5');
    const step6Checkbox = document.getElementById('step6');
    const step5Checked = step5Checkbox ? step5Checkbox.checked : false;
    const step6Checked = step6Checkbox ? step6Checkbox.checked : false;
    
    // Step 5: SMILES 입력 폼 (Step 5만 선택되었을 때만 표시)
    const step5Config = document.getElementById('step5-config');
    if (step5Config) {
        step5Config.style.display = step5Checked ? 'block' : 'none';
    }
    
    // Step 6: FooDB 옵션 폼 (Step 6이 선택되었을 때만 표시)
    const step6ConfigMode = document.getElementById('step6-config-mode');
    const step6ConfigTop = document.getElementById('step6-config-top');
    const step6ConfigThreshold = document.getElementById('step6-config-threshold');
    
    if (step6ConfigMode) step6ConfigMode.style.display = step6Checked ? 'block' : 'none';
    if (step6ConfigTop) step6ConfigTop.style.display = step6Checked ? 'block' : 'none';
    if (step6ConfigThreshold) step6ConfigThreshold.style.display = step6Checked ? 'block' : 'none';
}

// 단계 카드 업데이트
function updateStepCard(stepNum) {
    const checkbox = document.getElementById(`step${stepNum}`);
    const card = checkbox.closest('.step-card');
    
    if (checkbox.checked) {
        card.classList.add('selected');
    } else {
        card.classList.remove('selected');
    }
}

// 워크플로우 실행
async function executeWorkflow() {
    // 선택된 단계 확인 (Step 1 제외 - 사전에 수집해야 함)
    const steps = [
        { id: 'step2', name: '특성 변환', fn: executeFeatureTransform },
        { id: 'step3', name: '모델 학습', fn: executeModelTraining },
        { id: 'step4', name: 'SHAP 분석', fn: executeSHAPAnalysis },
        { id: 'step5', name: 'SMILES 예측', fn: executeSMILESPrediction },
        { id: 'step6', name: 'FooDB 예측', fn: executeFooDBPrediction }
    ];

    const selectedSteps = steps.filter(step => document.getElementById(step.id).checked);

    if (selectedSteps.length === 0) {
        alert('실행할 단계를 최소 1개 이상 선택해주세요.');
        return;
    }

    // 유효성 검사
    const proteinName = document.getElementById('protein-name').value;
    if (!proteinName) {
        alert('먼저 데이터를 수집하고 단백질을 선택해주세요.');
        return;
    }

    // Step 5 선택 시 SMILES 필수
    if (selectedSteps.find(s => s.id === 'step5')) {
        const smiles = document.getElementById('smiles-input').value.trim();
        if (!smiles) {
            alert('SMILES 예측을 위해 SMILES 리스트를 입력해주세요.');
            return;
        }
    }

    // 진행 상황 섹션 표시
    document.getElementById('progress-section').classList.remove('hidden');
    document.getElementById('results-section').classList.add('hidden');
    
    const progressContainer = document.getElementById('progress-container');
    const resultsContainer = document.getElementById('results-container');
    progressContainer.innerHTML = '';
    resultsContainer.innerHTML = '';

    // 각 단계 실행
    for (let i = 0; i < selectedSteps.length; i++) {
        const step = selectedSteps[i];
        const stepNumber = steps.findIndex(s => s.id === step.id) + 1;

        // 진행 상황 아이템 생성
        const progressItem = createProgressItem(stepNumber, step.name, 'running');
        progressContainer.appendChild(progressItem);

        // 상태 업데이트
        updateStepStatus(step.id, 'running');

        try {
            // 단계 실행
            const result = await step.fn();
            
            // 성공 처리
            updateProgressItem(progressItem, 'completed', result.message || '완료');
            updateStepStatus(step.id, 'success');

            // 결과 카드 추가
            const resultCard = createResultCard(stepNumber, step.name, result);
            resultsContainer.appendChild(resultCard);
            document.getElementById('results-section').classList.remove('hidden');

        } catch (error) {
            // 실패 처리
            console.error(`${step.name} 오류:`, error);
            
            let errorMsg = error.message || '오류 발생';
            
            // error.response?.data?.detail이 객체인 경우 JSON으로 변환
            if (error.response?.data?.detail) {
                const detail = error.response.data.detail;
                if (typeof detail === 'object') {
                    errorMsg = JSON.stringify(detail, null, 2);
                } else {
                    errorMsg = detail;
                }
            }
            
            updateProgressItem(progressItem, 'failed', errorMsg);
            updateStepStatus(step.id, 'error');
            
            alert(`${step.name} 실패:\n${errorMsg}`);
            break; // 실패 시 중단
        }
    }
}

// 진행 상황 아이템 생성
function createProgressItem(stepNumber, stepName, status) {
    const item = document.createElement('div');
    item.className = `progress-item ${status}`;
    
    let iconClass = 'fa-circle';
    if (status === 'running') iconClass = 'fa-spinner fa-spin';
    if (status === 'completed') iconClass = 'fa-check-circle';
    if (status === 'failed') iconClass = 'fa-times-circle';
    
    item.innerHTML = `
        <div class="progress-icon ${status}">
            <i class="fas ${iconClass}"></i>
        </div>
        <div class="progress-content">
            <div class="progress-title">${stepNumber}. ${stepName}</div>
            <div class="progress-detail">${status === 'running' ? '실행 중...' : ''}</div>
        </div>
    `;
    
    return item;
}

// 진행 상황 업데이트
function updateProgressItem(item, status, message) {
    item.className = `progress-item ${status}`;
    
    const icon = item.querySelector('.progress-icon');
    icon.className = `progress-icon ${status}`;
    
    let iconClass = 'fa-circle';
    if (status === 'running') iconClass = 'fa-spinner fa-spin';
    if (status === 'completed') iconClass = 'fa-check-circle';
    if (status === 'failed') iconClass = 'fa-times-circle';
    
    icon.innerHTML = `<i class="fas ${iconClass}"></i>`;
    
    const detail = item.querySelector('.progress-detail');
    detail.textContent = message;
}

// 단계 상태 업데이트
function updateStepStatus(stepId, status) {
    const stepNum = stepId.replace('step', '');
    const statusElement = document.getElementById(`status${stepNum}`);
    statusElement.className = `step-status ${status}`;
    
    let iconClass = 'fa-circle';
    let text = '대기중';
    
    if (status === 'running') {
        iconClass = 'fa-spinner fa-spin';
        text = '실행중';
    } else if (status === 'success') {
        iconClass = 'fa-check-circle';
        text = '완료';
    } else if (status === 'error') {
        iconClass = 'fa-times-circle';
        text = '실패';
    }
    
    statusElement.innerHTML = `<i class="fas ${iconClass}"></i><span>${text}</span>`;
}

// 결과 카드 생성
function createResultCard(stepNumber, stepName, result) {
    const card = document.createElement('div');
    card.className = 'result-card fade-in';
    
    // Step별 특화 UI
    if (stepNumber === 6 || stepNumber === 5) {
        // FooDB 또는 SMILES 예측 결과 - 특별한 UI
        card.innerHTML = createPredictionResultUI(stepNumber, stepName, result);
    } else if (stepNumber === 4) {
        // SHAP 분석 결과 - 특별한 UI
        card.innerHTML = createSHAPResultUI(stepNumber, stepName, result);
    } else if (stepNumber === 3) {
        // 모델 학습 결과 - 특별한 UI
        card.innerHTML = createModelResultUI(stepNumber, stepName, result);
    } else {
        // 기본 UI
        card.innerHTML = createDefaultResultUI(stepNumber, stepName, result);
    }
    
    return card;
}

// 기본 결과 UI
function createDefaultResultUI(stepNumber, stepName, result) {
    let content = `<h4><i class="fas fa-check-circle" style="color: #10b981; margin-right: 0.5rem;"></i>${stepNumber}. ${stepName} 결과</h4>`;
    content += '<div class="result-grid">';
    
    for (const [key, value] of Object.entries(result)) {
        if (key !== 'message') {
            const label = formatLabel(key);
            content += `
                <div class="result-item">
                    <div class="result-label">${label}</div>
                    <div class="result-value">${formatValue(key, value)}</div>
                </div>
            `;
        }
    }
    
    content += '</div>';
    return content;
}

// 예측 결과 UI (FooDB, SMILES)
function createPredictionResultUI(stepNumber, stepName, result) {
    const totalCount = result.total_count || 0;
    const activeCount = result.active_count || 0;
    const inactiveCount = result.inactive_count || 0;
    const filteredCount = result.filtered_count || result.top_predictions?.length || 0;
    const threshold = result.threshold || 0.5;
    const predictions = result.top_predictions || result.predictions || [];
    
    const activePercent = totalCount > 0 ? ((activeCount / totalCount) * 100).toFixed(1) : 0;
    const inactivePercent = totalCount > 0 ? ((inactiveCount / totalCount) * 100).toFixed(1) : 0;
    
    return `
        <div style="padding: 1.5rem;">
            <!-- 헤더 -->
            <div style="display: flex; align-items: center; gap: 1rem; margin-bottom: 1.5rem;">
                <div style="width: 48px; height: 48px; background: linear-gradient(135deg, #6366f1 0%, #8b5cf6 100%); border-radius: 12px; display: flex; align-items: center; justify-content: center;">
                    <i class="fas fa-flask" style="color: white; font-size: 24px;"></i>
                </div>
                <div>
                    <h3 style="margin: 0; font-size: 1.5rem; font-weight: 700; color: #111827;">${stepNumber}. ${stepName}</h3>
                    <p style="margin: 0.25rem 0 0 0; color: #6b7280; font-size: 0.875rem;">모델: ${result.model_id || 'N/A'}</p>
                </div>
            </div>
            
            <!-- 통계 카드 -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin-bottom: 1.5rem;">
                <!-- 전체 화합물 -->
                <div style="background: linear-gradient(135deg, #f3f4f6 0%, #e5e7eb 100%); padding: 1.25rem; border-radius: 12px; border: 1px solid #d1d5db;">
                    <div style="color: #6b7280; font-size: 0.875rem; font-weight: 500; margin-bottom: 0.5rem;">전체 화합물</div>
                    <div style="font-size: 2rem; font-weight: 700; color: #111827;">${totalCount.toLocaleString()}</div>
                </div>
                
                <!-- 활성 예측 -->
                <div style="background: linear-gradient(135deg, #dcfce7 0%, #bbf7d0 100%); padding: 1.25rem; border-radius: 12px; border: 1px solid #86efac;">
                    <div style="color: #166534; font-size: 0.875rem; font-weight: 500; margin-bottom: 0.5rem;">활성 예측</div>
                    <div style="display: flex; align-items: baseline; gap: 0.5rem;">
                        <span style="font-size: 2rem; font-weight: 700; color: #15803d;">${activeCount.toLocaleString()}</span>
                        <span style="font-size: 1rem; color: #166534;">(${activePercent}%)</span>
                    </div>
                </div>
                
                <!-- 비활성 예측 -->
                <div style="background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%); padding: 1.25rem; border-radius: 12px; border: 1px solid #fca5a5;">
                    <div style="color: #991b1b; font-size: 0.875rem; font-weight: 500; margin-bottom: 0.5rem;">비활성 예측</div>
                    <div style="display: flex; align-items: baseline; gap: 0.5rem;">
                        <span style="font-size: 2rem; font-weight: 700; color: #b91c1c;">${inactiveCount.toLocaleString()}</span>
                        <span style="font-size: 1rem; color: #991b1b;">(${inactivePercent}%)</span>
                    </div>
                </div>
                
                <!-- 임계값 이상 -->
                <div style="background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%); padding: 1.25rem; border-radius: 12px; border: 1px solid #93c5fd;">
                    <div style="color: #1e40af; font-size: 0.875rem; font-weight: 500; margin-bottom: 0.5rem;">임계값 ≥ ${threshold}</div>
                    <div style="font-size: 2rem; font-weight: 700; color: #1d4ed8;">${filteredCount.toLocaleString()}</div>
                </div>
            </div>
            
            <!-- 비율 차트 -->
            <div style="margin-bottom: 1.5rem; padding: 1.25rem; background: white; border-radius: 12px; border: 1px solid #e5e7eb;">
                <div style="font-weight: 600; margin-bottom: 1rem; color: #374151;">활성/비활성 비율</div>
                <div style="display: flex; height: 40px; border-radius: 8px; overflow: hidden; border: 1px solid #e5e7eb;">
                    <div style="background: linear-gradient(135deg, #10b981 0%, #059669 100%); width: ${activePercent}%; display: flex; align-items: center; justify-content: center; color: white; font-weight: 600; font-size: 0.875rem;">
                        ${activePercent > 10 ? `활성 ${activePercent}%` : ''}
                    </div>
                    <div style="background: linear-gradient(135deg, #ef4444 0%, #dc2626 100%); width: ${inactivePercent}%; display: flex; align-items: center; justify-content: center; color: white; font-weight: 600; font-size: 0.875rem;">
                        ${inactivePercent > 10 ? `비활성 ${inactivePercent}%` : ''}
                    </div>
                </div>
            </div>
            
            <!-- 상위 예측 결과 테이블 -->
            ${predictions.length > 0 ? `
                <div style="background: white; border-radius: 12px; border: 1px solid #e5e7eb; overflow: hidden;">
                    <div style="padding: 1rem 1.25rem; background: linear-gradient(135deg, #f9fafb 0%, #f3f4f6 100%); border-bottom: 1px solid #e5e7eb;">
                        <h4 style="margin: 0; font-weight: 600; color: #374151;">🏆 상위 ${predictions.length}개 화합물</h4>
                    </div>
                    <div style="overflow-x: auto;">
                        <table style="width: 100%; border-collapse: collapse;">
                            <thead>
                                <tr style="background: #f9fafb; border-bottom: 2px solid #e5e7eb;">
                                    <th style="padding: 0.75rem 1rem; text-align: left; font-weight: 600; color: #6b7280; font-size: 0.875rem;">순위</th>
                                    <th style="padding: 0.75rem 1rem; text-align: left; font-weight: 600; color: #6b7280; font-size: 0.875rem;">SMILES</th>
                                    <th style="padding: 0.75rem 1rem; text-align: center; font-weight: 600; color: #6b7280; font-size: 0.875rem;">예측</th>
                                    <th style="padding: 0.75rem 1rem; text-align: right; font-weight: 600; color: #6b7280; font-size: 0.875rem;">활성 확률</th>
                                    <th style="padding: 0.75rem 1rem; text-align: right; font-weight: 600; color: #6b7280; font-size: 0.875rem;">신뢰도</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${predictions.map((pred, i) => {
                                    const prob = (pred.active_probability || pred.probability_active || 0) * 100;
                                    const isActive = pred.prediction === 1 || prob >= 50;
                                    const rankBadge = i < 3 ? ['🥇', '🥈', '🥉'][i] : `#${i + 1}`;
                                    
                                    return `
                                        <tr style="border-bottom: 1px solid #f3f4f6; transition: background 0.2s;" onmouseover="this.style.background='#f9fafb'" onmouseout="this.style.background='white'">
                                            <td style="padding: 0.75rem 1rem;">
                                                <span style="font-size: 1.25rem;">${rankBadge}</span>
                                            </td>
                                            <td style="padding: 0.75rem 1rem; font-family: 'Courier New', monospace; font-size: 0.813rem; color: #374151; max-width: 300px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;" title="${pred.smiles || 'N/A'}">
                                                ${pred.smiles || 'N/A'}
                                            </td>
                                            <td style="padding: 0.75rem 1rem; text-align: center;">
                                                <span style="display: inline-block; padding: 0.375rem 0.75rem; border-radius: 6px; font-size: 0.813rem; font-weight: 600; 
                                                    ${isActive ? 'background: linear-gradient(135deg, #dcfce7 0%, #bbf7d0 100%); color: #15803d; border: 1px solid #86efac;' : 'background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%); color: #b91c1c; border: 1px solid #fca5a5;'}">
                                                    ${isActive ? '✓ 활성' : '✗ 비활성'}
                                                </span>
                                            </td>
                                            <td style="padding: 0.75rem 1rem; text-align: right;">
                                                <div style="display: flex; align-items: center; justify-content: flex-end; gap: 0.5rem;">
                                                    <span style="font-family: 'Courier New', monospace; font-weight: 600; font-size: 0.938rem; color: ${prob >= 70 ? '#15803d' : prob >= 50 ? '#ca8a04' : '#6b7280'};">
                                                        ${prob.toFixed(2)}%
                                                    </span>
                                                </div>
                                            </td>
                                            <td style="padding: 0.75rem 1rem; text-align: right;">
                                                <div style="width: 80px; height: 8px; background: #e5e7eb; border-radius: 4px; overflow: hidden;">
                                                    <div style="width: ${prob}%; height: 100%; background: linear-gradient(90deg, 
                                                        ${prob >= 70 ? '#10b981, #059669' : prob >= 50 ? '#f59e0b, #d97706' : '#6b7280, #4b5563'}); 
                                                        border-radius: 4px; transition: width 0.3s ease;"></div>
                                                </div>
                                            </td>
                                        </tr>
                                    `;
                                }).join('')}
                            </tbody>
                        </table>
                    </div>
                </div>
            ` : '<p style="color: #6b7280; text-align: center; padding: 2rem;">예측 결과가 없습니다.</p>'}
        </div>
    `;
}

// SHAP 분석 결과 UI
function createSHAPResultUI(stepNumber, stepName, result) {
    const topFeatures = result.top_features || [];
    const summary = result.shap_values_summary || {};
    const plotPath = result.plot_path || '';
    
    return `
        <div style="padding: 1.5rem;">
            <!-- 헤더 -->
            <div style="display: flex; align-items: center; gap: 1rem; margin-bottom: 1.5rem;">
                <div style="width: 48px; height: 48px; background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); border-radius: 12px; display: flex; align-items: center; justify-content: center;">
                    <i class="fas fa-lightbulb" style="color: white; font-size: 24px;"></i>
                </div>
                <div>
                    <h3 style="margin: 0; font-size: 1.5rem; font-weight: 700; color: #111827;">${stepNumber}. ${stepName}</h3>
                    <p style="margin: 0.25rem 0 0 0; color: #6b7280; font-size: 0.875rem;">모델: ${result.model_id || 'N/A'}</p>
                </div>
            </div>
            
            <!-- SHAP 통계 -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 1rem; margin-bottom: 1.5rem;">
                <div style="background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%); padding: 1rem; border-radius: 12px; border: 1px solid #fbbf24;">
                    <div style="color: #92400e; font-size: 0.75rem; font-weight: 500; margin-bottom: 0.5rem;">평균 |SHAP|</div>
                    <div style="font-size: 1.5rem; font-weight: 700; color: #b45309; font-family: monospace;">${(summary.mean_abs_shap || 0).toFixed(4)}</div>
                </div>
                <div style="background: linear-gradient(135deg, #fed7aa 0%, #fdba74 100%); padding: 1rem; border-radius: 12px; border: 1px solid #fb923c;">
                    <div style="color: #7c2d12; font-size: 0.75rem; font-weight: 500; margin-bottom: 0.5rem;">최대 |SHAP|</div>
                    <div style="font-size: 1.5rem; font-weight: 700; color: #9a3412; font-family: monospace;">${(summary.max_abs_shap || 0).toFixed(4)}</div>
                </div>
                <div style="background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); padding: 1rem; border-radius: 12px; border: 1px solid #a78bfa;">
                    <div style="color: #5b21b6; font-size: 0.75rem; font-weight: 500; margin-bottom: 0.5rem;">분석 샘플</div>
                    <div style="font-size: 1.5rem; font-weight: 700; color: #6d28d9; font-family: monospace;">${summary.samples_analyzed || 0}</div>
                </div>
            </div>
            
            ${topFeatures.length > 0 ? `
                <div style="background: white; border-radius: 12px; border: 1px solid #e5e7eb; overflow: hidden; margin-bottom: 1.5rem;">
                    <div style="padding: 1rem 1.25rem; background: linear-gradient(135deg, #f9fafb 0%, #f3f4f6 100%); border-bottom: 1px solid #e5e7eb;">
                        <h4 style="margin: 0; font-weight: 600; color: #374151;">🔍 상위 20개 중요 특성</h4>
                    </div>
                    <div style="overflow-x: auto; max-height: 400px;">
                        <table style="width: 100%; border-collapse: collapse;">
                            <thead style="position: sticky; top: 0; background: #f9fafb; z-index: 1;">
                                <tr style="border-bottom: 2px solid #e5e7eb;">
                                    <th style="padding: 0.75rem 1rem; text-align: left; font-weight: 600; color: #6b7280; font-size: 0.875rem;">순위</th>
                                    <th style="padding: 0.75rem 1rem; text-align: left; font-weight: 600; color: #6b7280; font-size: 0.875rem;">특성</th>
                                    <th style="padding: 0.75rem 1rem; text-align: right; font-weight: 600; color: #6b7280; font-size: 0.875rem;">중요도</th>
                                    <th style="padding: 0.75rem 1rem; text-align: right; font-weight: 600; color: #6b7280; font-size: 0.875rem;">비중</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${topFeatures.map((feat, i) => {
                                    const importance = feat.importance || 0;
                                    const maxImportance = topFeatures[0]?.importance || 1;
                                    const percent = (importance / maxImportance) * 100;
                                    
                                    return `
                                        <tr style="border-bottom: 1px solid #f3f4f6;" onmouseover="this.style.background='#f9fafb'" onmouseout="this.style.background='white'">
                                            <td style="padding: 0.75rem 1rem;">
                                                <span style="font-weight: 700; color: ${i < 3 ? '#d97706' : '#6b7280'};">
                                                    ${i < 3 ? ['🥇', '🥈', '🥉'][i] : `#${i + 1}`}
                                                </span>
                                            </td>
                                            <td style="padding: 0.75rem 1rem; font-weight: 600; color: #374151;">
                                                ${feat.feature || 'N/A'}
                                            </td>
                                            <td style="padding: 0.75rem 1rem; text-align: right; font-family: 'Courier New', monospace; font-weight: 600; color: #111827;">
                                                ${importance.toFixed(6)}
                                            </td>
                                            <td style="padding: 0.75rem 1rem;">
                                                <div style="display: flex; align-items: center; justify-content: flex-end; gap: 0.5rem;">
                                                    <div style="width: 100px; height: 8px; background: #e5e7eb; border-radius: 4px; overflow: hidden;">
                                                        <div style="width: ${percent}%; height: 100%; background: linear-gradient(90deg, #f59e0b, #d97706); border-radius: 4px;"></div>
                                                    </div>
                                                </div>
                                            </td>
                                        </tr>
                                    `;
                                }).join('')}
                            </tbody>
                        </table>
                    </div>
                </div>
            ` : ''}
            
            ${plotPath ? `
                <div style="background: white; border-radius: 12px; border: 1px solid #e5e7eb; overflow: hidden; padding: 1.25rem;">
                    <h4 style="margin: 0 0 1rem 0; font-weight: 600; color: #374151;">📊 SHAP Summary Plot</h4>
                    <div style="border-radius: 8px; overflow: hidden; border: 1px solid #e5e7eb;">
                        <img src="/${plotPath}" alt="SHAP Plot" style="width: 100%; height: auto; display: block;">
                    </div>
                </div>
            ` : ''}
        </div>
    `;
}

// 모델 학습 결과 UI
function createModelResultUI(stepNumber, stepName, result) {
    const metrics = {
        accuracy: result.accuracy || 0,
        precision: result.precision || 0,
        recall: result.recall || 0,
        f1_score: result.f1_score || 0,
        roc_auc: result.roc_auc || 0
    };
    
    return `
        <div style="padding: 1.5rem;">
            <!-- 헤더 -->
            <div style="display: flex; align-items: center; gap: 1rem; margin-bottom: 1.5rem;">
                <div style="width: 48px; height: 48px; background: linear-gradient(135deg, #ec4899 0%, #be185d 100%); border-radius: 12px; display: flex; align-items: center; justify-content: center;">
                    <i class="fas fa-brain" style="color: white; font-size: 24px;"></i>
                </div>
                <div>
                    <h3 style="margin: 0; font-size: 1.5rem; font-weight: 700; color: #111827;">${stepNumber}. ${stepName}</h3>
                    <p style="margin: 0.25rem 0 0 0; color: #6b7280; font-size: 0.875rem;">모델: ${result.model_id || 'N/A'}</p>
                </div>
            </div>
            
            <!-- 성능 지표 -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(120px, 1fr)); gap: 1rem;">
                ${Object.entries(metrics).map(([key, value]) => {
                    const percent = (value * 100).toFixed(1);
                    const color = percent >= 90 ? '#10b981' : percent >= 70 ? '#f59e0b' : '#ef4444';
                    const bgColor = percent >= 90 ? 'linear-gradient(135deg, #dcfce7 0%, #bbf7d0 100%)' : 
                                   percent >= 70 ? 'linear-gradient(135deg, #fef3c7 0%, #fde68a 100%)' : 
                                   'linear-gradient(135deg, #fee2e2 0%, #fecaca 100%)';
                    
                    return `
                        <div style="background: ${bgColor}; padding: 1.25rem; border-radius: 12px; border: 1px solid ${color}40;">
                            <div style="color: #6b7280; font-size: 0.75rem; font-weight: 500; margin-bottom: 0.5rem; text-transform: uppercase;">
                                ${formatLabel(key)}
                            </div>
                            <div style="font-size: 2rem; font-weight: 700; color: ${color};">
                                ${percent}%
                            </div>
                            <div style="width: 100%; height: 4px; background: #fff; border-radius: 2px; margin-top: 0.5rem; overflow: hidden;">
                                <div style="width: ${percent}%; height: 100%; background: ${color}; border-radius: 2px;"></div>
                            </div>
                        </div>
                    `;
                }).join('')}
            </div>
            
            <!-- 기타 정보 -->
            <div style="margin-top: 1.5rem; padding: 1.25rem; background: #f9fafb; border-radius: 12px; border: 1px solid #e5e7eb;">
                <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem;">
                    ${Object.entries(result).filter(([key]) => !['accuracy', 'precision', 'recall', 'f1_score', 'roc_auc', 'message', 'model_id'].includes(key)).map(([key, value]) => `
                        <div>
                            <div style="color: #6b7280; font-size: 0.813rem; font-weight: 500; margin-bottom: 0.25rem;">
                                ${formatLabel(key)}
                            </div>
                            <div style="color: #111827; font-weight: 600;">
                                ${formatValue(key, value)}
                            </div>
                        </div>
                    `).join('')}
                </div>
            </div>
        </div>
    `;
}

// 라벨 포맷
function formatLabel(key) {
    const labels = {
        'protein_name': '단백질',
        'total_compounds': '총 화합물',
        'active_compounds': '활성 화합물',
        'inactive_compounds': '비활성 화합물',
        'fingerprint_size': 'Fingerprint 크기',
        'descriptor_count': 'Descriptor 수',
        'feature_count': '총 특성 수',
        'model_id': '모델 ID',
        'accuracy': '정확도',
        'precision': '정밀도',
        'recall': '재현율',
        'f1_score': 'F1 Score',
        'roc_auc': 'ROC AUC',
        'count': '예측 수',
        'active_count': '활성 예측',
        'inactive_count': '비활성 예측',
        'top_features': '중요 특성 Top 20',
        'shap_values_summary': 'SHAP 요약',
        'plot_path': 'SHAP 플롯',
        'predictions': '예측 결과',
        'top_predictions': '상위 예측',
        'total_count': '전체 수',
        'threshold': '임계값'
    };
    return labels[key] || key;
}

// 값 포맷
function formatValue(key, value) {
    // SHAP top_features 배열
    if (key === 'top_features' && Array.isArray(value)) {
        return `
            <div style="max-height: 300px; overflow-y: auto; margin-top: 0.5rem;">
                <table style="width: 100%; border-collapse: collapse; font-size: 0.875rem;">
                    <thead>
                        <tr style="background: var(--gray-100); border-bottom: 2px solid var(--gray-300);">
                            <th style="padding: 0.5rem; text-align: left;">#</th>
                            <th style="padding: 0.5rem; text-align: left;">특성</th>
                            <th style="padding: 0.5rem; text-align: right;">중요도</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${value.map((feat, i) => `
                            <tr style="border-bottom: 1px solid var(--gray-200);">
                                <td style="padding: 0.5rem;">${i + 1}</td>
                                <td style="padding: 0.5rem; font-weight: 500;">${feat.feature || 'N/A'}</td>
                                <td style="padding: 0.5rem; text-align: right; font-family: monospace;">${(feat.importance || 0).toFixed(6)}</td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        `;
    }
    
    // SHAP summary 객체
    if (key === 'shap_values_summary' && typeof value === 'object') {
        return `
            <div style="margin-top: 0.5rem;">
                <div style="display: grid; grid-template-columns: auto 1fr; gap: 0.5rem; font-size: 0.875rem;">
                    <span style="font-weight: 500;">평균 |SHAP|:</span>
                    <span style="font-family: monospace;">${(value.mean_abs_shap || 0).toFixed(6)}</span>
                    <span style="font-weight: 500;">최대 |SHAP|:</span>
                    <span style="font-family: monospace;">${(value.max_abs_shap || 0).toFixed(6)}</span>
                    <span style="font-weight: 500;">분석 샘플 수:</span>
                    <span style="font-family: monospace;">${value.samples_analyzed || 0}</span>
                </div>
            </div>
        `;
    }
    
    // SHAP plot_path
    if (key === 'plot_path' && value) {
        return `
            <div style="margin-top: 0.5rem;">
                <a href="/${value}" target="_blank" style="color: var(--primary); text-decoration: underline;">
                    <i class="fas fa-image"></i> ${value}
                </a>
                <div style="margin-top: 1rem; border: 1px solid var(--gray-300); border-radius: 8px; overflow: hidden;">
                    <img src="/${value}" alt="SHAP Plot" style="width: 100%; height: auto; display: block;">
                </div>
            </div>
        `;
    }
    
    // Predictions 배열
    if ((key === 'predictions' || key === 'top_predictions') && Array.isArray(value)) {
        return `
            <div style="max-height: 400px; overflow-y: auto; margin-top: 0.5rem;">
                <table style="width: 100%; border-collapse: collapse; font-size: 0.875rem;">
                    <thead>
                        <tr style="background: var(--gray-100); border-bottom: 2px solid var(--gray-300);">
                            <th style="padding: 0.5rem; text-align: left;">#</th>
                            <th style="padding: 0.5rem; text-align: left;">SMILES</th>
                            <th style="padding: 0.5rem; text-align: center;">예측</th>
                            <th style="padding: 0.5rem; text-align: right;">활성 확률</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${value.map((pred, i) => `
                            <tr style="border-bottom: 1px solid var(--gray-200);">
                                <td style="padding: 0.5rem;">${i + 1}</td>
                                <td style="padding: 0.5rem; font-family: monospace; font-size: 0.75rem; max-width: 300px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;">${pred.smiles || 'N/A'}</td>
                                <td style="padding: 0.5rem; text-align: center;">
                                    <span style="display: inline-block; padding: 0.25rem 0.5rem; border-radius: 4px; font-size: 0.75rem; font-weight: 600; 
                                        ${pred.prediction === 1 ? 'background: #dcfce7; color: #166534;' : 'background: #fee2e2; color: #991b1b;'}">
                                        ${pred.prediction === 1 ? '활성' : '비활성'}
                                    </span>
                                </td>
                                <td style="padding: 0.5rem; text-align: right; font-family: monospace;">
                                    ${((pred.active_probability || 0) * 100).toFixed(2)}%
                                </td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        `;
    }
    
    // 기존 포맷팅
    if (key.includes('accuracy') || key.includes('precision') || key.includes('recall') || 
        key.includes('f1') || key.includes('roc') || key.includes('auc')) {
        return `${(value * 100).toFixed(2)}%`;
    }
    
    if (key === 'threshold' && typeof value === 'number') {
        return value.toFixed(2);
    }
    
    if (typeof value === 'number') {
        return value.toLocaleString();
    }
    
    // 객체/배열이 그대로 넘어온 경우
    if (typeof value === 'object') {
        return JSON.stringify(value, null, 2);
    }
    
    return value;
}

// Step 2: 특성 변환
async function executeFeatureTransform() {
    const proteinName = document.getElementById('protein-name').value;
    const fingerprintType = document.getElementById('fingerprint-type').value;
    const datasetRatio = document.getElementById('dataset-ratio').value;
    const ignore3D = document.getElementById('ignore3d').value === 'true';
    const posThreshold = parseFloat(document.getElementById('pos-threshold').value);
    const negThreshold = parseFloat(document.getElementById('neg-threshold').value);
    
    const response = await axios.post(`${API_BASE_URL}/api/features/transform`, {
        protein_name: proteinName,
        fingerprint_type: fingerprintType,
        dataset_type: 'dataset',
        dataset_ratio: datasetRatio,
        ignore3D: ignore3D,
        pos_threshold: posThreshold,
        neg_threshold: negThreshold
    });
    
    return response.data;
}

// Step 3: 모델 학습
async function executeModelTraining() {
    const proteinName = document.getElementById('protein-name').value;
    const modelType = document.getElementById('model-type').value;
    
    const response = await axios.post(`${API_BASE_URL}/api/models/train`, {
        protein_name: proteinName,
        model_type: modelType,
        feature_type: 'fingerprint',
        test_size: 0.2,
        random_state: 42
    });
    
    // 학습한 모델 ID 저장
    lastTrainedModel = response.data.model_id;
    
    // 모델 선택 드롭다운 업데이트
    await loadModelList();
    
    // AUTO 선택 (방금 학습한 모델 사용)
    document.getElementById('model-select').value = 'AUTO';
    
    return response.data;
}

// Step 4: SHAP 분석
async function executeSHAPAnalysis() {
    // 모델 선택
    let modelId = document.getElementById('model-select').value;
    
    // AUTO이거나 Step 3에서 방금 학습한 경우
    if (modelId === 'AUTO' || !modelId) {
        if (lastTrainedModel) {
            modelId = lastTrainedModel;
        } else {
            throw new Error('사용할 모델이 없습니다. Step 3을 먼저 실행하거나 기존 모델을 선택하세요.');
        }
    }
    
    const response = await axios.post(`${API_BASE_URL}/api/shap/analyze`, {
        model_id: modelId,
        feature_type: 'fingerprint',
        top_n: 20
    });
    
    return response.data;
}

// Step 5: SMILES 직접 입력 예측
async function executeSMILESPrediction() {
    const smilesInput = document.getElementById('smiles-input').value.trim();
    const smilesList = smilesInput.split('\n').map(s => s.trim()).filter(s => s);
    
    // 모델 선택
    let modelId = document.getElementById('model-select').value;
    
    // AUTO이거나 Step 3에서 방금 학습한 경우
    if (modelId === 'AUTO' || !modelId) {
        if (lastTrainedModel) {
            modelId = lastTrainedModel;
        } else {
            throw new Error('사용할 모델이 없습니다. Step 3을 먼저 실행하거나 기존 모델을 선택하세요.');
        }
    }
    
    const response = await axios.post(`${API_BASE_URL}/api/foodb/predict-smiles`, {
        smiles_list: smilesList,
        model_id: modelId
    });
    
    return {
        message: `${response.data.predictions.length}개 화합물 예측 완료`,
        model_id: modelId,
        count: response.data.predictions.length,
        active_count: response.data.predictions.filter(p => p.prediction === 1).length,
        inactive_count: response.data.predictions.filter(p => p.prediction === 0).length,
        predictions: response.data.predictions.slice(0, 10) // 상위 10개만 표시
    };
}

// Step 6: FooDB 전체 예측
async function executeFooDBPrediction() {
    // 모델 선택
    let modelId = document.getElementById('model-select').value;
    
    // AUTO이거나 모델이 없는 경우
    if (modelId === 'AUTO' || !modelId || modelId === '') {
        // 1. Step 3에서 방금 학습한 모델 사용
        if (lastTrainedModel) {
            modelId = lastTrainedModel;
            console.log('Using last trained model:', modelId);
        } else {
            // 2. API에서 가장 최근 모델 가져오기
            try {
                const modelsResponse = await axios.get(`${API_BASE_URL}/api/models/list`);
                const models = modelsResponse.data.models || [];
                
                if (models.length === 0) {
                    throw new Error('사용 가능한 모델이 없습니다. Step 3을 먼저 실행해주세요.');
                }
                
                // 가장 최근 모델 사용 (첫 번째 모델)
                const firstModel = models[0];
                modelId = typeof firstModel === 'string' ? firstModel : firstModel.model_id;
                console.log('Using latest available model:', modelId);
            } catch (error) {
                throw new Error('모델을 찾을 수 없습니다. Step 3을 먼저 실행하거나 기존 모델을 선택해주세요.');
            }
        }
    }
    
    console.log('Final model_id for FooDB prediction:', modelId);
    
    const foodbMode = document.getElementById('foodb-mode').value;
    const topN = parseInt(document.getElementById('foodb-top-n').value) || 100;
    const threshold = parseFloat(document.getElementById('foodb-threshold').value) || 0.5;
    
    // FooDB CSV가 업로드되어 있는지 확인
    // API는 query parameters를 사용함
    
    const response = await axios.post(`${API_BASE_URL}/api/foodb/predict`, null, {
        params: {
            model_id: modelId,
            foodb_file: 'saved_data/FooDB/Compound.csv',  // FooDB 공식 파일명
            batch_size: 100,
            top_n: foodbMode === 'top' ? topN : 10000  // 전체 모드는 큰 값
        }
    });
    
    console.log('FooDB API response:', response.data);
    
    // 백엔드는 top_active_compounds를 리턴함
    const predictions = response.data.top_active_compounds || [];
    console.log('Total predictions received:', predictions.length);
    
    // threshold로 필터링
    const activePredictions = predictions.filter(p => {
        const prob = p.probability_active || p.active_probability || 0;
        return prob >= threshold;
    });
    
    console.log('After threshold filtering:', activePredictions.length);
    
    // SMILES 예측 형식으로 변환
    const formattedPredictions = activePredictions.map(p => ({
        smiles: p.smiles,
        prediction: p.prediction,
        active_probability: p.probability_active || p.active_probability
    }));
    
    // 백엔드 응답에서 전체 통계 가져오기
    const totalCompounds = response.data.total_compounds || predictions.length;
    const activeCount = response.data.active_compounds || 0;
    const inactiveCount = response.data.inactive_compounds || 0;
    
    console.log('Statistics:', {
        total: totalCompounds,
        active: activeCount,
        inactive: inactiveCount,
        filtered: activePredictions.length
    });
    
    return {
        message: `FooDB 예측 완료 (전체: ${totalCompounds}개, 임계값 ${threshold} 이상: ${activePredictions.length}개)`,
        model_id: modelId,
        total_count: totalCompounds,
        active_count: activeCount,
        inactive_count: inactiveCount,
        threshold: threshold,
        filtered_count: activePredictions.length,
        top_predictions: formattedPredictions.slice(0, 20) // 상위 20개만 표시
    };
}

// 모델 목록 로드
async function loadModelList() {
    try {
        const response = await axios.get(`${API_BASE_URL}/api/models/list`);
        const models = response.data.models || [];
        
        const select = document.getElementById('model-select');
        
        // 기본 옵션: AUTO (방금 학습한 모델)
        let options = '<option value="AUTO">자동 (방금 학습한 모델 사용)</option>';
        
        if (models.length > 0) {
            options += '<option value="" disabled>--- 또는 기존 모델 선택 ---</option>';
            // 모델이 객체인 경우 model_id 추출
            options += models.map(m => {
                const modelId = typeof m === 'string' ? m : m.model_id;
                const modelType = typeof m === 'object' && m.model_type ? ` (${m.model_type})` : '';
                const accuracy = typeof m === 'object' && m.metrics?.accuracy ? ` - ${(m.metrics.accuracy * 100).toFixed(1)}%` : '';
                return `<option value="${modelId}">${modelId}${modelType}${accuracy}</option>`;
            }).join('');
        }
        
        select.innerHTML = options;
        
        // 방금 학습한 모델이 있으면 하이라이트
        if (lastTrainedModel) {
            const modelOption = Array.from(select.options).find(opt => opt.value === lastTrainedModel);
            if (modelOption) {
                modelOption.text = `${modelOption.text} ✅ (방금 학습됨)`;
            }
        }
        
    } catch (error) {
        console.error('Failed to load models:', error);
    }
}

// 페이지 로드 시 초기화
document.addEventListener('DOMContentLoaded', function() {
    loadModelList();
    
    // 체크박스 변경 시 카드 스타일 업데이트
    document.querySelectorAll('.step-checkbox').forEach((checkbox, index) => {
        checkbox.addEventListener('change', function() {
            const stepNum = index + 1;
            updateStepCard(stepNum);
            updateStepConfig(stepNum);
            
            // 상태 초기화
            const stepId = this.id;
            const statusElement = document.getElementById(`status${stepNum}`);
            statusElement.className = 'step-status pending';
            statusElement.innerHTML = '<i class="fas fa-circle"></i><span>대기중</span>';
        });
    });
    
    // FooDB 모드 변경 시 상위 N개 옵션 표시/숨김
    const foodbModeSelect = document.getElementById('foodb-mode');
    if (foodbModeSelect) {
        foodbModeSelect.addEventListener('change', function() {
            const topNGroup = document.getElementById('step6-config-top');
            if (topNGroup) {
                topNGroup.style.display = this.value === 'top' ? 'block' : 'none';
            }
        });
    }
});
