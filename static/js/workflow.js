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
    const checkbox = document.getElementById(`step${stepNum}`);
    const isChecked = checkbox.checked;
    
    // Step 5: SMILES 입력 폼
    if (stepNum === 5) {
        const step5Config = document.getElementById('step5-config');
        if (step5Config) {
            step5Config.style.display = isChecked ? 'block' : 'none';
        }
    }
    
    // Step 6: FooDB 옵션 폼
    if (stepNum === 6) {
        const step6ConfigMode = document.getElementById('step6-config-mode');
        const step6ConfigTop = document.getElementById('step6-config-top');
        const step6ConfigThreshold = document.getElementById('step6-config-threshold');
        
        if (step6ConfigMode) step6ConfigMode.style.display = isChecked ? 'block' : 'none';
        if (step6ConfigTop) step6ConfigTop.style.display = isChecked ? 'block' : 'none';
        if (step6ConfigThreshold) step6ConfigThreshold.style.display = isChecked ? 'block' : 'none';
    }
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
            const errorMsg = error.response?.data?.detail || error.message || '오류 발생';
            updateProgressItem(progressItem, 'failed', errorMsg);
            updateStepStatus(step.id, 'error');
            
            alert(`${step.name} 실패: ${errorMsg}`);
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
    
    let content = `<h4><i class="fas fa-check-circle text-green-600 mr-2"></i>${stepNumber}. ${stepName} 결과</h4>`;
    content += '<div class="result-grid">';
    
    // 결과 데이터 표시
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
    card.innerHTML = content;
    
    return card;
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
        'inactive_count': '비활성 예측'
    };
    return labels[key] || key;
}

// 값 포맷
function formatValue(key, value) {
    if (key.includes('accuracy') || key.includes('precision') || key.includes('recall') || 
        key.includes('f1') || key.includes('roc') || key.includes('auc')) {
        return `${(value * 100).toFixed(2)}%`;
    }
    if (typeof value === 'number') {
        return value.toLocaleString();
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
    
    // AUTO이거나 Step 3에서 방금 학습한 경우
    if (modelId === 'AUTO' || !modelId) {
        if (lastTrainedModel) {
            modelId = lastTrainedModel;
        } else {
            throw new Error('사용할 모델이 없습니다. Step 3을 먼저 실행하거나 기존 모델을 선택하세요.');
        }
    }
    
    const foodbMode = document.getElementById('foodb-mode').value;
    const topN = parseInt(document.getElementById('foodb-top-n').value) || 100;
    const threshold = parseFloat(document.getElementById('foodb-threshold').value) || 0.5;
    
    // FooDB CSV가 업로드되어 있는지 확인
    // (여기서는 간단히 API 호출만 수행)
    
    const response = await axios.post(`${API_BASE_URL}/api/foodb/predict`, {
        model_id: modelId,
        top_n: foodbMode === 'top' ? topN : null,
        threshold: threshold
    });
    
    const predictions = response.data.predictions || [];
    const activePredictions = predictions.filter(p => p.prediction === 1 || p.active_probability >= threshold);
    
    return {
        message: `FooDB 예측 완료 (${foodbMode === 'top' ? `상위 ${topN}개` : '전체'})`,
        model_id: modelId,
        total_count: predictions.length,
        active_count: activePredictions.length,
        inactive_count: predictions.length - activePredictions.length,
        threshold: threshold,
        top_predictions: activePredictions.slice(0, 20) // 상위 20개만 표시
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
            options += models.map(m => `<option value="${m}">${m}</option>`).join('');
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
