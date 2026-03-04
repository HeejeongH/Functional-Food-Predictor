// API Base URL
const API_BASE_URL = '';

// 워크플로우 실행
async function executeWorkflow() {
    // 선택된 단계 확인
    const steps = [
        { id: 'step1', name: '데이터 수집', fn: executeDataCollection },
        { id: 'step2', name: '특성 변환', fn: executeFeatureTransform },
        { id: 'step3', name: '모델 학습', fn: executeModelTraining },
        { id: 'step4', name: '예측', fn: executePrediction },
        { id: 'step5', name: 'SHAP 분석', fn: executeSHAPAnalysis }
    ];

    const selectedSteps = steps.filter(step => document.getElementById(step.id).checked);

    if (selectedSteps.length === 0) {
        alert('실행할 단계를 최소 1개 이상 선택해주세요.');
        return;
    }

    // 유효성 검사
    const proteinName = document.getElementById('protein-name').value.trim();
    if (!proteinName) {
        alert('단백질 이름을 입력해주세요.');
        return;
    }

    // Step 1 선택 시 타겟 유전자 필수
    if (selectedSteps.find(s => s.id === 'step1')) {
        const targetGenes = document.getElementById('target-genes').value.trim();
        if (!targetGenes) {
            alert('데이터 수집을 위해 타겟 유전자 리스트를 입력해주세요.');
            return;
        }
    }

    // Step 4 선택 시 SMILES 필수
    if (selectedSteps.find(s => s.id === 'step4')) {
        const smiles = document.getElementById('smiles-input').value.trim();
        if (!smiles) {
            alert('예측을 위해 SMILES 리스트를 입력해주세요.');
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
    const statusElement = document.getElementById(`status${stepId.replace('step', '')}`);
    statusElement.className = `step-status ${status}`;
    
    let iconClass = 'fa-circle';
    if (status === 'running') iconClass = 'fa-spinner fa-spin';
    if (status === 'success') iconClass = 'fa-check-circle';
    if (status === 'error') iconClass = 'fa-times-circle';
    
    statusElement.innerHTML = `<i class="fas ${iconClass}"></i>`;
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

// Step 1: 데이터 수집
async function executeDataCollection() {
    const targetGenes = document.getElementById('target-genes').value
        .split(',')
        .map(s => s.trim())
        .filter(s => s);
    
    const response = await axios.post(`${API_BASE_URL}/api/data/collect`, {
        target_list: targetGenes,
        standard_type: 'IC50'
    });
    
    return response.data;
}

// Step 2: 특성 변환
async function executeFeatureTransform() {
    const proteinName = document.getElementById('protein-name').value.trim();
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
    const proteinName = document.getElementById('protein-name').value.trim();
    const modelType = document.getElementById('model-type').value;
    
    const response = await axios.post(`${API_BASE_URL}/api/models/train`, {
        protein_name: proteinName,
        model_type: modelType,
        feature_type: 'fingerprint',
        test_size: 0.2,
        random_state: 42
    });
    
    // 모델 목록 새로고침
    await loadModelList();
    
    return response.data;
}

// Step 4: 예측
async function executePrediction() {
    const smilesInput = document.getElementById('smiles-input').value.trim();
    const smilesList = smilesInput.split('\n').map(s => s.trim()).filter(s => s);
    const modelId = document.getElementById('model-select').value;
    
    if (!modelId) {
        throw new Error('모델을 선택해주세요');
    }
    
    const response = await axios.post(`${API_BASE_URL}/api/models/predict`, {
        smiles_list: smilesList,
        model_id: modelId,
        feature_type: 'fingerprint'
    });
    
    return response.data;
}

// Step 5: SHAP 분석
async function executeSHAPAnalysis() {
    const modelId = document.getElementById('model-select').value;
    
    if (!modelId) {
        throw new Error('모델을 선택해주세요');
    }
    
    const response = await axios.post(`${API_BASE_URL}/api/shap/analyze`, {
        model_id: modelId,
        feature_type: 'fingerprint',
        top_n: 20
    });
    
    return response.data;
}

// 모델 목록 로드
async function loadModelList() {
    try {
        const response = await axios.get(`${API_BASE_URL}/api/models/list`);
        const models = response.data.models || [];
        
        const select = document.getElementById('model-select');
        select.innerHTML = models.length === 0 
            ? '<option value="">학습된 모델이 없습니다</option>'
            : models.map(m => `<option value="${m}">${m}</option>`).join('');
        
    } catch (error) {
        console.error('Failed to load models:', error);
    }
}

// 페이지 로드 시 초기화
document.addEventListener('DOMContentLoaded', function() {
    loadModelList();
    
    // 체크박스 변경 시 상태 초기화
    document.querySelectorAll('.step-checkbox').forEach(checkbox => {
        checkbox.addEventListener('change', function() {
            const stepId = this.id;
            const statusElement = document.getElementById(`status${stepId.replace('step', '')}`);
            statusElement.className = 'step-status pending';
            statusElement.innerHTML = '<i class="fas fa-circle text-gray-300"></i>';
        });
    });
});
