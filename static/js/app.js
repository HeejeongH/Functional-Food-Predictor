// API 기본 URL
const API_BASE_URL = window.location.origin;

// 페이지 로드 시 초기화
document.addEventListener('DOMContentLoaded', function() {
    // 시스템 상태 확인
    checkSystemStatus();
    
    // 네비게이션 설정
    setupNavigation();
    
    // 모델 목록 로드
    loadModelList();
});

// 네비게이션 설정
function setupNavigation() {
    const navLinks = document.querySelectorAll('.nav-link');
    navLinks.forEach(link => {
        link.addEventListener('click', function(e) {
            e.preventDefault();
            const targetId = this.getAttribute('href').substring(1);
            showSection(targetId);
        });
    });
}

// 섹션 표시
function showSection(sectionId) {
    // 모든 섹션 숨기기
    const sections = document.querySelectorAll('section');
    sections.forEach(section => {
        section.classList.add('hidden');
    });
    
    // 선택된 섹션만 표시
    const targetSection = document.getElementById(sectionId);
    if (targetSection) {
        targetSection.classList.remove('hidden');
        targetSection.classList.add('fade-in');
    }
}

// 시스템 상태 확인
async function checkSystemStatus() {
    try {
        const response = await axios.get(`${API_BASE_URL}/health`);
        const data = response.data;
        
        document.getElementById('status-health').textContent = data.status === 'healthy' ? '정상' : '오류';
        document.getElementById('status-datasets').textContent = data.data.collected_datasets || 0;
        document.getElementById('status-models').textContent = data.data.trained_models || 0;
        
    } catch (error) {
        console.error('Error checking system status:', error);
        document.getElementById('status-health').textContent = '오류';
    }
}

// 데이터 수집
async function collectData() {
    const targetListText = document.getElementById('target-list').value;
    const standardType = document.getElementById('standard-type').value;
    
    if (!targetListText.trim()) {
        showAlert('collection-result', '타겟 유전자 목록을 입력해주세요.', 'error');
        return;
    }
    
    const targetList = targetListText.split('\n').map(t => t.trim()).filter(t => t);
    
    const resultDiv = document.getElementById('collection-result');
    resultDiv.innerHTML = '<div class="text-center"><div class="spinner mx-auto"></div><p class="mt-4">데이터 수집 중...</p></div>';
    
    try {
        const response = await axios.post(`${API_BASE_URL}/api/data/collect`, {
            target_list: targetList,
            standard_type: standardType
        });
        
        const data = response.data;
        resultDiv.innerHTML = `
            <div class="alert alert-success">
                <i class="fas fa-check-circle mr-2"></i>
                ${data.message}
            </div>
            <div class="mt-4 space-y-2">
                <p><strong>ChEMBL 데이터:</strong> ${data.chembl_count}개</p>
                <p><strong>BindingDB 데이터:</strong> ${data.bindingdb_count}개</p>
                <p><strong>총 화합물:</strong> ${data.total_compounds}개</p>
            </div>
        `;
        
        checkSystemStatus();
        
    } catch (error) {
        console.error('Error collecting data:', error);
        resultDiv.innerHTML = `
            <div class="alert alert-error">
                <i class="fas fa-exclamation-circle mr-2"></i>
                데이터 수집 실패: ${error.response?.data?.detail || error.message}
            </div>
        `;
    }
}

// 특성 변환
async function transformFeatures() {
    const proteinName = document.getElementById('protein-name').value;
    const fingerprintType = document.getElementById('fingerprint-type').value;
    const datasetType = document.getElementById('dataset-type').value;
    const datasetRatio = document.getElementById('dataset-ratio').value;
    const ignore3D = document.getElementById('ignore3d').value === 'true';
    const posThreshold = parseFloat(document.getElementById('pos-threshold').value);
    const negThreshold = parseFloat(document.getElementById('neg-threshold').value);
    
    if (!proteinName.trim()) {
        showAlert('transform-result', '단백질 이름을 입력해주세요.', 'error');
        return;
    }
    
    const resultDiv = document.getElementById('transform-result');
    resultDiv.innerHTML = '<div class="text-center"><div class="spinner mx-auto"></div><p class="mt-4">특성 변환 중...</p></div>';
    
    try {
        const response = await axios.post(`${API_BASE_URL}/api/features/transform`, {
            protein_name: proteinName,
            fingerprint_type: fingerprintType,
            dataset_type: datasetType,
            dataset_ratio: datasetRatio,
            ignore3D: ignore3D,
            pos_threshold: posThreshold,
            neg_threshold: negThreshold
        });
        
        const data = response.data;
        
        // Fingerprint 크기와 descriptor 수 계산
        const fpSize = data.fingerprint_size || (fingerprintType === 'MACCS' ? 167 : 1024);
        const descriptorCount = data.descriptor_count || (data.feature_count - fpSize);
        
        resultDiv.innerHTML = `
            <div class="alert alert-success">
                <i class="fas fa-check-circle mr-2"></i>
                ${data.message}
            </div>
            <div class="mt-4 space-y-2">
                <p><strong>단백질:</strong> ${data.protein_name}</p>
                <p><strong>Fingerprint:</strong> ${data.fingerprint_type || fingerprintType} (${fpSize}비트)</p>
                <p><strong>데이터셋:</strong> ${data.dataset_type || datasetType} / ${data.dataset_ratio || datasetRatio} / ${data.ignore3D !== undefined ? (data.ignore3D ? '2D' : '3D') : (ignore3D ? '2D' : '3D')}</p>
                <p><strong>총 화합물:</strong> ${data.total_compounds}개</p>
                <p><strong>활성 화합물:</strong> ${data.active_compounds}개</p>
                <p><strong>비활성 화합물:</strong> ${data.inactive_compounds}개</p>
                <p class="text-lg font-bold text-indigo-600">
                    <strong>총 특성 수:</strong> ${data.feature_count}개 
                    (Fingerprint ${fpSize} + Descriptor ${descriptorCount})
                </p>
            </div>
        `;
        
    } catch (error) {
        console.error('Error transforming features:', error);
        resultDiv.innerHTML = `
            <div class="alert alert-error">
                <i class="fas fa-exclamation-circle mr-2"></i>
                특성 변환 실패: ${error.response?.data?.detail || error.message}
            </div>
        `;
    }
}

// 모델 학습
async function trainModel() {
    const proteinName = document.getElementById('train-protein-name').value;
    const modelType = document.getElementById('model-type').value;
    const featureType = document.getElementById('feature-type').value;
    
    if (!proteinName.trim()) {
        showAlert('train-result', '단백질 이름을 입력해주세요.', 'error');
        return;
    }
    
    const resultDiv = document.getElementById('train-result');
    resultDiv.innerHTML = '<div class="text-center"><div class="spinner mx-auto"></div><p class="mt-4">모델 학습 중... (시간이 걸릴 수 있습니다)</p></div>';
    
    try {
        const response = await axios.post(`${API_BASE_URL}/api/models/train`, {
            protein_name: proteinName,
            model_type: modelType,
            feature_type: featureType,
            test_size: 0.2,
            random_state: 42
        });
        
        const data = response.data;
        resultDiv.innerHTML = `
            <div class="alert alert-success">
                <i class="fas fa-check-circle mr-2"></i>
                ${data.message}
            </div>
            <div class="mt-4">
                <p class="font-bold mb-2">모델 ID: ${data.model_id}</p>
                <div class="grid grid-cols-2 gap-2 text-sm">
                    <p><strong>정확도:</strong> ${(data.accuracy * 100).toFixed(2)}%</p>
                    <p><strong>정밀도:</strong> ${(data.precision * 100).toFixed(2)}%</p>
                    <p><strong>재현율:</strong> ${(data.recall * 100).toFixed(2)}%</p>
                    <p><strong>F1 Score:</strong> ${(data.f1_score * 100).toFixed(2)}%</p>
                    <p><strong>ROC AUC:</strong> ${(data.roc_auc * 100).toFixed(2)}%</p>
                </div>
            </div>
        `;
        
        loadModelList();
        checkSystemStatus();
        
    } catch (error) {
        console.error('Error training model:', error);
        resultDiv.innerHTML = `
            <div class="alert alert-error">
                <i class="fas fa-exclamation-circle mr-2"></i>
                모델 학습 실패: ${error.response?.data?.detail || error.message}
            </div>
        `;
    }
}

// 모델 목록 로드
async function loadModelList() {
    try {
        const response = await axios.get(`${API_BASE_URL}/api/models/list`);
        const data = response.data;
        
        const modelListDiv = document.getElementById('model-list');
        
        if (data.total_models === 0) {
            modelListDiv.innerHTML = `
                <p class="text-gray-500 text-center py-8">
                    <i class="fas fa-box-open mr-2"></i>
                    학습된 모델이 없습니다
                </p>
            `;
            return;
        }
        
        modelListDiv.innerHTML = data.models.map(model => `
            <div class="model-card" onclick="selectModel('${model.model_id}')">
                <div class="flex justify-between items-start mb-2">
                    <h4 class="font-bold text-gray-800">${model.protein_name}</h4>
                    <span class="badge badge-info">${model.model_type}</span>
                </div>
                <p class="text-sm text-gray-600 mb-2">ID: ${model.model_id}</p>
                <div class="grid grid-cols-3 gap-2 text-xs text-gray-500">
                    <div>정확도: ${(model.metrics.accuracy * 100).toFixed(1)}%</div>
                    <div>F1: ${(model.metrics.f1_score * 100).toFixed(1)}%</div>
                    <div>AUC: ${(model.metrics.roc_auc * 100).toFixed(1)}%</div>
                </div>
            </div>
        `).join('');
        
    } catch (error) {
        console.error('Error loading model list:', error);
    }
}

// 모델 선택
function selectModel(modelId) {
    document.getElementById('predict-model-id').value = modelId;
    document.getElementById('shap-model-id').value = modelId;
    showSection('prediction');
}

// 화합물 예측
async function predictCompounds() {
    const modelId = document.getElementById('predict-model-id').value;
    const smilesListText = document.getElementById('smiles-list').value;
    
    if (!modelId.trim()) {
        showAlert('prediction-result', '모델 ID를 입력해주세요.', 'error');
        return;
    }
    
    if (!smilesListText.trim()) {
        showAlert('prediction-result', 'SMILES 목록을 입력해주세요.', 'error');
        return;
    }
    
    const smilesList = smilesListText.split('\n').map(s => s.trim()).filter(s => s);
    
    const resultDiv = document.getElementById('prediction-result');
    resultDiv.innerHTML = '<div class="text-center"><div class="spinner mx-auto"></div><p class="mt-4">예측 중...</p></div>';
    
    try {
        const response = await axios.post(`${API_BASE_URL}/api/models/predict`, {
            smiles_list: smilesList,
            model_id: modelId,
            feature_type: 'fingerprint'
        });
        
        const data = response.data;
        
        resultDiv.innerHTML = `
            <div class="alert alert-success mb-4">
                <i class="fas fa-check-circle mr-2"></i>
                예측 완료!
            </div>
            <div class="space-y-3 max-h-96 overflow-y-auto">
                ${data.predictions.map((pred, idx) => `
                    <div class="prediction-item">
                        <div class="flex justify-between items-start mb-2">
                            <p class="text-sm font-mono text-gray-600">${pred.smiles}</p>
                            <span class="badge ${pred.prediction === 1 ? 'badge-success' : 'badge-danger'}">
                                ${pred.prediction_label}
                            </span>
                        </div>
                        ${pred.probability_active !== undefined ? `
                            <div class="flex items-center space-x-2 text-xs">
                                <span>활성 확률:</span>
                                <div class="flex-1 bg-gray-200 rounded-full h-2">
                                    <div class="bg-green-500 h-2 rounded-full" style="width: ${(pred.probability_active * 100).toFixed(0)}%"></div>
                                </div>
                                <span class="font-bold">${(pred.probability_active * 100).toFixed(1)}%</span>
                            </div>
                        ` : ''}
                    </div>
                `).join('')}
            </div>
        `;
        
    } catch (error) {
        console.error('Error predicting compounds:', error);
        resultDiv.innerHTML = `
            <div class="alert alert-error">
                <i class="fas fa-exclamation-circle mr-2"></i>
                예측 실패: ${error.response?.data?.detail || error.message}
            </div>
        `;
    }
}

// SHAP 분석
async function analyzeSHAP() {
    const modelId = document.getElementById('shap-model-id').value;
    const topN = parseInt(document.getElementById('top-n').value);
    
    if (!modelId.trim()) {
        showAlert('shap-result', '모델 ID를 입력해주세요.', 'error');
        return;
    }
    
    const resultDiv = document.getElementById('shap-result');
    resultDiv.innerHTML = '<div class="text-center"><div class="spinner mx-auto"></div><p class="mt-4">SHAP 분석 중... (시간이 걸릴 수 있습니다)</p></div>';
    
    try {
        const response = await axios.post(`${API_BASE_URL}/api/shap/analyze`, {
            model_id: modelId,
            feature_type: 'fingerprint',
            top_n: topN
        });
        
        const data = response.data;
        
        resultDiv.innerHTML = `
            <div class="alert alert-success mb-4">
                <i class="fas fa-check-circle mr-2"></i>
                SHAP 분석 완료!
            </div>
            <h4 class="font-bold mb-2">상위 ${topN}개 중요 특성</h4>
            <div class="space-y-2 max-h-96 overflow-y-auto">
                ${data.top_features.map((feature, idx) => `
                    <div class="flex items-center space-x-2 text-sm">
                        <span class="font-bold text-gray-600">${idx + 1}.</span>
                        <span class="flex-1 font-mono text-xs">${feature.feature}</span>
                        <div class="flex-1 bg-gray-200 rounded-full h-2">
                            <div class="bg-orange-500 h-2 rounded-full" style="width: ${(feature.importance / data.top_features[0].importance * 100).toFixed(0)}%"></div>
                        </div>
                        <span class="font-bold">${feature.importance.toFixed(4)}</span>
                    </div>
                `).join('')}
            </div>
        `;
        
    } catch (error) {
        console.error('Error analyzing SHAP:', error);
        resultDiv.innerHTML = `
            <div class="alert alert-error">
                <i class="fas fa-exclamation-circle mr-2"></i>
                SHAP 분석 실패: ${error.response?.data?.detail || error.message}
            </div>
        `;
    }
}

// 알림 표시 헬퍼 함수
function showAlert(elementId, message, type = 'info') {
    const element = document.getElementById(elementId);
    const alertClass = `alert-${type}`;
    const icon = {
        success: 'fa-check-circle',
        error: 'fa-exclamation-circle',
        info: 'fa-info-circle',
        warning: 'fa-exclamation-triangle'
    }[type] || 'fa-info-circle';
    
    element.innerHTML = `
        <div class="alert ${alertClass}">
            <i class="fas ${icon} mr-2"></i>
            ${message}
        </div>
    `;
}
