"""
PCI Prediction Platform - API Test Script
백엔드 API가 제대로 작동하는지 테스트합니다.
"""

import requests
import json
import time

BASE_URL = 'http://localhost:5000/api'

def print_section(title):
    print(f"\n{'='*60}")
    print(f"🧪 {title}")
    print('='*60)

def test_health():
    """서버 상태 확인"""
    print_section("Health Check")
    response = requests.get(f'{BASE_URL}/health')
    print(f"Status: {response.status_code}")
    print(f"Response: {json.dumps(response.json(), indent=2)}")
    return response.status_code == 200

def test_collect_data():
    """데이터 수집 테스트"""
    print_section("Data Collection")
    payload = {
        'gene_names': ['FTO'],
        'ic50_threshold': 10000
    }
    print(f"Request: {json.dumps(payload, indent=2)}")
    
    response = requests.post(f'{BASE_URL}/collect-data', json=payload)
    print(f"Status: {response.status_code}")
    
    if response.status_code == 200:
        data = response.json()
        print(f"✅ Session ID: {data['session_id']}")
        print(f"✅ Total Compounds: {data['total_compounds']}")
        print(f"✅ ChEMBL: {data['chembl_count']}")
        print(f"✅ BindingDB: {data['bindingdb_count']}")
        return data['session_id']
    else:
        print(f"❌ Error: {response.text}")
        return None

def test_prepare_features(session_id):
    """특성 변환 테스트"""
    print_section("Feature Preparation")
    payload = {
        'session_id': session_id,
        'fingerprint_type': 'ECFP4',
        'fingerprint_size': 1024,
        'active_threshold': 10000,
        'inactive_threshold': 20000,
        'use_decoys': True,
        'decoy_ratio': 50.0,
        'decoy_method': 'dude',
        'include_descriptors': True
    }
    print(f"Request: {json.dumps(payload, indent=2)}")
    
    response = requests.post(f'{BASE_URL}/prepare-features', json=payload)
    print(f"Status: {response.status_code}")
    
    if response.status_code == 200:
        data = response.json()
        print(f"✅ Total Samples: {data['total_samples']}")
        print(f"✅ Actives: {data['n_actives']}")
        print(f"✅ Inactives: {data['n_inactives']}")
        print(f"✅ Features: {data['n_features']}")
        return True
    else:
        print(f"❌ Error: {response.text}")
        return False

def test_train_model(session_id):
    """모델 학습 테스트"""
    print_section("Model Training")
    payload = {
        'session_id': session_id,
        'model_type': 'TabPFN',
        'feature_selection': 'mutual_info',
        'n_features': 100,
        'test_size': 0.2
    }
    print(f"Request: {json.dumps(payload, indent=2)}")
    
    response = requests.post(f'{BASE_URL}/train-model', json=payload)
    print(f"Status: {response.status_code}")
    
    if response.status_code == 200:
        data = response.json()
        print(f"✅ Model: {data['model_type']}")
        print(f"✅ F1 Score: {data['training_metrics']['f1_score']:.4f}")
        print(f"✅ AUC: {data['training_metrics']['auc']:.4f}")
        return True
    else:
        print(f"❌ Error: {response.text}")
        return False

def main():
    print("\n" + "🚀"*30)
    print("PCI Prediction Platform - API Test Suite")
    print("🚀"*30)
    
    # 1. Health Check
    if not test_health():
        print("\n❌ Server is not running! Please start the backend first:")
        print("   cd backend && python app.py")
        return
    
    time.sleep(1)
    
    # 2. Data Collection
    session_id = test_collect_data()
    if not session_id:
        print("\n❌ Data collection failed!")
        return
    
    time.sleep(2)
    
    # 3. Feature Preparation
    if not test_prepare_features(session_id):
        print("\n❌ Feature preparation failed!")
        return
    
    time.sleep(2)
    
    # 4. Model Training
    if not test_train_model(session_id):
        print("\n❌ Model training failed!")
        return
    
    print("\n" + "✅"*30)
    print("All tests passed! 🎉")
    print("✅"*30 + "\n")

if __name__ == '__main__':
    try:
        main()
    except requests.exceptions.ConnectionError:
        print("\n❌ Cannot connect to the server!")
        print("Please start the backend first:")
        print("   cd backend && python app.py")
    except KeyboardInterrupt:
        print("\n\n🛑 Test interrupted by user")
