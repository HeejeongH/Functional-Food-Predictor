"""
PCI Prediction API 테스트 예제

이 스크립트는 API의 전체 워크플로우를 보여줍니다.
"""

import requests
import json

# API 베이스 URL
BASE_URL = "http://localhost:3000"
# 또는 공개 URL 사용: 
# BASE_URL = "https://3000-itvs5ceiquncckj2qhd7g-2b54fc91.sandbox.novita.ai"

def test_health_check():
    """서버 상태 확인"""
    print("=" * 60)
    print("1. 서버 상태 확인")
    print("=" * 60)
    
    response = requests.get(f"{BASE_URL}/health")
    print(f"Status Code: {response.status_code}")
    print(f"Response: {json.dumps(response.json(), indent=2)}")
    print()


def test_data_status():
    """데이터 상태 확인"""
    print("=" * 60)
    print("2. 데이터 상태 확인")
    print("=" * 60)
    
    response = requests.get(f"{BASE_URL}/api/data/status")
    print(f"Status Code: {response.status_code}")
    print(f"Response: {json.dumps(response.json(), indent=2)}")
    print()


def test_fingerprint_generation():
    """간단한 Fingerprint 생성 테스트"""
    print("=" * 60)
    print("3. SMILES를 Fingerprint로 변환")
    print("=" * 60)
    
    # 간단한 SMILES 예제
    test_smiles = [
        "CCO",           # 에탄올
        "CC(=O)O",       # 아세트산
        "c1ccccc1"       # 벤젠
    ]
    
    data = {
        "smiles_list": test_smiles,
        "fp_type": "ECFP4"
    }
    
    response = requests.post(
        f"{BASE_URL}/api/features/fingerprint",
        json=data
    )
    
    print(f"Status Code: {response.status_code}")
    result = response.json()
    
    # Fingerprint 결과는 너무 길어서 요약만 출력
    print(f"Fingerprint Type: {result.get('fingerprint_type')}")
    print(f"Count: {result.get('count')}")
    print(f"Fingerprint Size: {result.get('fingerprint_size')}")
    print(f"First Fingerprint (처음 10개 값): {result.get('fingerprints', [[]])[0][:10]}")
    print()


def test_model_list():
    """학습된 모델 리스트 확인"""
    print("=" * 60)
    print("4. 학습된 모델 리스트 확인")
    print("=" * 60)
    
    response = requests.get(f"{BASE_URL}/api/models/list")
    print(f"Status Code: {response.status_code}")
    print(f"Response: {json.dumps(response.json(), indent=2)}")
    print()


def example_full_workflow():
    """
    전체 워크플로우 예제 (주석으로 설명)
    실제 데이터가 필요하므로 코드만 제공
    """
    print("=" * 60)
    print("5. 전체 워크플로우 예제 (실행 안됨 - 참고용)")
    print("=" * 60)
    print("""
    # 1단계: 데이터 수집 (ChEMBL/BindingDB)
    collect_data = {
        "target_list": ["PDE4", "PDE5"],
        "standard_type": "IC50"
    }
    response = requests.post(f"{BASE_URL}/api/data/collect", json=collect_data)
    
    # 2단계: 특성 변환
    transform_data = {
        "protein_name": "PDE4",
        "fingerprint_type": "ECFP4",
        "dataset_type": "transfer",
        "pos_threshold": 10000,
        "neg_threshold": 20000
    }
    response = requests.post(f"{BASE_URL}/api/features/transform", json=transform_data)
    
    # 3단계: 모델 학습
    train_data = {
        "protein_name": "PDE4",
        "model_type": "XGBoost",
        "feature_type": "fingerprint",
        "test_size": 0.2,
        "random_state": 42
    }
    response = requests.post(f"{BASE_URL}/api/models/train", json=train_data)
    model_id = response.json()["model_id"]
    
    # 4단계: 예측
    predict_data = {
        "smiles_list": ["CCO", "CC(=O)Oc1ccccc1C(=O)O"],
        "model_id": model_id,
        "feature_type": "fingerprint"
    }
    response = requests.post(f"{BASE_URL}/api/models/predict", json=predict_data)
    
    # 5단계: SHAP 분석
    shap_data = {
        "model_id": model_id,
        "feature_type": "fingerprint",
        "top_n": 20
    }
    response = requests.post(f"{BASE_URL}/api/shap/analyze", json=shap_data)
    """)
    print()


if __name__ == "__main__":
    print("\n" + "🚀 " * 20)
    print("PCI Prediction API 테스트 시작")
    print("🚀 " * 20 + "\n")
    
    try:
        # 실제 실행 가능한 테스트들
        test_health_check()
        test_data_status()
        test_fingerprint_generation()
        test_model_list()
        
        # 참고용 예제
        example_full_workflow()
        
        print("=" * 60)
        print("✅ 모든 테스트 완료!")
        print("=" * 60)
        print("\n📚 더 많은 테스트를 하려면:")
        print("   브라우저에서 API 문서를 열어보세요:")
        print(f"   {BASE_URL}/docs")
        print()
        
    except requests.exceptions.ConnectionError:
        print("❌ 서버에 연결할 수 없습니다.")
        print("   서버가 실행 중인지 확인하세요:")
        print("   pm2 list")
        print()
    except Exception as e:
        print(f"❌ 에러 발생: {str(e)}")
        print()
