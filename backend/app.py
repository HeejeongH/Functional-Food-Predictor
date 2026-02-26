"""
PCI Prediction Platform - Flask REST API Backend
전문적인 웹 인터페이스를 위한 백엔드 서버
"""

from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
import sys
import os
import json
import traceback
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime

# 프로젝트 루트를 Python 경로에 추가
sys.path.insert(0, str(Path(__file__).parent.parent))

from modules.data_collector import collect_pci_data
from modules.feature_extractor import prepare_training_data
from modules.model_trainer import train_model, evaluate_model
from modules.shap_analyzer import SHAPAnalyzer
from modules.foodb_predictor import FooDBPredictor

app = Flask(__name__)
CORS(app)  # React 프론트엔드와의 통신을 위해 CORS 활성화

# 데이터 저장 디렉터리
DATA_DIR = Path(__file__).parent.parent / 'data'
MODELS_DIR = Path(__file__).parent.parent / 'models'
RESULTS_DIR = Path(__file__).parent.parent / 'results'

for directory in [DATA_DIR, MODELS_DIR, RESULTS_DIR]:
    directory.mkdir(exist_ok=True)

# 세션 데이터 저장 (실제 프로덕션에서는 Redis 등 사용)
session_data = {}


@app.route('/api/health', methods=['GET'])
def health_check():
    """서버 상태 확인"""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.now().isoformat(),
        'version': '2.0.0'
    })


@app.route('/api/collect-data', methods=['POST'])
def collect_data():
    """
    Step 1: ChEMBL 및 BindingDB에서 PCI 데이터 수집
    
    Request Body:
    {
        "gene_names": ["FTO", "PPARG"],
        "ic50_threshold": 10000
    }
    """
    try:
        data = request.json
        gene_names = data.get('gene_names', [])
        ic50_threshold = data.get('ic50_threshold', 10000)
        
        if not gene_names:
            return jsonify({'error': 'Gene names are required'}), 400
        
        print(f"🔍 Collecting data for genes: {gene_names}")
        
        # 데이터 수집
        chembl_df, bindingdb_df = collect_pci_data(
            gene_names=gene_names,
            ic50_threshold=ic50_threshold
        )
        
        # 데이터 통합
        all_data = pd.concat([chembl_df, bindingdb_df], ignore_index=True)
        
        # 중복 제거 (SMILES 기준)
        all_data = all_data.drop_duplicates(subset=['canonical_smiles'], keep='first')
        
        # 세션에 저장
        session_id = f"session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        session_data[session_id] = {
            'raw_data': all_data,
            'gene_names': gene_names,
            'ic50_threshold': ic50_threshold
        }
        
        # CSV로 저장
        output_path = DATA_DIR / f"{session_id}_raw_data.csv"
        all_data.to_csv(output_path, index=False)
        
        result = {
            'session_id': session_id,
            'total_compounds': len(all_data),
            'chembl_count': len(chembl_df),
            'bindingdb_count': len(bindingdb_df),
            'unique_smiles': all_data['canonical_smiles'].nunique(),
            'genes': gene_names,
            'preview': all_data.head(10).to_dict('records'),
            'ic50_stats': {
                'min': float(all_data['standard_value'].min()),
                'max': float(all_data['standard_value'].max()),
                'mean': float(all_data['standard_value'].mean()),
                'median': float(all_data['standard_value'].median())
            }
        }
        
        return jsonify(result)
        
    except Exception as e:
        print(f"❌ Error in collect_data: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e), 'traceback': traceback.format_exc()}), 500


@app.route('/api/prepare-features', methods=['POST'])
def prepare_features():
    """
    Step 2: 분자 특성 변환 (Fingerprints + Descriptors)
    
    Request Body:
    {
        "session_id": "session_20240226_123456",
        "fingerprint_type": "ECFP4",
        "fingerprint_size": 1024,
        "active_threshold": 10000,
        "inactive_threshold": 20000,
        "use_decoys": true,
        "decoy_ratio": 50.0,
        "decoy_method": "dude",
        "include_descriptors": true
    }
    """
    try:
        data = request.json
        session_id = data.get('session_id')
        
        if not session_id or session_id not in session_data:
            return jsonify({'error': 'Invalid session ID'}), 400
        
        # 세션 데이터 가져오기
        session = session_data[session_id]
        raw_data = session['raw_data']
        gene_names = session['gene_names']
        
        # 파라미터 추출
        fingerprint_type = data.get('fingerprint_type', 'ECFP4')
        fingerprint_size = data.get('fingerprint_size', 1024)
        active_threshold = data.get('active_threshold', 10000)
        inactive_threshold = data.get('inactive_threshold', 20000)
        use_decoys = data.get('use_decoys', True)
        decoy_ratio = data.get('decoy_ratio', 50.0)
        decoy_method = data.get('decoy_method', 'dude')
        include_descriptors = data.get('include_descriptors', True)
        
        print(f"🧬 Preparing features with {fingerprint_type} (size: {fingerprint_size})")
        
        # 특성 변환
        prepared_df = prepare_training_data(
            chembl_df=raw_data,
            protein_name=gene_names[0] if gene_names else 'Target',
            fingerprint_type=fingerprint_type,
            fingerprint_size=fingerprint_size,
            active_threshold=active_threshold,
            inactive_threshold=inactive_threshold,
            use_decoys=use_decoys,
            decoy_ratio=decoy_ratio,
            decoy_method=decoy_method,
            include_descriptors=include_descriptors
        )
        
        # 세션에 저장
        session['prepared_data'] = prepared_df
        session['feature_params'] = {
            'fingerprint_type': fingerprint_type,
            'fingerprint_size': fingerprint_size,
            'active_threshold': active_threshold,
            'inactive_threshold': inactive_threshold,
            'use_decoys': use_decoys,
            'decoy_ratio': decoy_ratio,
            'decoy_method': decoy_method,
            'include_descriptors': include_descriptors
        }
        
        # CSV로 저장
        output_path = DATA_DIR / f"{session_id}_prepared_features.csv"
        prepared_df.to_csv(output_path, index=False)
        
        # 통계 계산
        n_actives = (prepared_df['Y'] == 1).sum()
        n_inactives = (prepared_df['Y'] == 0).sum()
        n_features = len([col for col in prepared_df.columns if col not in ['smiles', 'Y']])
        
        result = {
            'session_id': session_id,
            'total_samples': len(prepared_df),
            'n_actives': int(n_actives),
            'n_inactives': int(n_inactives),
            'n_features': n_features,
            'class_balance': {
                'active': float(n_actives / len(prepared_df)),
                'inactive': float(n_inactives / len(prepared_df))
            },
            'feature_params': session['feature_params'],
            'preview': prepared_df.head(5).to_dict('records')
        }
        
        return jsonify(result)
        
    except Exception as e:
        print(f"❌ Error in prepare_features: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e), 'traceback': traceback.format_exc()}), 500


@app.route('/api/train-model', methods=['POST'])
def train_model_endpoint():
    """
    Step 3: 머신러닝 모델 학습
    
    Request Body:
    {
        "session_id": "session_20240226_123456",
        "model_type": "TabPFN",
        "feature_selection": "mutual_info",
        "n_features": 100,
        "test_size": 0.2
    }
    """
    try:
        data = request.json
        session_id = data.get('session_id')
        
        if not session_id or session_id not in session_data:
            return jsonify({'error': 'Invalid session ID'}), 400
        
        session = session_data[session_id]
        prepared_df = session.get('prepared_data')
        
        if prepared_df is None:
            return jsonify({'error': 'No prepared data found. Run prepare-features first.'}), 400
        
        # 파라미터 추출
        model_type = data.get('model_type', 'TabPFN')
        feature_selection = data.get('feature_selection', 'mutual_info')
        n_features = data.get('n_features', 100)
        test_size = data.get('test_size', 0.2)
        
        print(f"🤖 Training {model_type} model with {feature_selection} feature selection")
        
        # 특성과 레이블 분리
        X = prepared_df.drop(['smiles', 'Y'], axis=1)
        y = prepared_df['Y']
        
        # 모델 학습
        model, selected_features, train_metrics = train_model(
            X=X,
            y=y,
            model_type=model_type,
            feature_selection_method=feature_selection,
            n_features=n_features,
            test_size=test_size
        )
        
        # 모델 평가
        eval_results = evaluate_model(model, X[selected_features], y, test_size=test_size)
        
        # 세션에 저장
        session['model'] = model
        session['selected_features'] = selected_features
        session['train_metrics'] = train_metrics
        session['eval_results'] = eval_results
        session['model_type'] = model_type
        
        # 모델 저장
        import joblib
        model_path = MODELS_DIR / f"{session_id}_{model_type}_model.pkl"
        joblib.dump({
            'model': model,
            'selected_features': selected_features,
            'feature_params': session['feature_params']
        }, model_path)
        
        result = {
            'session_id': session_id,
            'model_type': model_type,
            'n_selected_features': len(selected_features),
            'training_metrics': {
                'accuracy': float(eval_results['accuracy']),
                'precision': float(eval_results['precision']),
                'recall': float(eval_results['recall']),
                'f1_score': float(eval_results['f1_score']),
                'auc': float(eval_results['auc'])
            },
            'confusion_matrix': eval_results['confusion_matrix'].tolist(),
            'feature_importance': train_metrics.get('feature_importance', [])[:20]  # Top 20
        }
        
        return jsonify(result)
        
    except Exception as e:
        print(f"❌ Error in train_model: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e), 'traceback': traceback.format_exc()}), 500


@app.route('/api/shap-analysis', methods=['POST'])
def shap_analysis():
    """
    Step 4: SHAP 분석으로 중요한 화학 특성 도출
    
    Request Body:
    {
        "session_id": "session_20240226_123456"
    }
    """
    try:
        data = request.json
        session_id = data.get('session_id')
        
        if not session_id or session_id not in session_data:
            return jsonify({'error': 'Invalid session ID'}), 400
        
        session = session_data[session_id]
        model = session.get('model')
        prepared_df = session.get('prepared_data')
        selected_features = session.get('selected_features')
        
        if model is None or prepared_df is None:
            return jsonify({'error': 'Model not trained yet'}), 400
        
        print(f"🔍 Running SHAP analysis...")
        
        # SHAP 분석기 초기화
        X = prepared_df.drop(['smiles', 'Y'], axis=1)[selected_features]
        analyzer = SHAPAnalyzer(model, X)
        
        # SHAP 값 계산
        shap_values = analyzer.calculate_shap_values()
        
        # Feature importance
        feature_importance = analyzer.get_feature_importance()
        
        # Summary plot 생성
        summary_plot_path = RESULTS_DIR / f"{session_id}_shap_summary.png"
        analyzer.plot_summary(save_path=str(summary_plot_path))
        
        # Bar plot 생성
        bar_plot_path = RESULTS_DIR / f"{session_id}_shap_bar.png"
        analyzer.plot_bar(save_path=str(bar_plot_path))
        
        # 세션에 저장
        session['shap_values'] = shap_values
        session['feature_importance'] = feature_importance
        
        # Feature importance를 CSV로 저장
        importance_df = pd.DataFrame({
            'feature': selected_features,
            'importance': feature_importance
        }).sort_values('importance', ascending=False)
        
        importance_path = RESULTS_DIR / f"{session_id}_feature_importance.csv"
        importance_df.to_csv(importance_path, index=False)
        
        result = {
            'session_id': session_id,
            'top_features': importance_df.head(20).to_dict('records'),
            'summary_plot': f"/api/results/{session_id}_shap_summary.png",
            'bar_plot': f"/api/results/{session_id}_shap_bar.png",
            'csv_download': f"/api/results/{session_id}_feature_importance.csv"
        }
        
        return jsonify(result)
        
    except Exception as e:
        print(f"❌ Error in shap_analysis: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e), 'traceback': traceback.format_exc()}), 500


@app.route('/api/predict-foodb', methods=['POST'])
def predict_foodb():
    """
    Step 5: FooDB 데이터로 예측
    
    Request Body:
    {
        "session_id": "session_20240226_123456",
        "prediction_threshold": 0.5
    }
    """
    try:
        data = request.json
        session_id = data.get('session_id')
        prediction_threshold = data.get('prediction_threshold', 0.5)
        
        if not session_id or session_id not in session_data:
            return jsonify({'error': 'Invalid session ID'}), 400
        
        session = session_data[session_id]
        model = session.get('model')
        selected_features = session.get('selected_features')
        feature_params = session.get('feature_params')
        
        if model is None:
            return jsonify({'error': 'Model not trained yet'}), 400
        
        print(f"🍎 Predicting on FooDB...")
        
        # FooDB predictor 초기화
        predictor = FooDBPredictor(
            model=model,
            selected_features=selected_features,
            fingerprint_type=feature_params['fingerprint_type'],
            fingerprint_size=feature_params['fingerprint_size'],
            include_descriptors=feature_params['include_descriptors']
        )
        
        # 예측 수행
        predictions_df = predictor.predict_foodb()
        
        # 활성 화합물 필터링
        active_compounds = predictions_df[
            predictions_df['prediction_proba'] >= prediction_threshold
        ].sort_values('prediction_proba', ascending=False)
        
        # 결과 저장
        output_path = RESULTS_DIR / f"{session_id}_foodb_predictions.csv"
        predictions_df.to_csv(output_path, index=False)
        
        active_path = RESULTS_DIR / f"{session_id}_active_compounds.csv"
        active_compounds.to_csv(active_path, index=False)
        
        result = {
            'session_id': session_id,
            'total_compounds': len(predictions_df),
            'active_compounds': len(active_compounds),
            'active_percentage': float(len(active_compounds) / len(predictions_df) * 100),
            'top_predictions': active_compounds.head(20).to_dict('records'),
            'prediction_distribution': {
                'min': float(predictions_df['prediction_proba'].min()),
                'max': float(predictions_df['prediction_proba'].max()),
                'mean': float(predictions_df['prediction_proba'].mean()),
                'median': float(predictions_df['prediction_proba'].median())
            },
            'csv_download': f"/api/results/{session_id}_foodb_predictions.csv",
            'active_csv_download': f"/api/results/{session_id}_active_compounds.csv"
        }
        
        return jsonify(result)
        
    except Exception as e:
        print(f"❌ Error in predict_foodb: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e), 'traceback': traceback.format_exc()}), 500


@app.route('/api/results/<filename>', methods=['GET'])
def download_result(filename):
    """결과 파일 다운로드"""
    try:
        file_path = RESULTS_DIR / filename
        if file_path.exists():
            return send_file(file_path, as_attachment=True)
        else:
            return jsonify({'error': 'File not found'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/sessions', methods=['GET'])
def list_sessions():
    """모든 세션 목록"""
    sessions = []
    for session_id, data in session_data.items():
        sessions.append({
            'session_id': session_id,
            'genes': data.get('gene_names', []),
            'status': 'completed' if 'model' in data else 'in_progress'
        })
    return jsonify({'sessions': sessions})


if __name__ == '__main__':
    print("🚀 Starting PCI Prediction Platform Backend Server...")
    print("📡 API Server running at: http://localhost:5000")
    print("📚 API Documentation:")
    print("  - Health Check: GET /api/health")
    print("  - Collect Data: POST /api/collect-data")
    print("  - Prepare Features: POST /api/prepare-features")
    print("  - Train Model: POST /api/train-model")
    print("  - SHAP Analysis: POST /api/shap-analysis")
    print("  - Predict FooDB: POST /api/predict-foodb")
    app.run(host='0.0.0.0', port=5000, debug=True)
