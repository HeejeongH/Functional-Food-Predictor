import pandas as pd
import numpy as np
import pickle
import os
from datetime import datetime
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, 
    f1_score, roc_auc_score, confusion_matrix
)
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier

try:
    from tabpfn import TabPFNClassifier
    TABPFN_AVAILABLE = True
except ImportError:
    TABPFN_AVAILABLE = False
    print("Warning: TabPFN not available")

class ModelTrainingService:
    
    def __init__(self):
        self.models = {}
        self.model_dir = "models_trained"
        os.makedirs(self.model_dir, exist_ok=True)
    
    def get_model(self, model_type: str):
        """
        모델 타입에 따라 ML 모델 반환
        """
        if model_type == "TabPFN":
            if not TABPFN_AVAILABLE:
                raise ValueError("TabPFN is not available. Please install tabpfn package.")
            return TabPFNClassifier(device='cpu', N_ensemble_configurations=32)
        
        elif model_type == "XGBoost":
            return XGBClassifier(
                n_estimators=100,
                max_depth=6,
                learning_rate=0.1,
                random_state=42
            )
        
        elif model_type == "LightGBM":
            return LGBMClassifier(
                n_estimators=100,
                max_depth=6,
                learning_rate=0.1,
                random_state=42,
                verbose=-1
            )
        
        elif model_type == "CatBoost":
            return CatBoostClassifier(
                iterations=100,
                depth=6,
                learning_rate=0.1,
                random_state=42,
                verbose=0
            )
        
        elif model_type == "RandomForest":
            return RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                random_state=42
            )
        
        else:
            raise ValueError(f"Unknown model type: {model_type}")
    
    def train_model(self, data_path: str, model_type: str, protein_name: str,
                   test_size: float = 0.2, random_state: int = 42,
                   fingerprint_type: str = "ECFP4", dataset_ratio: str = "20x", ignore3D: bool = True):
        """
        모델 학습 (fingerprint 메타데이터 저장)
        """
        # 데이터 로드
        df = pd.read_csv(data_path)
        
        # X, y 분리
        y = df['Y'].values
        X = df.drop(['SMILES', 'Y', 'potency'], axis=1, errors='ignore').values
        
        # Train/Test 분리
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_state, stratify=y
        )
        
        # 모델 학습
        model = self.get_model(model_type)
        print(f"Training {model_type} model...")
        model.fit(X_train, y_train)
        
        # 예측
        y_pred = model.predict(X_test)
        y_pred_proba = model.predict_proba(X_test)[:, 1] if hasattr(model, 'predict_proba') else y_pred
        
        # 평가 지표 계산
        metrics = {
            'accuracy': accuracy_score(y_test, y_pred),
            'precision': precision_score(y_test, y_pred, zero_division=0),
            'recall': recall_score(y_test, y_pred, zero_division=0),
            'f1_score': f1_score(y_test, y_pred, zero_division=0),
            'roc_auc': roc_auc_score(y_test, y_pred_proba) if len(np.unique(y_test)) > 1 else 0.0
        }
        
        # 모델 저장
        model_id = f"{protein_name}_{model_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        model_path = os.path.join(self.model_dir, f"{model_id}.pkl")
        
        model_info = {
            'model': model,
            'model_type': model_type,
            'protein_name': protein_name,
            'metrics': metrics,
            'feature_count': X.shape[1],
            'train_size': len(X_train),
            'test_size': len(X_test),
            'fingerprint_type': fingerprint_type,
            'dataset_ratio': dataset_ratio,
            'ignore3D': ignore3D
        }
        
        with open(model_path, 'wb') as f:
            pickle.dump(model_info, f)
        
        print(f"Model saved: {model_path}")
        print(f"Metrics: {metrics}")
        
        return model_id, metrics
    
    def load_model(self, model_id: str):
        """
        저장된 모델 로드
        """
        model_path = os.path.join(self.model_dir, f"{model_id}.pkl")
        
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model not found: {model_id}")
        
        with open(model_path, 'rb') as f:
            model_info = pickle.load(f)
        
        return model_info
    
    def predict(self, model_id: str, X: np.ndarray):
        """
        예측 수행 (feature shape 검증 포함)
        """
        model_info = self.load_model(model_id)
        model = model_info['model']
        expected_features = model_info['feature_count']
        
        # Feature shape 검증
        if X.shape[1] != expected_features:
            raise ValueError(
                f"Feature shape mismatch, expected: {expected_features}, got {X.shape[1]}. "
                f"Please use the same fingerprint_type and settings as training: "
                f"fingerprint_type={model_info.get('fingerprint_type', 'unknown')}, "
                f"dataset_ratio={model_info.get('dataset_ratio', 'unknown')}, "
                f"ignore3D={model_info.get('ignore3D', 'unknown')}"
            )
        
        predictions = model.predict(X)
        
        if hasattr(model, 'predict_proba'):
            probabilities = model.predict_proba(X)
            return predictions, probabilities
        else:
            return predictions, None
