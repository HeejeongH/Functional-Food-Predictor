"""
머신러닝 모델 학습 모듈 - TabPFN, AutoML 등 다양한 모델 지원
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, 
    f1_score, roc_auc_score, confusion_matrix
)
from sklearn.feature_selection import (
    mutual_info_classif, RFE, SelectKBest, f_classif
)
from sklearn.decomposition import PCA
import joblib

try:
    from tabpfn import TabPFNClassifier
    from tabpfn.constants import ModelVersion
    TABPFN_AVAILABLE = True
except ImportError:
    TABPFN_AVAILABLE = False
    print("Warning: TabPFN not available")


class ModelTrainer:
    """머신러닝 모델 학습 클래스"""
    
    AVAILABLE_MODELS = {
        'TabPFN': 'TabPFN Classifier (v2.5)',
        'RandomForest': 'Random Forest Classifier',
        'GradientBoosting': 'Gradient Boosting Classifier',
        'SVM': 'Support Vector Machine',
        'LogisticRegression': 'Logistic Regression'
    }
    
    FEATURE_SELECTION_METHODS = {
        'none': 'Use all features',
        'mi': 'Mutual Information',
        'rfe': 'Recursive Feature Elimination',
        'pca': 'Principal Component Analysis',
        'univariate': 'Univariate Feature Selection'
    }
    
    def __init__(
        self, 
        model_type: str = 'TabPFN',
        n_features: Optional[int] = None,
        feature_selection_method: str = 'none',
        random_state: int = 42
    ):
        """
        Args:
            model_type: 사용할 모델 타입
            n_features: 선택할 피처 수
            feature_selection_method: 피처 선택 방법
            random_state: 랜덤 시드
        """
        self.model_type = model_type
        self.n_features = n_features
        self.feature_selection_method = feature_selection_method
        self.random_state = random_state
        self.model = None
        self.feature_selector = None
        self.feature_names = None
        self.metrics = {}
        
    def _create_model(self):
        """모델 생성"""
        if self.model_type == 'TabPFN':
            if not TABPFN_AVAILABLE:
                raise ImportError("TabPFN is not available")
            return TabPFNClassifier.create_default_for_version(ModelVersion.V2_5)
        elif self.model_type == 'RandomForest':
            return RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                random_state=self.random_state,
                n_jobs=-1
            )
        elif self.model_type == 'GradientBoosting':
            return GradientBoostingClassifier(
                n_estimators=100,
                max_depth=5,
                random_state=self.random_state
            )
        elif self.model_type == 'SVM':
            return SVC(
                kernel='rbf',
                probability=True,
                random_state=self.random_state
            )
        elif self.model_type == 'LogisticRegression':
            return LogisticRegression(
                max_iter=1000,
                random_state=self.random_state
            )
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")
    
    def _apply_feature_selection(
        self, 
        X_train: np.ndarray, 
        y_train: np.ndarray,
        X_test: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """피처 선택 적용"""
        if self.feature_selection_method == 'none' or self.n_features is None:
            return X_train, X_test
        
        if self.n_features >= X_train.shape[1]:
            return X_train, X_test
        
        print(f"Applying {self.feature_selection_method} feature selection...")
        
        if self.feature_selection_method == 'mi':
            mi_scores = mutual_info_classif(
                X_train, y_train, 
                random_state=self.random_state, 
                n_jobs=-1
            )
            idx = np.argsort(mi_scores)[-self.n_features:]
            return X_train[:, idx], X_test[:, idx]
        
        elif self.feature_selection_method == 'rfe':
            rf_base = RandomForestClassifier(
                n_estimators=50, 
                max_depth=10, 
                random_state=self.random_state, 
                n_jobs=-1
            )
            self.feature_selector = RFE(
                estimator=rf_base, 
                n_features_to_select=self.n_features, 
                step=50
            )
            X_train_sel = self.feature_selector.fit_transform(X_train, y_train)
            X_test_sel = self.feature_selector.transform(X_test)
            return X_train_sel, X_test_sel
        
        elif self.feature_selection_method == 'pca':
            self.feature_selector = PCA(
                n_components=self.n_features, 
                random_state=self.random_state
            )
            X_train_sel = self.feature_selector.fit_transform(X_train)
            X_test_sel = self.feature_selector.transform(X_test)
            return X_train_sel, X_test_sel
        
        elif self.feature_selection_method == 'univariate':
            self.feature_selector = SelectKBest(
                score_func=f_classif, 
                k=self.n_features
            )
            X_train_sel = self.feature_selector.fit_transform(X_train, y_train)
            X_test_sel = self.feature_selector.transform(X_test)
            return X_train_sel, X_test_sel
        
        return X_train, X_test
    
    def train(
        self, 
        df: pd.DataFrame, 
        test_size: float = 0.2,
        feature_columns: Optional[List[str]] = None
    ) -> Dict:
        """
        모델 학습
        
        Args:
            df: 학습 데이터프레임 (Y 컬럼 포함)
            test_size: 테스트 데이터 비율
            feature_columns: 사용할 피처 컬럼 리스트
        
        Returns:
            학습 결과 메트릭
        """
        # 데이터 준비
        if feature_columns is None:
            feature_columns = [col for col in df.columns if col.startswith(('FP_', 'DESC_'))]
        
        X = df[feature_columns].values
        y = df['Y'].values
        
        self.feature_names = feature_columns
        
        # Train/Test 분할
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, 
            random_state=self.random_state, 
            stratify=y
        )
        
        print(f"Training set: {X_train.shape}, Test set: {X_test.shape}")
        
        # 피처 선택
        X_train_sel, X_test_sel = self._apply_feature_selection(
            X_train, y_train, X_test
        )
        
        print(f"After feature selection: {X_train_sel.shape}")
        
        # 모델 학습
        print(f"Training {self.model_type} model...")
        self.model = self._create_model()
        self.model.fit(X_train_sel, y_train)
        
        # 예측
        y_pred = self.model.predict(X_test_sel)
        y_pred_proba = self.model.predict_proba(X_test_sel)
        
        # 메트릭 계산
        self.metrics = {
            'accuracy': accuracy_score(y_test, y_pred),
            'precision': precision_score(y_test, y_pred, average='weighted'),
            'recall': recall_score(y_test, y_pred, average='weighted'),
            'f1': f1_score(y_test, y_pred, average='weighted'),
            'auc': roc_auc_score(y_test, y_pred_proba[:, 1]),
            'confusion_matrix': confusion_matrix(y_test, y_pred).tolist(),
            'n_features': X_train_sel.shape[1],
            'n_train_samples': X_train.shape[0],
            'n_test_samples': X_test.shape[0]
        }
        
        print(f"Training complete - F1: {self.metrics['f1']:.4f}, AUC: {self.metrics['auc']:.4f}")
        
        return self.metrics
    
    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        예측 수행
        
        Args:
            X: 입력 데이터
        
        Returns:
            (predictions, probabilities)
        """
        if self.model is None:
            raise ValueError("Model not trained yet")
        
        # 피처 선택 적용
        if self.feature_selector is not None:
            X = self.feature_selector.transform(X)
        
        predictions = self.model.predict(X)
        probabilities = self.model.predict_proba(X)
        
        return predictions, probabilities
    
    def save_model(self, filepath: str):
        """모델 저장"""
        model_data = {
            'model': self.model,
            'feature_selector': self.feature_selector,
            'model_type': self.model_type,
            'n_features': self.n_features,
            'feature_selection_method': self.feature_selection_method,
            'feature_names': self.feature_names,
            'metrics': self.metrics
        }
        joblib.dump(model_data, filepath)
        print(f"Model saved to {filepath}")
    
    @classmethod
    def load_model(cls, filepath: str) -> 'ModelTrainer':
        """모델 로드"""
        model_data = joblib.load(filepath)
        
        trainer = cls(
            model_type=model_data['model_type'],
            n_features=model_data['n_features'],
            feature_selection_method=model_data['feature_selection_method']
        )
        trainer.model = model_data['model']
        trainer.feature_selector = model_data['feature_selector']
        trainer.feature_names = model_data['feature_names']
        trainer.metrics = model_data['metrics']
        
        return trainer
