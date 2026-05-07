import pandas as pd
import numpy as np
import joblib
import os
from datetime import datetime
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE, mutual_info_classif
from sklearn.decomposition import PCA
from sklearn.metrics import (
    f1_score, roc_auc_score, matthews_corrcoef,
    accuracy_score, precision_score, recall_score
)
from app.utils.fp_utils import compute_shap_importance

try:
    from tabpfn import TabPFNClassifier
    from tabpfn.constants import ModelVersion
    TABPFN_AVAILABLE = True
except ImportError:
    TABPFN_AVAILABLE = False
    print("Warning: TabPFN not available")

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False

try:
    from lightgbm import LGBMClassifier
    LIGHTGBM_AVAILABLE = True
except ImportError:
    LIGHTGBM_AVAILABLE = False

try:
    from catboost import CatBoostClassifier
    CATBOOST_AVAILABLE = True
except ImportError:
    CATBOOST_AVAILABLE = False


class ModelTrainingService:

    def __init__(self):
        self.models = {}
        self.model_dir = "models_trained"
        os.makedirs(self.model_dir, exist_ok=True)

    # ── 모델 팩토리 ──────────────────────────────────────────────
    def get_model(self, model_type: str):
        """모델 타입에 따라 ML 모델 반환"""
        if model_type == "TabPFN":
            if not TABPFN_AVAILABLE:
                raise ValueError("TabPFN is not installed.")
            # 노트북과 동일: TabPFN V2.5
            return TabPFNClassifier.create_default_for_version(
                ModelVersion.V2_5, ignore_pretraining_limits=True)

        elif model_type == "XGBoost":
            if not XGBOOST_AVAILABLE:
                raise ValueError("XGBoost is not installed.")
            return XGBClassifier(
                n_estimators=100, max_depth=6,
                learning_rate=0.1, random_state=42)

        elif model_type == "LightGBM":
            if not LIGHTGBM_AVAILABLE:
                raise ValueError("LightGBM is not installed.")
            return LGBMClassifier(
                n_estimators=100, max_depth=6,
                learning_rate=0.1, random_state=42, verbose=-1)

        elif model_type == "CatBoost":
            if not CATBOOST_AVAILABLE:
                raise ValueError("CatBoost is not installed.")
            return CatBoostClassifier(
                iterations=100, depth=6,
                learning_rate=0.1, random_state=42, verbose=0)

        elif model_type == "RandomForest":
            return RandomForestClassifier(
                n_estimators=100, max_depth=10, random_state=42)

        else:
            raise ValueError(f"Unknown model type: {model_type}")

    # ── Feature Selection ────────────────────────────────────────
    def apply_feature_selection(
        self,
        X_tr: pd.DataFrame,
        X_te: pd.DataFrame,
        y_tr,
        fs_method: str,
        n_features: int,
        random_state: int = 42,
    ):
        """
        노트북 Step3/Step4 와 동일한 FS 로직.

        Returns
        -------
        Xtr_t, Xte_t  : 변환된 numpy array
        transformer   : PCA / RFE 인스턴스 (없으면 None)
        features      : 선택된 컬럼 이름 리스트 (PCA면 None)
        """
        transformer = None
        features    = None

        if fs_method == 'mi':
            scores   = mutual_info_classif(X_tr, y_tr,
                                           random_state=random_state, n_jobs=-1)
            features = X_tr.columns[np.argsort(scores)[-n_features:]].tolist()
            Xtr_t    = X_tr[features].values
            Xte_t    = X_te[features].values

        elif fs_method == 'rfe':
            rf  = RandomForestClassifier(n_estimators=50, max_depth=10,
                                         random_state=random_state, n_jobs=-1)
            rfe = RFE(estimator=rf, n_features_to_select=n_features, step=50)
            rfe.fit(X_tr, y_tr)
            features    = X_tr.columns[rfe.support_].tolist()
            transformer = rfe
            Xtr_t       = X_tr[features].values
            Xte_t       = X_te[features].values

        elif fs_method == 'shap':
            imp      = compute_shap_importance(X_tr, y_tr,
                                               random_state=random_state)
            features = X_tr.columns[np.argsort(imp)[-n_features:]].tolist()
            Xtr_t    = X_tr[features].values
            Xte_t    = X_te[features].values

        elif fs_method == 'pca':
            pca         = PCA(n_components=n_features, random_state=random_state)
            Xtr_t       = pca.fit_transform(X_tr)
            Xte_t       = pca.transform(X_te)
            transformer = pca

        elif fs_method == 'random':
            np.random.seed(random_state)
            idx      = np.random.choice(X_tr.shape[1], n_features, replace=False)
            features = X_tr.columns[idx].tolist()
            Xtr_t    = X_tr[features].values
            Xte_t    = X_te[features].values

        elif fs_method == 'none' or not fs_method:
            # FS 없음 — 전체 feature 사용
            Xtr_t    = X_tr.values
            Xte_t    = X_te.values
            features = X_tr.columns.tolist()

        else:
            raise ValueError(f"Unknown fs_method: {fs_method}")

        return Xtr_t, Xte_t, transformer, features

    # ── 모델 학습 ────────────────────────────────────────────────
    def train_model(
        self,
        data_path: str,
        model_type: str,
        protein_name: str,
        test_size: float = 0.1,          # 노트북과 동일: 0.1
        random_state: int = 42,
        fingerprint_type: str = "ecfp4",
        dataset_ratio: str = "20x",
        ignore3D: bool = True,
        fs_method: str = "pca",          # 노트북 기본: pca
        n_features: int = 100,           # 노트북 기본: 100
    ):
        """
        모델 학습 — 노트북 Step4 final 셀과 동일한 파이프라인.
        저장 포맷: joblib (pickle 대신)
        """
        df = pd.read_csv(data_path)

        y = df['potency'].values if 'potency' in df.columns else df['Y'].values
        drop_cols = ['SMILES', 'canonical_SMILES', 'Y', 'potency', 'source']
        X = df.drop([c for c in drop_cols if c in df.columns], axis=1)

        if len(X) < 100:
            print(f"⚠ 소규모 데이터셋 ({len(X)}개) — 과적합 주의")
        if X.shape[1] > len(X):
            print(f"⚠ feature({X.shape[1]}) > sample({len(X)}) — 과적합 위험")

        X_train_df, X_test_df, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_state, stratify=y)

        # Feature Selection
        Xtr, Xte, transformer, features = self.apply_feature_selection(
            X_train_df, X_test_df, y_train,
            fs_method=fs_method,
            n_features=n_features,
            random_state=random_state,
        )

        # 학습
        model = self.get_model(model_type)
        print(f"학습 중: {model_type} | FS={fs_method} | n_feat={n_features}")
        model.fit(Xtr, y_train)

        # 평가 — 노트북과 동일: F1(weighted), AUC, MCC
        y_pred  = model.predict(Xte)
        y_proba = model.predict_proba(Xte)[:, 1] if hasattr(model, 'predict_proba') else y_pred

        metrics = {
            'f1' : float(f1_score(y_test, y_pred, average='weighted')),
            'auc': float(roc_auc_score(y_test, y_proba)
                         if len(np.unique(y_test)) > 1 else 0.0),
            'mcc': float(matthews_corrcoef(y_test, y_pred)),
            # 보조 지표
            'accuracy' : float(accuracy_score(y_test, y_pred)),
            'precision': float(precision_score(y_test, y_pred, zero_division=0)),
            'recall'   : float(recall_score(y_test, y_pred, zero_division=0)),
        }

        # 저장 (joblib)
        model_id   = f"{protein_name}_{model_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        model_path = os.path.join(self.model_dir, f"{model_id}.joblib")

        model_info = {
            'model'            : model,
            'model_type'       : model_type,
            'protein_name'     : protein_name,
            'metrics'          : metrics,
            'feature_count'    : Xtr.shape[1],
            'train_size'       : len(Xtr),
            'test_size_n'      : len(Xte),
            'fingerprint_type' : fingerprint_type.lower(),
            'dataset_ratio'    : dataset_ratio,
            'ignore3D'         : ignore3D,
            'fs_method'        : fs_method,
            'n_features'       : n_features,
            'selected_features': features,        # mi/shap/random/rfe → 이름 리스트
            'transformer'      : transformer,     # pca/rfe → 객체
            'X_train_columns'  : X.columns.tolist(),
        }

        joblib.dump(model_info, model_path)
        print(f"모델 저장: {model_path}")
        print(f"지표: F1={metrics['f1']:.4f}  AUC={metrics['auc']:.4f}  MCC={metrics['mcc']:.4f}")

        return model_id, metrics

    # ── 모델 로드 ────────────────────────────────────────────────
    def load_model(self, model_id: str) -> dict:
        """저장된 모델 로드 (joblib 우선, 구버전 pkl 폴백)"""
        joblib_path = os.path.join(self.model_dir, f"{model_id}.joblib")
        pkl_path    = os.path.join(self.model_dir, f"{model_id}.pkl")

        if os.path.exists(joblib_path):
            return joblib.load(joblib_path)
        elif os.path.exists(pkl_path):
            import pickle
            with open(pkl_path, 'rb') as f:
                return pickle.load(f)
        else:
            raise FileNotFoundError(f"모델 없음: {model_id}")

    # ── 예측 ─────────────────────────────────────────────────────
    def predict(self, model_id: str, X: np.ndarray):
        """
        예측 수행.
        X는 이미 FS 변환이 완료된 배열이어야 함.
        feature shape 검증 포함.
        """
        model_info = self.load_model(model_id)
        model      = model_info['model']
        expected   = model_info['feature_count']

        if X.shape[1] != expected:
            raise ValueError(
                f"Feature shape 불일치: 기대={expected}, 입력={X.shape[1]}. "
                f"학습 설정과 동일한 fs_method={model_info.get('fs_method')}, "
                f"n_features={model_info.get('n_features')} 를 사용하세요."
            )

        predictions = model.predict(X)
        probabilities = model.predict_proba(X) if hasattr(model, 'predict_proba') else None
        return predictions, probabilities
