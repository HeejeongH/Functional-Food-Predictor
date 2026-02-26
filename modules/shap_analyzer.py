"""
SHAP 분석 모듈 - 모델 설명가능성 분석
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, List
import shap


class SHAPAnalyzer:
    """SHAP 분석 클래스"""
    
    def __init__(self, model, X_train: np.ndarray, feature_names: Optional[List[str]] = None):
        """
        Args:
            model: 학습된 모델
            X_train: 학습 데이터
            feature_names: 피처 이름 리스트
        """
        self.model = model
        self.X_train = X_train
        self.feature_names = feature_names
        self.explainer = None
        self.shap_values = None
        
    def calculate_shap_values(
        self, 
        X_test: np.ndarray, 
        max_background_samples: int = 100,
        nsamples: int = 100
    ):
        """
        SHAP 값 계산
        
        Args:
            X_test: 테스트 데이터
            max_background_samples: 배경 데이터 샘플 수
            nsamples: SHAP 계산 샘플 수
        """
        print("Calculating SHAP values...")
        
        # 배경 데이터 선택
        n_background = min(max_background_samples, self.X_train.shape[0])
        X_background = self.X_train[:n_background]
        
        # Explainer 생성
        def model_predict(X):
            return self.model.predict_proba(X)
        
        self.explainer = shap.KernelExplainer(model_predict, X_background)
        
        # SHAP 값 계산
        self.shap_values = self.explainer.shap_values(X_test, nsamples=nsamples)
        
        # 다차원 배열 처리
        sv = np.array(self.shap_values)
        if sv.ndim == 3:
            self.shap_values_class1 = sv[:, :, 1]
        else:
            self.shap_values_class1 = sv[1]
        
        print(f"SHAP values calculated for {X_test.shape[0]} samples")
        
    def plot_summary(
        self, 
        X_test: np.ndarray, 
        save_path: Optional[str] = None,
        max_display: int = 20
    ):
        """SHAP Summary Plot 생성"""
        if self.shap_values_class1 is None:
            raise ValueError("SHAP values not calculated yet")
        
        plt.figure(figsize=(12, 8))
        shap.summary_plot(
            self.shap_values_class1, 
            X_test, 
            feature_names=self.feature_names,
            show=False,
            max_display=max_display
        )
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Summary plot saved to {save_path}")
        
        return plt.gcf()
    
    def plot_bar(
        self, 
        X_test: np.ndarray, 
        save_path: Optional[str] = None,
        max_display: int = 20
    ):
        """SHAP Bar Plot 생성"""
        if self.shap_values_class1 is None:
            raise ValueError("SHAP values not calculated yet")
        
        plt.figure(figsize=(10, 8))
        shap.summary_plot(
            self.shap_values_class1, 
            X_test, 
            feature_names=self.feature_names,
            plot_type='bar',
            show=False,
            max_display=max_display
        )
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Bar plot saved to {save_path}")
        
        return plt.gcf()
    
    def get_feature_importance(self) -> pd.DataFrame:
        """피처 중요도 DataFrame 반환"""
        if self.shap_values_class1 is None:
            raise ValueError("SHAP values not calculated yet")
        
        importance_df = pd.DataFrame({
            'feature': self.feature_names if self.feature_names else [f'Feature_{i}' for i in range(self.shap_values_class1.shape[1])],
            'mean_abs_shap': np.abs(self.shap_values_class1).mean(axis=0)
        }).sort_values('mean_abs_shap', ascending=False)
        
        return importance_df
    
    def get_top_features(self, n: int = 20) -> List[str]:
        """상위 N개 중요 피처 반환"""
        importance_df = self.get_feature_importance()
        return importance_df.head(n)['feature'].tolist()


def analyze_model_with_shap(
    model,
    X_train: np.ndarray,
    X_test: np.ndarray,
    feature_names: Optional[List[str]] = None,
    save_dir: Optional[str] = None,
    max_display: int = 20
) -> dict:
    """
    모델의 SHAP 분석 수행
    
    Args:
        model: 학습된 모델
        X_train: 학습 데이터
        X_test: 테스트 데이터
        feature_names: 피처 이름 리스트
        save_dir: 저장 디렉토리
        max_display: 표시할 최대 피처 수
    
    Returns:
        분석 결과 딕셔너리
    """
    analyzer = SHAPAnalyzer(model, X_train, feature_names)
    
    # SHAP 값 계산
    analyzer.calculate_shap_values(X_test)
    
    # 플롯 생성
    if save_dir:
        import os
        os.makedirs(save_dir, exist_ok=True)
        
        summary_plot = analyzer.plot_summary(
            X_test, 
            save_path=f"{save_dir}/shap_summary.png",
            max_display=max_display
        )
        
        bar_plot = analyzer.plot_bar(
            X_test, 
            save_path=f"{save_dir}/shap_bar.png",
            max_display=max_display
        )
    else:
        summary_plot = analyzer.plot_summary(X_test, max_display=max_display)
        bar_plot = analyzer.plot_bar(X_test, max_display=max_display)
    
    # 피처 중요도 추출
    importance_df = analyzer.get_feature_importance()
    top_features = analyzer.get_top_features(max_display)
    
    if save_dir:
        importance_df.to_csv(f"{save_dir}/feature_importance.csv", index=False)
        print(f"Feature importance saved to {save_dir}/feature_importance.csv")
    
    return {
        'importance_df': importance_df,
        'top_features': top_features,
        'summary_plot': summary_plot,
        'bar_plot': bar_plot
    }
