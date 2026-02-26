import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime

try:
    import shap
    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False
    print("Warning: SHAP not available")

class SHAPAnalysisService:
    
    def __init__(self):
        self.output_dir = "shap_outputs"
        os.makedirs(self.output_dir, exist_ok=True)
    
    def analyze_features(self, model, X: np.ndarray, feature_names: list = None, 
                        top_n: int = 20):
        """
        SHAP을 사용한 특성 중요도 분석
        """
        if not SHAP_AVAILABLE:
            raise ValueError("SHAP is not available. Please install shap package.")
        
        print("Computing SHAP values...")
        
        # SHAP Explainer 생성
        try:
            # Tree-based 모델용
            explainer = shap.TreeExplainer(model)
            shap_values = explainer.shap_values(X)
        except:
            # 일반 모델용
            explainer = shap.KernelExplainer(model.predict, X[:100])  # 샘플링
            shap_values = explainer.shap_values(X)
        
        # SHAP 값이 리스트인 경우 (이진 분류)
        if isinstance(shap_values, list):
            shap_values = shap_values[1]  # Positive class
        
        # 절대값 평균으로 특성 중요도 계산
        feature_importance = np.abs(shap_values).mean(axis=0)
        
        # Feature names 설정
        if feature_names is None:
            feature_names = [f"Feature_{i}" for i in range(X.shape[1])]
        
        # DataFrame 생성
        importance_df = pd.DataFrame({
            'feature': feature_names,
            'importance': feature_importance
        }).sort_values('importance', ascending=False)
        
        # Top N 특성
        top_features = importance_df.head(top_n)
        
        # SHAP Summary Plot 저장
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        plot_path = os.path.join(self.output_dir, f"shap_summary_{timestamp}.png")
        
        try:
            plt.figure(figsize=(12, 8))
            
            # Top N 특성만 선택
            top_indices = importance_df.head(top_n).index.tolist()
            X_top = X[:, top_indices]
            shap_values_top = shap_values[:, top_indices]
            feature_names_top = [feature_names[i] for i in top_indices]
            
            shap.summary_plot(
                shap_values_top, 
                X_top, 
                feature_names=feature_names_top,
                show=False
            )
            plt.tight_layout()
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"SHAP plot saved: {plot_path}")
        except Exception as e:
            print(f"Error saving SHAP plot: {str(e)}")
            plot_path = None
        
        return {
            'top_features': top_features.to_dict('records'),
            'shap_values': shap_values,
            'plot_path': plot_path
        }
    
    def get_feature_contributions(self, shap_values: np.ndarray, 
                                 feature_names: list, 
                                 sample_idx: int = 0):
        """
        특정 샘플에 대한 특성 기여도
        """
        contributions = pd.DataFrame({
            'feature': feature_names,
            'shap_value': shap_values[sample_idx]
        }).sort_values('shap_value', key=abs, ascending=False)
        
        return contributions
