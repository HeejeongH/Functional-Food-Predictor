import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import gc
import signal
from datetime import datetime
from contextlib import contextmanager

try:
    import shap
    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False
    print("Warning: SHAP not available")

class TimeoutException(Exception):
    pass

@contextmanager
def timeout(seconds):
    """타임아웃 컨텍스트 매니저"""
    def signal_handler(signum, frame):
        raise TimeoutException(f"Operation timed out after {seconds} seconds")
    
    # Windows에서는 signal.SIGALRM이 없으므로 조건부 처리
    if hasattr(signal, 'SIGALRM'):
        signal.signal(signal.SIGALRM, signal_handler)
        signal.alarm(seconds)
        try:
            yield
        finally:
            signal.alarm(0)
    else:
        # Windows: 타임아웃 없이 진행 (또는 threading.Timer 사용)
        print("⚠️  Timeout not supported on Windows. Proceeding without timeout.")
        yield

class SHAPAnalysisService:
    
    def __init__(self):
        self.output_dir = "shap_outputs"
        os.makedirs(self.output_dir, exist_ok=True)
    
    def analyze_features(self, model, X: np.ndarray, feature_names: list = None, 
                        top_n: int = 20, max_samples: int = 50, timeout_seconds: int = 300):
        """
        SHAP을 사용한 특성 중요도 분석
        
        Args:
            model: 학습된 모델
            X: 입력 데이터
            feature_names: 특성 이름 리스트
            top_n: 상위 N개 특성
            max_samples: 최대 샘플 개수 (기본값: 50, 메모리 절약)
            timeout_seconds: 타임아웃 (초, 기본값: 300초 = 5분)
        """
        if not SHAP_AVAILABLE:
            raise ValueError("SHAP is not available. Please install shap package.")
        
        print(f"\n{'='*60}")
        print(f"🔍 Starting SHAP Analysis")
        print(f"{'='*60}")
        
        # 샘플 제한 (메모리 부족 방지) - 기본값을 100→50으로 축소
        original_samples = X.shape[0]
        if X.shape[0] > max_samples:
            print(f"⚠️  Reducing samples: {X.shape[0]} → {max_samples} (memory optimization)")
            X = X[:max_samples]
        
        print(f"📊 Dataset: {X.shape[0]} samples × {X.shape[1]} features")
        print(f"⏱️  Timeout: {timeout_seconds}s")
        print(f"")
        
        try:
            # Tree-based 모델용 (빠르고 정확)
            print("[1/3] Initializing TreeExplainer...")
            with timeout(timeout_seconds):
                explainer = shap.TreeExplainer(model, feature_perturbation='interventional')
                print("[2/3] Computing SHAP values...")
                shap_values = explainer.shap_values(X, check_additivity=False)
                print("✅ TreeExplainer succeeded")
            
        except TimeoutException as te:
            print(f"\n❌ SHAP computation timed out: {str(te)}")
            raise ValueError(f"SHAP analysis exceeded {timeout_seconds}s timeout. Try reducing max_samples or simplifying the model.")
        
        except Exception as e:
            print(f"\n⚠️  TreeExplainer failed: {str(e)}")
            print("[1/3] Falling back to LinearExplainer...")
            
            try:
                with timeout(timeout_seconds):
                    explainer = shap.LinearExplainer(model, X)
                    print("[2/3] Computing SHAP values (Linear)...")
                    shap_values = explainer.shap_values(X)
                    print("✅ LinearExplainer succeeded")
            
            except TimeoutException as te:
                print(f"\n❌ LinearExplainer timed out: {str(te)}")
                raise ValueError(f"SHAP analysis exceeded timeout. Try reducing max_samples.")
            
            except Exception as e2:
                print(f"\n⚠️  LinearExplainer failed: {str(e2)}")
                print("[1/3] Using KernelExplainer (slow, last resort)...")
                
                try:
                    # 일반 모델용 (매우 느림, 샘플 더 줄이기)
                    kernel_samples = min(30, X.shape[0])  # 50→30으로 축소
                    background_samples = min(25, X.shape[0] // 2)  # 50→25로 축소
                    print(f"   Using {kernel_samples} samples (background: {background_samples})")
                    
                    with timeout(timeout_seconds):
                        explainer = shap.KernelExplainer(
                            model.predict, 
                            X[:background_samples]
                        )
                        print("[2/3] Computing SHAP values (Kernel, very slow)...")
                        shap_values = explainer.shap_values(X[:kernel_samples])
                        print("✅ KernelExplainer succeeded")
                        X = X[:kernel_samples]  # 샘플 수 조정
                
                except TimeoutException as te:
                    print(f"\n❌ KernelExplainer timed out: {str(te)}")
                    raise ValueError(f"All SHAP methods exceeded timeout. Dataset too large or complex.")
                
                except Exception as e3:
                    print(f"\n❌ All SHAP methods failed: {str(e3)}")
                    raise ValueError(f"SHAP analysis failed with all explainers. Error: {str(e3)}")
        
        # SHAP 값이 리스트인 경우 (이진 분류)
        if isinstance(shap_values, list):
            print("\n📊 Binary classification detected, using positive class (index 1)")
            shap_values = shap_values[1]  # Positive class
        
        print(f"✅ SHAP values computed: shape {shap_values.shape}")
        
        # 메모리 정리
        print("[3/3] Cleaning up memory...")
        gc.collect()
        
        # 절대값 평균으로 특성 중요도 계산
        print("\n📈 Computing feature importance...")
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
        
        print(f"\n🏆 Top {top_n} Features:")
        for i, row in enumerate(top_features.head(5).itertuples(), 1):
            print(f"   {i}. {row.feature}: {row.importance:.6f}")
        if top_n > 5:
            print(f"   ... (showing top 5 of {top_n})")
        
        # SHAP Summary Plot 저장
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        plot_path = os.path.join(self.output_dir, f"shap_summary_{timestamp}.png")
        
        try:
            print("\n🎨 Generating SHAP summary plot...")
            plt.figure(figsize=(10, 6))  # 크기 축소 (12×8 → 10×6)
            
            # Top N 특성만 선택
            top_indices = importance_df.head(top_n).index.tolist()
            X_top = X[:, top_indices]
            shap_values_top = shap_values[:, top_indices]
            feature_names_top = [feature_names[i] for i in top_indices]
            
            shap.summary_plot(
                shap_values_top, 
                X_top, 
                feature_names=feature_names_top,
                show=False,
                max_display=top_n
            )
            plt.tight_layout()
            plt.savefig(plot_path, dpi=100, bbox_inches='tight')  # DPI 낮춤 (150→100)
            plt.close()
            print(f"✅ Plot saved: {plot_path}")
            
        except Exception as e:
            print(f"⚠️  Plot generation failed: {str(e)}")
            print("   (Analysis results are still valid)")
            plot_path = None
        
        # 메모리 정리
        print("\n🧹 Final memory cleanup...")
        try:
            del explainer
        except:
            pass
        gc.collect()
        
        print(f"\n{'='*60}")
        print(f"✅ SHAP Analysis Completed Successfully")
        print(f"{'='*60}\n")
        
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
