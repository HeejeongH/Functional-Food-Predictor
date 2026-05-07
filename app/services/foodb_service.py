import pandas as pd
import numpy as np
from typing import List, Dict

from app.services.feature_transform import FeatureTransformService
from app.utils.fp_utils import get_fp_config, compute_fingerprints


class FooDBService:
    """
    FooDB 식품 화합물 예측 서비스.
    노트북(TabPFN_F_v3) FooDB/Phyto 예측 셀과 동일한 FS 분기 로직 적용.
    """

    def __init__(self):
        self.feature_service = FeatureTransformService()

    # ── 배치 예측 ────────────────────────────────────────────────
    def predict_batch(
        self,
        model_info: dict,
        smiles_list: List[str],
        batch_size: int = 500,
    ) -> Dict:
        """
        FooDB/Phyto 화합물 배치 예측.

        노트북과 동일한 FS 분기:
          pca    → transformer.transform(X_raw)
          mi / shap / random → X_raw[selected_features]
          rfe    → transformer.transform(X_raw)
          none   → X_raw 전체

        Parameters
        ----------
        model_info  : ModelTrainingService.load_model() 반환값
        smiles_list : SMILES 리스트
        batch_size  : 배치 크기
        """
        model             = model_info['model']
        X_train_columns   = model_info['X_train_columns']    # 학습 시 전체 컬럼 순서
        transformer       = model_info.get('transformer')    # PCA / RFE or None
        features          = model_info.get('selected_features')  # 컬럼 이름 리스트 or None
        fs_method         = model_info.get('fs_method', 'none')
        fingerprint_type  = model_info.get('fingerprint_type', 'ecfp4')
        dataset_ratio     = model_info.get('dataset_ratio', '20x')
        fp_bits           = 1024

        total      = len(smiles_list)
        n_batches  = (total + batch_size - 1) // batch_size
        all_preds  = []
        all_probas = []

        print(f"배치 예측 시작: {total}개 화합물, {n_batches}배치")

        # FP 컬럼명 — col_prefix='FP' (학습 시와 통일)
        _, fp_cols_all, _ = get_fp_config(
            fingerprint_type.lower(), n_bits=fp_bits, col_prefix='FP', start=0)

        for i in range(n_batches):
            s = i * batch_size
            e = min(s + batch_size, total)
            batch = smiles_list[s:e]
            print(f"  배치 {i+1}/{n_batches}: {s}~{e-1}")

            try:
                # ── feature 생성 ───────────────────────────────
                features_df = self.feature_service.transform_to_fingerprint_with_descriptors(
                    smiles_list=batch,
                    fp_type=fingerprint_type,
                    dataset_ratio=dataset_ratio,
                )

                # 학습 컬럼 순서에 맞게 정렬 (없는 컬럼은 0 채움)
                X_raw = features_df.reindex(columns=X_train_columns, fill_value=0).fillna(0)

                # ── FS 변환 (노트북과 동일) ─────────────────────
                X_pred = self._apply_fs_transform(
                    X_raw, fs_method, transformer, features)

                # ── 예측 ──────────────────────────────────────
                preds  = model.predict(X_pred)
                probas = (model.predict_proba(X_pred)
                          if hasattr(model, 'predict_proba')
                          else np.column_stack([1 - preds, preds]))

                all_preds.extend(preds.tolist())
                all_probas.extend(probas.tolist())

            except Exception as exc:
                print(f"  [오류] 배치 {i+1}: {exc}")
                all_preds.extend([0] * len(batch))
                all_probas.extend([[1.0, 0.0]] * len(batch))

        probas_arr = np.array(all_probas)
        return {
            'predictions'          : all_preds,
            'probabilities_inactive': probas_arr[:, 0].tolist(),
            'probabilities_active'  : probas_arr[:, 1].tolist(),
        }

    # ── FS 변환 헬퍼 ─────────────────────────────────────────────
    @staticmethod
    def _apply_fs_transform(
        X_raw: pd.DataFrame,
        fs_method: str,
        transformer,
        features,
    ) -> np.ndarray:
        """
        노트북 FooDB/Phyto 예측 셀과 동일한 FS 분기.
        """
        if fs_method == 'pca':
            if transformer is None:
                raise ValueError("PCA transformer가 저장되지 않았습니다.")
            return transformer.transform(X_raw.values)

        elif fs_method in ('mi', 'shap', 'random'):
            if features is None:
                raise ValueError("selected_features가 저장되지 않았습니다.")
            return X_raw[features].values

        elif fs_method == 'rfe':
            if transformer is None:
                raise ValueError("RFE transformer가 저장되지 않았습니다.")
            return transformer.transform(X_raw.values)

        elif fs_method in ('none', None, ''):
            return X_raw.values

        else:
            raise ValueError(f"알 수 없는 fs_method: {fs_method}")
