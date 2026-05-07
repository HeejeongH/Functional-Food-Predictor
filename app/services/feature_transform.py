import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.utils.fp_utils import (
    get_fp_config,
    get_fingerprint,
    compute_fingerprints,
    compute_mordred_descriptors,
)
from app.utils.molecular_descriptor import (
    calc_descriptors,
    calc_specific_descriptors,
    remove_invalid_descriptors,
    remove_low_variance_descriptors,
    remove_correlated_descriptors
)
from app.utils.utils import convert_protein_data
import pandas as pd
import numpy as np


# ──────────────────────────────────────────────
# 노트북(TabPFN_F_v3) 기준 컬럼명 헬퍼
# col_prefix='FP', start=0  →  FP0, FP1, ...  (언더바 없음)
# MACCS: FP0 제외 → FP1~FP166
# ──────────────────────────────────────────────

def get_fp_cols(fp_type: str, n_bits: int = 1024) -> list:
    """
    노트북 기준 FP 사용 컬럼명 반환.
    - 비-MACCS : FP0 ~ FP{n_bits-1}
    - MACCS    : FP1 ~ FP166  (bit0은 항상 0이므로 제외)
    """
    _, _, fp_cols = get_fp_config(fp_type.lower(), n_bits=n_bits,
                                   col_prefix='FP', start=0)
    return fp_cols


class FeatureTransformService:

    def __init__(self):
        self.descriptor_selection = self._load_descriptor_selection()

    def _load_descriptor_selection(self) -> dict:
        """
        descriptor_selection.csv → {ratio_key: [col1, col2, ...]} 딕셔너리.
        키는 소문자 ratio (예: '10x', '20x') — 노트북 DESC_BY_RATIO 와 동일.
        """
        csv_path = "descriptor_selection.csv"
        if not os.path.exists(csv_path):
            # 루트가 아닌 경우 상위 탐색
            for alt in ["../descriptor_selection.csv",
                        os.path.join(os.path.dirname(__file__), "../../descriptor_selection.csv")]:
                if os.path.exists(alt):
                    csv_path = alt
                    break

        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            # 노트북과 동일: 컬럼명 소문자로 → ratio key
            return {
                col.lower(): df[col].dropna().tolist()
                for col in df.columns
                if not df[col].dropna().empty
            }
        return {}

    # ── Fingerprint ──────────────────────────────────────────────

    def transform_to_fingerprint(self, smiles_list: list,
                                  fp_type: str = "ecfp4",
                                  fp_size: int = 1024) -> pd.DataFrame:
        """
        SMILES → Fingerprint DataFrame.
        컬럼명: FP0~FP{n-1}  (MACCS는 FP1~FP166)
        fp_utils.compute_fingerprints 를 사용해 chunk 처리 + uint8 최적화.
        """
        fp_type_lower = fp_type.lower()
        n_fp_bits, fp_cols_all, fp_cols = get_fp_config(
            fp_type_lower, n_bits=fp_size, col_prefix='FP', start=0)

        smiles_series = pd.Series(smiles_list)
        fp_df_all = compute_fingerprints(smiles_series, fp_type_lower,
                                         fp_cols_all, n_bits=fp_size)
        return fp_df_all[fp_cols]   # MACCS: FP0 제외

    def get_selected_descriptors(self, ratio: str) -> list:
        """ratio(소문자) 에 해당하는 Weka 선별 descriptor 목록 반환."""
        return self.descriptor_selection.get(ratio.lower(), [])

    def transform_to_fingerprint_with_descriptors(
        self,
        smiles_list: list,
        fp_type: str = "ecfp4",
        dataset_ratio: str = "5x",
        ignore3D: bool = True
    ) -> pd.DataFrame:
        """
        SMILES → Fingerprint + Weka 선별 Descriptor DataFrame.
        노트북 load_data() 와 동일한 feature 구성:
          FP0~FP{n-1}  +  descriptor_selection[dataset_ratio]
        """
        fp_type_lower = fp_type.lower()

        # 1. Fingerprint (FP0~FPn)
        fp_df = self.transform_to_fingerprint(smiles_list, fp_type=fp_type_lower)

        # 2. Weka 선별 descriptor
        selected_descriptors = self.get_selected_descriptors(dataset_ratio)
        if not selected_descriptors:
            print(f"[경고] {dataset_ratio} descriptor 없음 → FP만 반환")
            return fp_df

        # 3. Mordred descriptor 계산 (fp_utils — chunk + uint8)
        try:
            desc_df = compute_mordred_descriptors(
                pd.Series(smiles_list),
                descriptor_names=selected_descriptors,
            )
            # NaN → 0
            desc_df = desc_df.fillna(0)
        except Exception as e:
            print(f"[오류] compute_mordred_descriptors: {e}")
            return fp_df

        # 4. 병합
        result_df = pd.concat(
            [fp_df.reset_index(drop=True), desc_df.reset_index(drop=True)],
            axis=1
        )
        print(f"feature 구성: {fp_df.shape[1]} FP + {desc_df.shape[1]} MD "
              f"= {result_df.shape[1]} 합계")
        return result_df

    def transform_to_descriptors(self, smiles_list: list,
                                   descriptor_type: str = "MORDRED_2D",
                                   descriptor_list: list = None):
        ignore3D = descriptor_type != "MORDRED_3D"
        if descriptor_list:
            desc_df = calc_specific_descriptors(smiles_list,
                                                descriptor_list=descriptor_list,
                                                ignore3D=ignore3D)
        else:
            desc_df = calc_descriptors(smiles_list, ignore3D=ignore3D)

        desc_df, _ = remove_invalid_descriptors(desc_df)
        desc_df_numeric = desc_df.select_dtypes(include=[np.number])
        desc_df_cleaned, _ = remove_low_variance_descriptors(desc_df_numeric)
        desc_df_final, _ = remove_correlated_descriptors(desc_df_cleaned)
        return desc_df_final

    def prepare_training_data(self, chembl_df: pd.DataFrame,
                               bindingdb_df: pd.DataFrame,
                               protein_name: str,
                               dataset_type: str = "fewshot",
                               pos_threshold: float = 10000,
                               neg_threshold: float = 20000):
        return convert_protein_data(
            chembl_df=chembl_df,
            BDB_df=bindingdb_df,
            protein_name=protein_name,
            dataset_type=dataset_type,
            pos_threshold=pos_threshold,
            neg_threshold=neg_threshold
        )
