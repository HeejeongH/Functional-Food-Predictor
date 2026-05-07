import pandas as pd
import numpy as np
import os
from typing import List, Tuple
from rdkit import Chem

from app.services.feature_transform import FeatureTransformService, get_fp_cols
from app.utils.fp_utils import (
    get_fp_config,
    compute_fingerprints,
    compute_mordred_descriptors,
)


class FooDBPreprocessingService:
    """
    FooDB 전처리 서비스 — 노트북(TabPFN_F_v3) FooDB 예측 셀과 동일한 파이프라인.

    변경 사항 (구버전 대비):
    - Fingerprint 컬럼명: X1~X1024 → FP0~FP1023  (fp_utils.get_fp_config 사용)
    - Mordred 계산: 수동 구현 → fp_utils.compute_mordred_descriptors (청크+체크포인트)
    - 3D conformer 생성: fp_utils 내부에서 처리 (2D/3D 자동 분기)
    """

    RENAME_COLS = {'MATS1pe': 'MATS1p', 'ATSC6pe': 'ATSC6p'}

    def __init__(self):
        self.feature_service = FeatureTransformService()

    # ── FooDB 전처리 전체 파이프라인 ─────────────────────────────
    def preprocess_foodb(
        self,
        foodb_csv_path: str,
        protein_name: str,
        fingerprint_type: str,
        dataset_ratio: str,
        ignore3D: bool,
        descriptor_names: List[str],
        output_path: str = None,
        cache_3d_path: str = None,          # 미사용 (fp_utils 내부 처리)
        descriptor_checkpoint_path: str = None,
    ) -> pd.DataFrame:
        """
        FooDB 전처리 전체 파이프라인.

        Parameters
        ----------
        foodb_csv_path            : FooDB 원본 CSV 경로
        protein_name              : 단백질 이름
        fingerprint_type          : 'ecfp4' | 'ecfp6' | 'fcfp4' | 'fcfp6' | 'maccs' | 'rdkit'
        dataset_ratio             : '5x' | '10x' | '20x' 등
        ignore3D                  : True → 2D descriptor만 계산
        descriptor_names          : Weka 선별 descriptor 이름 리스트
        output_path               : 결과 CSV 저장 경로 (None이면 저장 안 함)
        descriptor_checkpoint_path: mordred 계산 중간 저장 경로 (대용량 시 권장)
        """

        print(f'\n[FooDB 전처리] 시작')
        print(f'  입력: {foodb_csv_path}')
        print(f'  Fingerprint: {fingerprint_type}  Ratio: {dataset_ratio}  ignore3D: {ignore3D}')

        # ── 원본 로드 ──────────────────────────────────────────
        foodb_raw = pd.read_csv(foodb_csv_path)

        smiles_col = None
        for col in ['canonical_SMILES', 'canonical_smiles', 'smiles', 'SMILES']:
            if col in foodb_raw.columns:
                smiles_col = col
                break
        if not smiles_col:
            raise ValueError("FooDB CSV에 SMILES 컬럼이 없습니다.")

        print(f'  화합물 수: {len(foodb_raw)}  SMILES 컬럼: {smiles_col}')

        smiles_series = foodb_raw[smiles_col].reset_index(drop=True)

        # ── Step 1. Fingerprint (FP0~FPn) ─────────────────────
        print('\n[Step 1] Fingerprint 계산')
        fp_type_lower = fingerprint_type.lower()
        n_fp_bits, fp_cols_all, fp_cols = get_fp_config(
            fp_type_lower, n_bits=1024, col_prefix='FP', start=0)

        fp_df = compute_fingerprints(smiles_series, fp_type_lower,
                                     fp_cols_all, n_bits=1024)
        fp_df = fp_df[fp_cols]   # MACCS: FP0 제외
        print(f'  FP 컬럼: {fp_cols[0]}~{fp_cols[-1]} ({len(fp_cols)}개)')

        # ── Step 2. Mordred Descriptor (fp_utils — 청크+체크포인트) ──
        print('\n[Step 2] Mordred Descriptor 계산')
        desc_df = compute_mordred_descriptors(
            smiles_series,
            descriptor_names=descriptor_names,
            checkpoint_path=descriptor_checkpoint_path,
        )
        # Mordred 버전 차이 컬럼명 수정
        desc_df.rename(columns=self.RENAME_COLS, inplace=True)
        desc_df = desc_df.fillna(0)
        print(f'  Descriptor: {desc_df.shape[1]}개')

        # ── Step 3. 병합 ──────────────────────────────────────
        print('\n[Step 3] 데이터 병합')
        meta_cols = [c for c in foodb_raw.columns
                     if c not in fp_cols_all and c != smiles_col]
        final_df = pd.concat([
            foodb_raw[meta_cols].reset_index(drop=True),
            foodb_raw[[smiles_col]].reset_index(drop=True),
            fp_df.reset_index(drop=True),
            desc_df.reset_index(drop=True),
        ], axis=1)

        # 중복 컬럼 제거
        dupes = final_df.columns[final_df.columns.duplicated()].tolist()
        if dupes:
            print(f'  중복 컬럼 제거: {dupes}')
            final_df = final_df.loc[:, ~final_df.columns.duplicated()]

        # 수치 변환 및 결측 제거
        feature_cols = fp_cols + descriptor_names
        for col in feature_cols:
            if col in final_df.columns:
                final_df[col] = pd.to_numeric(final_df[col], errors='coerce')

        if smiles_col in final_df.columns:
            final_df = final_df.drop_duplicates(subset=[smiles_col])
        final_df = final_df.dropna(subset=feature_cols).reset_index(drop=True)

        print(f'  최종 화합물 수: {len(final_df)}')

        # ── Step 4. 저장 ──────────────────────────────────────
        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            final_df.to_csv(output_path, index=False)
            print(f'  저장 완료: {output_path}')

        return final_df
