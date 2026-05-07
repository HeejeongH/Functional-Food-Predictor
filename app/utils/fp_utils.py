"""
fp_utils.py — FTO 프로젝트 공통 유틸리티
TabPFN_F_v2.py, Pycaret_F_v2.py에서 공유하는 함수 모음
"""

import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys, RDKFingerprint, rdFingerprintGenerator
from sklearn.ensemble import RandomForestClassifier
import shap


# ──────────────────────────────────────────────
# Fingerprint
# ──────────────────────────────────────────────

def get_fp_config(fp_type: str, n_bits: int = 1024,
                  col_prefix: str = 'FP', start: int = 0):
    """
    FP 관련 설정값 반환.

    Parameters
    ----------
    fp_type    : 'maccs' | 'ecfp4' | 'ecfp6' | 'fcfp4' | 'fcfp6' | 'rdkit'
    n_bits     : hash FP 비트 수 (MACCS는 무시)
    col_prefix : 컬럼 이름 접두어 ('FP' → FP0, FP1, ... / 'X' → X1, X2, ...)
    start      : 컬럼 번호 시작 인덱스 (0-based: 'FP' → FP0; 1-based: 'X' → X1)

    Returns
    -------
    n_fp_bits  : 전체 비트 수 (MACCS=167, 나머지=n_bits)
    fp_cols_all: 전체 컬럼 이름 리스트 (길이 = n_fp_bits)
    fp_cols    : 사용 컬럼 이름 리스트 (MACCS는 첫 번째 비트 제거)
    """
    n_fp_bits   = 167 if fp_type == 'maccs' else n_bits
    fp_cols_all = [f'{col_prefix}{i + start}' for i in range(n_fp_bits)]
    fp_cols     = fp_cols_all[1:] if fp_type == 'maccs' else fp_cols_all
    return n_fp_bits, fp_cols_all, fp_cols


def get_fingerprint(smiles: str, fp_type: str, n_bits: int = 1024):
    """
    SMILES → fingerprint 비트 벡터 변환.
    분자 파싱 실패 시 None 리스트 반환.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [None] * (167 if fp_type == 'maccs' else n_bits)

    if fp_type == 'maccs':
        return list(MACCSkeys.GenMACCSKeys(mol))

    if fp_type == 'ecfp4':
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    elif fp_type == 'ecfp6':
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=n_bits)
    elif fp_type == 'fcfp4':
        invgen = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
        gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=n_bits, atomInvariantsGenerator=invgen)
    elif fp_type == 'fcfp6':
        invgen = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
        gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=3, fpSize=n_bits, atomInvariantsGenerator=invgen)
    elif fp_type == 'rdkit':
        return list(RDKFingerprint(mol, fpSize=n_bits))
    else:
        raise ValueError(f'Unknown fp_type: {fp_type}')

    return list(gen.GetFingerprintAsNumPy(mol))


def compute_fingerprints(smiles_series: pd.Series, fp_type: str,
                         fp_cols_all: list, n_bits: int = 1024,
                         chunk_size: int = 2000) -> pd.DataFrame:
    """SMILES 시리즈 전체 → FP DataFrame 변환.

    chunk_size 단위로 나눠 처리하고 uint8로 저장해 메모리를 최소화한다.
    """
    import gc
    chunks = []
    for start in range(0, len(smiles_series), chunk_size):
        sl = smiles_series.iloc[start:start + chunk_size]
        rows = sl.apply(lambda s: get_fingerprint(s, fp_type, n_bits)).tolist()
        chunk_df = pd.DataFrame(rows, columns=fp_cols_all, index=sl.index)
        # None이 있으면 0으로 채운 뒤 uint8로 변환 (int64 대비 8배 절약)
        chunk_df = chunk_df.fillna(0).astype('uint8')
        chunks.append(chunk_df)
        del rows, chunk_df
        gc.collect()
    return pd.concat(chunks)


# ──────────────────────────────────────────────
# Descriptor / Data 로딩
# ──────────────────────────────────────────────

def load_descriptor_selection(path: str = '../data/descriptor_selection.csv') -> dict:
    """
    descriptor_selection.csv → {filename: [col1, col2, ...]} 딕셔너리 반환
    """
    df = pd.read_csv(path)
    return {
        col: df[col].dropna().tolist()
        for col in df.columns
        if not df[col].dropna().empty
    }


def get_md_cols(file_md_list: dict, ratio: str,
                fallback_ratio: str = '5x') -> list:
    """ratio에 해당하는 MD 컬럼 목록 반환 (없으면 fallback_ratio 사용)"""
    key = f'descriptors_filtered_FTO_training_{ratio}_ignore3D_False.csv'
    if key not in file_md_list:
        fallback = f'descriptors_filtered_FTO_training_{fallback_ratio}_ignore3D_False.csv'
        print(f'  [경고] {ratio} descriptor 없음 → {fallback_ratio} 사용')
        key = fallback
    return file_md_list[key]


def build_ratio_datasets(fp_md: pd.DataFrame, ratio_multipliers: dict) -> dict:
    active_df = fp_md[fp_md['source'].isin(['active'])]
    decoy_df  = fp_md[fp_md['source'].isin(['decoy', 'assay_inactive'])]
    n_active  = len(active_df)
    n_decoy   = len(decoy_df)

    ratio_dfs = {}
    for ratio, mult in ratio_multipliers.items():
        n_sample = n_active * mult
        if n_sample > n_decoy:
            print(f'[건너뜀] {ratio}: 필요 {n_sample} > 보유 decoy {n_decoy}')
            continue
        ratio_dfs[ratio] = pd.concat(
            [active_df, decoy_df.sample(n=n_sample, random_state=42)]
        ).reset_index(drop=True)
        out_path = f'../data/preprocessed/filtered_FTO_training_{ratio}.csv'
        print(f'{ratio}: {len(ratio_dfs[ratio])}행')
        if os.path.exists(out_path):
            print(f'  [건너뜀] 이미 존재: {out_path}')
        else:
            ratio_dfs[ratio].to_csv(out_path, index=False)
    return ratio_dfs


# ──────────────────────────────────────────────
# SHAP (버그 수정 버전)
# ──────────────────────────────────────────────

def compute_shap_importance(X_train: pd.DataFrame, y_train,
                            n_sample: int = 1000,
                            random_state: int = 42) -> np.ndarray:
    """
    SHAP TreeExplainer 기반 feature 중요도 계산.
    shap 버전에 따른 반환 shape 차이를 안전하게 처리.

    반환값 shape: (n_features,)  — feature별 평균 |SHAP| 값
    """
    rf = RandomForestClassifier(
        n_estimators=100, max_depth=10,
        random_state=random_state, n_jobs=-1)
    rf.fit(X_train, y_train)

    X_sample  = X_train.sample(n=min(n_sample, X_train.shape[0]),
                               random_state=random_state)
    explainer = shap.TreeExplainer(rf)
    shap_vals = explainer.shap_values(X_sample)

    # shape 처리
    if isinstance(shap_vals, list):
        # 구버전 shap: list of [class0_arr, class1_arr]
        # 각 arr shape: (n_samples, n_features)
        # → class 1 (positive class) 기준으로 feature 중요도 계산
        importance = np.abs(shap_vals[1]).mean(axis=0)          # (n_features,)
    else:
        arr = np.array(shap_vals)
        if arr.ndim == 3:
            # (n_samples, n_features, n_classes)
            importance = arr.mean(axis=0).mean(axis=-1)         # (n_features,)
            importance = np.abs(importance)
        else:
            # (n_samples, n_features)
            importance = np.abs(arr).mean(axis=0)               # (n_features,)

    assert importance.shape == (X_train.shape[1],), (
        f'SHAP importance shape 불일치: {importance.shape} vs (n_features={X_train.shape[1]},)'
    )
    return importance


# ──────────────────────────────────────────────
# Mordred Descriptor Computation
# ──────────────────────────────────────────────

# mordred에서 3D 좌표가 필요한 descriptor 목록
_MORDRED_3D_DESCRIPTORS = frozenset({
    'GeomPetitjeanIndex', 'GeomRadius', 'GeomDiameter', 'GeomShapeIndex',
    'PBF', 'MOMI-X', 'MOMI-Y', 'MOMI-Z',
})


def compute_mordred_descriptors(
    smiles_series: pd.Series,
    descriptor_names: list,
    leave_cores: int = 2,
    checkpoint_path: str = None,
    chunk_size: int = 50000,
) -> pd.DataFrame:
    """
    지정된 mordred descriptor를 SMILES로부터 계산한다.

    2D descriptor는 ignore_3D=True로 빠르게 계산하고,
    3D descriptor(GeomPetitjeanIndex 등)는 3D 좌표 생성 후 별도 계산한다.
    계산 실패 분자는 NaN으로 반환 (호출 측에서 fillna(0) 권장).

    Parameters
    ----------
    smiles_series    : canonical SMILES Series (index 보존)
    descriptor_names : 계산할 descriptor 이름 목록
    leave_cores      : 남겨둘 CPU 코어 수
    checkpoint_path  : 중간 저장 파일 경로. 지정 시 chunk_size 단위로 저장하며,
                       재실행 시 이미 완료된 행을 자동으로 복원한다.
    chunk_size       : 한 번에 처리할 분자 수 (2D 기준, 기본 50000)

    Returns
    -------
    DataFrame(index=smiles_series.index, columns=descriptor_names)
    """
    from mordred import Calculator, descriptors as mordred_all
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from tqdm import tqdm
    import multiprocessing

    nproc = max(1, multiprocessing.cpu_count() - leave_cores)

    need_3d = [n for n in descriptor_names if n in _MORDRED_3D_DESCRIPTORS]
    need_2d = [n for n in descriptor_names if n not in _MORDRED_3D_DESCRIPTORS]

    # Calculator 인스턴스를 통해 str(d)가 'AATSC1c' 같은 정확한 이름을 반환함
    # ignore_3D=True/False 두 버전을 합쳐야 GeomPetitjeanIndex 등 3D descriptor도 포함됨
    _ref_2d = Calculator(mordred_all, ignore_3D=True)
    _ref_3d = Calculator(mordred_all, ignore_3D=False)
    all_descs_map = {str(d): d for d in _ref_2d.descriptors}
    all_descs_map.update({str(d): d for d in _ref_3d.descriptors})

    unknown = [n for n in descriptor_names if n not in all_descs_map]
    if unknown:
        print(f'  [경고] mordred에 없는 descriptor: {unknown}')

    result = pd.DataFrame(
        np.nan, index=smiles_series.index, columns=descriptor_names,
        dtype=float,
    )

    # ── 체크포인트 복원 ────────────────────────────────────────
    start_from = 0
    if checkpoint_path and os.path.exists(checkpoint_path):
        ckpt = pd.read_csv(checkpoint_path, index_col=0)
        common_idx = result.index.intersection(ckpt.index)
        for col in descriptor_names:
            if col in ckpt.columns:
                result.loc[common_idx, col] = ckpt.loc[common_idx, col]
        start_from = len(common_idx)
        print(f'  체크포인트 복원: {start_from:,}개 / 전체 {len(smiles_series):,}개')

    # ── 2D descriptors (청크 단위 계산) ───────────────────────
    if need_2d:
        selected = [all_descs_map[n] for n in need_2d if n in all_descs_map]
        calc_2d = Calculator(selected, ignore_3D=True)

        remaining = smiles_series.iloc[start_from:]
        total_valid_2d = 0

        for chunk_start in range(0, len(remaining), chunk_size):
            chunk = remaining.iloc[chunk_start:chunk_start + chunk_size]
            chunk_num = chunk_start // chunk_size + 1
            total_chunks = (len(remaining) + chunk_size - 1) // chunk_size
            print(f'  2D chunk {chunk_num}/{total_chunks} '
                  f'(행 {start_from + chunk_start:,}~{start_from + chunk_start + len(chunk) - 1:,})')

            mols, valid_idx = [], []
            for i, smi in enumerate(chunk):
                try:
                    mol = Chem.MolFromSmiles(str(smi))
                    if mol is not None:
                        mols.append(Chem.AddHs(mol))
                        valid_idx.append(i)
                except Exception:
                    pass

            if mols:
                tmp = calc_2d.pandas(mols, nproc=nproc)
                tmp.index = [chunk.index[j] for j in valid_idx]
                for col in need_2d:
                    if col in tmp.columns:
                        result.loc[tmp.index, col] = pd.to_numeric(
                            tmp[col], errors='coerce'
                        ).values
                total_valid_2d += len(mols)

            # 청크 완료 후 중간 저장
            if checkpoint_path:
                done_idx = list(smiles_series.index[:start_from + chunk_start + len(chunk)])
                result.loc[done_idx].to_csv(checkpoint_path)
                print(f'    → 체크포인트 저장 ({start_from + chunk_start + len(chunk):,}행)')

        print(f'  2D descriptor 완료: {len(need_2d)}개, 유효 분자 {total_valid_2d}/{len(smiles_series)}')

    # ── 3D descriptors (GeomPetitjeanIndex 등) ─────────────────
    if need_3d:
        selected = [all_descs_map[n] for n in need_3d if n in all_descs_map]
        calc_3d = Calculator(selected, ignore_3D=False)

        mols, valid_idx = [], []
        fail_3d = 0
        for i, smi in enumerate(tqdm(smiles_series, desc='3D embed')):
            try:
                mol = Chem.MolFromSmiles(str(smi))
                if mol is None:
                    continue
                mol = Chem.AddHs(mol)
                if AllChem.EmbedMolecule(mol, randomSeed=42) == 0:
                    try:
                        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
                    except Exception:
                        pass
                    # mordred internally calls RemoveHs+MolToSmiles to name the mol;
                    # skip molecules that fail kekulization at that step
                    try:
                        Chem.MolToSmiles(Chem.RemoveHs(mol, updateExplicitCount=True))
                    except Exception:
                        fail_3d += 1
                        continue
                    mols.append(mol)
                    valid_idx.append(i)
                else:
                    fail_3d += 1
            except Exception:
                pass

        if mols:
            tmp = calc_3d.pandas(mols, nproc=nproc)
            tmp.index = [smiles_series.index[i] for i in valid_idx]
            for col in need_3d:
                if col in tmp.columns:
                    result.loc[tmp.index, col] = pd.to_numeric(
                        tmp[col], errors='coerce'
                    ).values

        print(f'  3D descriptor 완료: {len(need_3d)}개, '
              f'유효 {len(mols)}/{len(smiles_series)} (3D 실패 {fail_3d}개 → NaN)')

    return result[descriptor_names]
