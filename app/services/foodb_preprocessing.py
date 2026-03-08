import pandas as pd
import numpy as np
import os
import re
import multiprocessing
from typing import List, Dict, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors as mordred_desc
from app.services.feature_transform import FeatureTransformService

class FooDBPreprocessingService:
    """
    FooDB 전처리 서비스
    - 3D Conformer 생성
    - Mordred 3D Descriptor 계산
    - Fingerprint + Descriptor 결합
    """
    
    def __init__(self):
        self.feature_service = FeatureTransformService()
        self.all_descs = list(Calculator(mordred_desc).descriptors)
        self.nproc = max(1, multiprocessing.cpu_count() - 2)
    
    def find_descriptor(self, name: str) -> any:
        """
        Descriptor 이름으로 Mordred descriptor 객체 찾기
        
        Args:
            name: Descriptor 이름 (예: 'RPCG', 'JGI8', 'MATS1p')
        
        Returns:
            Mordred descriptor 객체 또는 None
        """
        # 정확히 일치하는 이름 찾기
        for d in self.all_descs:
            s = str(d)
            if s.endswith(name) or name in s:
                return d
        
        # 패턴 매칭 (예: MATS1p → MATS(1p))
        m = re.match(r'([A-Za-z]+)(\d+)([a-z]*)$', name)
        if m:
            prefix, num, suffix = m.group(1), m.group(2), m.group(3)
            for d in self.all_descs:
                s = str(d)
                pattern = f'({num}{suffix})' if suffix else f'({num})'
                if prefix in s and pattern in s:
                    return d
        
        return None
    
    def generate_3d_conformers(
        self,
        smiles_list: List[str],
        cache_path: str = None
    ) -> Tuple[List, List[int]]:
        """
        SMILES로부터 3D conformer 생성
        
        Args:
            smiles_list: SMILES 리스트
            cache_path: SDF 캐시 파일 경로
        
        Returns:
            (mols, valid_idx) - 3D 분자 객체 리스트와 유효한 인덱스
        """
        
        # 캐시 로드
        if cache_path and os.path.exists(cache_path):
            print(f'  3D mol 캐시 로드: {cache_path}')
            suppl = Chem.SDMolSupplier(cache_path, removeHs=True)
            mols, valid_idx = [], []
            for mol in suppl:
                if mol is not None:
                    mols.append(mol)
                    valid_idx.append(int(mol.GetProp('orig_idx')))
            return mols, valid_idx
        
        # 3D Conformer 생성
        print(f'  3D conformer 생성 중... (시간 소요, {len(smiles_list)}개 분자)')
        mols, valid_idx = [], []
        
        for i, smi in enumerate(smiles_list):
            if i % 100 == 0:
                print(f'    진행: {i}/{len(smiles_list)} ({i/len(smiles_list)*100:.1f}%)')
            
            try:
                mol = Chem.MolFromSmiles(str(smi))
                if mol is None:
                    continue
                
                # 수소 추가 → 3D 생성 → 최적화 → 수소 제거
                mol = Chem.RemoveHs(mol)
                mol = Chem.AddHs(mol)
                result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                
                if result == -1:  # 실패
                    continue
                
                AllChem.MMFFOptimizeMolecule(mol)
                mol = Chem.RemoveHs(mol)
                mol.SetProp('orig_idx', str(i))
                
                mols.append(mol)
                valid_idx.append(i)
            except Exception as e:
                # 에러 무시하고 계속
                pass
        
        print(f'  3D conformer 생성 완료: {len(mols)}/{len(smiles_list)} ({len(mols)/len(smiles_list)*100:.1f}%)')
        
        # 캐시 저장
        if cache_path:
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            writer = Chem.SDWriter(cache_path)
            for mol in mols:
                writer.write(mol)
            writer.close()
            print(f'  3D mol 저장 완료: {cache_path}')
        
        return mols, valid_idx
    
    def calculate_mordred_descriptors(
        self,
        mols: List,
        descriptor_names: List[str],
        ignore3D: bool = False
    ) -> pd.DataFrame:
        """
        Mordred descriptor 계산
        
        Args:
            mols: RDKit 분자 객체 리스트
            descriptor_names: 계산할 descriptor 이름 리스트
            ignore3D: 3D descriptor 무시 여부
        
        Returns:
            Descriptor DataFrame
        """
        # Descriptor 매칭
        descriptor_list, found_names = [], []
        for name in descriptor_names:
            d = self.find_descriptor(name)
            if d:
                descriptor_list.append(d)
                found_names.append(name)
            else:
                print(f'  경고: {name} 매칭 실패')
        
        print(f'  MD {len(descriptor_list)}/{len(descriptor_names)}개 매칭됨')
        
        if not descriptor_list:
            raise ValueError("No descriptors matched")
        
        # Mordred 계산
        calc = Calculator(descriptor_list, ignore_3D=ignore3D)
        desc_df = calc.pandas(mols, nproc=self.nproc)
        desc_df.columns = found_names
        
        # 에러값 처리
        for col in desc_df.columns:
            desc_df[col] = pd.to_numeric(desc_df[col], errors='coerce')
        
        return desc_df
    
    def preprocess_foodb(
        self,
        foodb_csv_path: str,
        protein_name: str,
        fingerprint_type: str,
        dataset_ratio: str,
        ignore3D: bool,
        descriptor_names: List[str],
        output_path: str = None,
        cache_3d_path: str = None
    ) -> pd.DataFrame:
        """
        FooDB 전처리 전체 파이프라인
        
        Args:
            foodb_csv_path: FooDB 원본 CSV 경로
            protein_name: 단백질 이름
            fingerprint_type: Fingerprint 타입 (MACCS, ECFP4, MORGAN)
            dataset_ratio: 데이터셋 비율 (5x, 10x, 20x)
            ignore3D: 3D descriptor 무시 여부
            descriptor_names: Descriptor 이름 리스트
            output_path: 출력 CSV 경로
            cache_3d_path: 3D conformer 캐시 경로 (SDF)
        
        Returns:
            전처리된 DataFrame
        """
        
        print(f'\n[FooDB 전처리] 시작')
        print(f'  입력: {foodb_csv_path}')
        print(f'  Fingerprint: {fingerprint_type}')
        print(f'  Ratio: {dataset_ratio}')
        print(f'  Ignore3D: {ignore3D}')
        
        # FooDB 원본 로드
        foodb_raw = pd.read_csv(foodb_csv_path)
        
        # SMILES 컬럼 찾기
        smiles_col = None
        for col in ['canonical_smiles', 'smiles', 'SMILES', 'canonical_SMILES']:
            if col in foodb_raw.columns:
                smiles_col = col
                break
        
        if not smiles_col:
            raise ValueError("FooDB CSV must contain a SMILES column")
        
        print(f'  FooDB 로드: {len(foodb_raw)}개 화합물, SMILES 컬럼: {smiles_col}')
        
        # 1. Fingerprint 계산
        print('\n[Step 1] Fingerprint 계산')
        fp_df = self.feature_service.transform_to_fingerprint(
            smiles_list=foodb_raw[smiles_col].tolist(),
            fp_type=fingerprint_type
        )
        
        # 2. 3D Conformer 생성
        print('\n[Step 2] 3D Conformer 생성')
        if not cache_3d_path:
            cache_3d_path = f'saved_data/FooDB/3d_conformers/{protein_name}_{dataset_ratio}.sdf'
        
        mols, valid_idx = self.generate_3d_conformers(
            smiles_list=foodb_raw[smiles_col].tolist(),
            cache_path=cache_3d_path
        )
        
        print(f'  사용 가능 분자: {len(mols)}/{len(foodb_raw)}개 ({len(mols)/len(foodb_raw)*100:.1f}%)')
        
        # 3. Mordred Descriptor 계산
        print('\n[Step 3] Mordred Descriptor 계산')
        desc_df = self.calculate_mordred_descriptors(
            mols=mols,
            descriptor_names=descriptor_names,
            ignore3D=ignore3D
        )
        desc_df.index = valid_idx
        
        # 4. 병합
        print('\n[Step 4] 데이터 병합')
        
        # Fingerprint 컬럼 이름
        if fingerprint_type == 'MACCS':
            fp_cols = [f'FP{i}' for i in range(167)]
        else:
            fp_cols = [f'FP{i}' for i in range(1024)]
        
        # 메타 컬럼 (SMILES, id, name 등)
        meta_cols = [c for c in foodb_raw.columns if c not in fp_cols]
        
        # 병합
        final_data = pd.concat([
            foodb_raw[meta_cols].iloc[valid_idx].reset_index(drop=True),
            fp_df.iloc[valid_idx].reset_index(drop=True),
            desc_df.reset_index(drop=True)
        ], axis=1)
        
        # 중복 컬럼 제거
        dupes = final_data.columns[final_data.columns.duplicated()].tolist()
        if dupes:
            print(f'  경고: 중복 컬럼 {dupes} → 첫 번째만 유지')
            final_data = final_data.loc[:, ~final_data.columns.duplicated()]
        
        # 데이터 정제
        data_cols = fp_cols + descriptor_names
        for col in data_cols:
            if col in final_data.columns:
                final_data[col] = pd.to_numeric(final_data[col], errors='coerce')
        
        # 중복 및 결측값 제거
        if smiles_col in final_data.columns:
            final_data = final_data.drop_duplicates(subset=[smiles_col])
        
        final_data = final_data.dropna(subset=data_cols).reset_index(drop=True)
        
        print(f'  최종 데이터: {len(final_data)}개 화합물')
        
        # 5. 저장
        if output_path:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            final_data.to_csv(output_path, index=False)
            print(f'  저장 완료: {output_path}')
        
        return final_data
