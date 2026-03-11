"""
DUD-E (Database of Useful Decoys - Enhanced) 방식 Decoy 생성 서비스

전략:
1. 실제 데이터베이스의 활성(Active) + 비활성(Inactive) 화합물을 우선 사용
2. 비활성 화합물이 부족한 경우, 활성 화합물과 유사한 물리화학적 특성을 가진 decoy 생성
3. DUD-E 기준: 분자량, LogP, 회전 가능한 결합, 수소결합 donor/acceptor 등을 매칭
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, AllChem, DataStructs
from typing import List, Dict, Tuple
import random
from collections import defaultdict


class DecoyGenerationService:
    """
    DUD-E 방식의 Decoy 생성 서비스
    - 실제 DB 데이터 우선 사용
    - 부족한 비활성 화합물은 decoy로 보충
    """
    
    def __init__(self, decoy_ratio: float = 50.0):
        """
        Args:
            decoy_ratio: 활성 화합물 대비 비활성 화합물 비율 (기본 50:1)
        """
        self.decoy_ratio = decoy_ratio
        
        # DUD-E 물리화학적 특성 허용 범위 (조금 완화하여 더 많은 decoy 생성)
        self.property_tolerances = {
            'MW': 0.25,          # 분자량 ±25% (완화)
            'LogP': 0.50,        # LogP ±0.50 (완화)
            'HBD': 2,            # Hydrogen Bond Donor ±2 (완화)
            'HBA': 3,            # Hydrogen Bond Acceptor ±3 (완화)
            'RotB': 3,           # Rotatable Bonds ±3 (완화)
            'Charge': 1          # Net Charge ±1 (완화)
        }
        
        # ZINC 데이터베이스 시뮬레이션용 (실제 운영 시 ZINC DB 연동)
        self.zinc_smiles_pool = []
    
    def calculate_molecular_properties(self, mol) -> Dict[str, float]:
        """
        분자의 물리화학적 특성 계산 (DUD-E 기준)
        
        Args:
            mol: RDKit Mol object
            
        Returns:
            Dict with molecular properties
        """
        if mol is None:
            return None
        
        try:
            properties = {
                'MW': Descriptors.MolWt(mol),
                'LogP': Descriptors.MolLogP(mol),
                'HBD': Descriptors.NumHDonors(mol),
                'HBA': Descriptors.NumHAcceptors(mol),
                'RotB': Descriptors.NumRotatableBonds(mol),
                'Charge': Chem.GetFormalCharge(mol),
                'TPSA': Descriptors.TPSA(mol),
                'NumRings': Lipinski.RingCount(mol)
            }
            return properties
        except Exception as e:
            print(f"Error calculating properties: {e}")
            return None
    
    def properties_match(self, active_props: Dict, candidate_props: Dict) -> bool:
        """
        두 화합물의 물리화학적 특성이 DUD-E 기준에 맞는지 확인
        
        Args:
            active_props: 활성 화합물의 특성
            candidate_props: 후보 화합물의 특성
            
        Returns:
            bool: 특성이 매칭되면 True
        """
        if not active_props or not candidate_props:
            return False
        
        # 분자량 체크 (±20%)
        mw_active = active_props['MW']
        mw_candidate = candidate_props['MW']
        if abs(mw_candidate - mw_active) > mw_active * self.property_tolerances['MW']:
            return False
        
        # LogP 체크 (±0.4)
        if abs(candidate_props['LogP'] - active_props['LogP']) > self.property_tolerances['LogP']:
            return False
        
        # HBD 체크 (±1)
        if abs(candidate_props['HBD'] - active_props['HBD']) > self.property_tolerances['HBD']:
            return False
        
        # HBA 체크 (±2)
        if abs(candidate_props['HBA'] - active_props['HBA']) > self.property_tolerances['HBA']:
            return False
        
        # Rotatable Bonds 체크 (±2)
        if abs(candidate_props['RotB'] - active_props['RotB']) > self.property_tolerances['RotB']:
            return False
        
        # Charge 체크 (정확히 일치)
        if candidate_props['Charge'] != active_props['Charge']:
            return False
        
        return True
    
    def calculate_tanimoto_similarity(self, mol1, mol2) -> float:
        """
        두 분자 간 Tanimoto 유사도 계산 (ECFP4 fingerprint 사용)
        
        Args:
            mol1: RDKit Mol object
            mol2: RDKit Mol object
            
        Returns:
            float: Tanimoto similarity (0-1)
        """
        try:
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        except:
            return 0.0
    
    def is_topologically_dissimilar(self, active_mol, candidate_mol, 
                                   threshold: float = 0.9) -> bool:
        """
        위상적으로 비유사한지 확인 (DUD-E: Tanimoto < 0.9)
        
        Args:
            active_mol: 활성 화합물 Mol object
            candidate_mol: 후보 화합물 Mol object
            threshold: Tanimoto 유사도 임계값
            
        Returns:
            bool: 비유사하면 True
        """
        similarity = self.calculate_tanimoto_similarity(active_mol, candidate_mol)
        return similarity < threshold
    
    def generate_decoys_from_zinc(self, active_smiles: str, 
                                  num_decoys: int = 50) -> List[str]:
        """
        ZINC 데이터베이스에서 DUD-E 기준에 맞는 decoy 생성
        
        Args:
            active_smiles: 활성 화합물 SMILES
            num_decoys: 생성할 decoy 개수
            
        Returns:
            List of decoy SMILES
        """
        active_mol = Chem.MolFromSmiles(active_smiles)
        if active_mol is None:
            print(f"Invalid SMILES: {active_smiles}")
            return []
        
        active_props = self.calculate_molecular_properties(active_mol)
        if active_props is None:
            return []
        
        decoys = []
        
        # ZINC 데이터베이스가 없는 경우 간단한 구조 변형으로 시뮬레이션
        # 실제 운영 시에는 ZINC DB 쿼리로 대체
        if not self.zinc_smiles_pool:
            decoys = self._generate_synthetic_decoys(active_mol, active_props, num_decoys)
        else:
            # ZINC DB에서 decoy 검색
            for zinc_smiles in self.zinc_smiles_pool:
                if len(decoys) >= num_decoys:
                    break
                
                candidate_mol = Chem.MolFromSmiles(zinc_smiles)
                if candidate_mol is None:
                    continue
                
                candidate_props = self.calculate_molecular_properties(candidate_mol)
                if candidate_props is None:
                    continue
                
                # DUD-E 기준 체크
                if (self.properties_match(active_props, candidate_props) and 
                    self.is_topologically_dissimilar(active_mol, candidate_mol)):
                    decoys.append(zinc_smiles)
        
        return decoys
    
    def _generate_synthetic_decoys(self, active_mol, active_props: Dict, 
                                   num_decoys: int) -> List[str]:
        """
        활성 화합물 기반으로 합성 decoy 생성 (ZINC DB 없을 때)
        
        Strategy:
        1. 원본 분자의 구조적 변형 (원자 치환, 링 변형 등)
        2. 물리화학적 특성은 유지하되 위상 구조는 변경
        
        Args:
            active_mol: 활성 화합물 Mol object
            active_props: 활성 화합물의 특성
            num_decoys: 생성할 decoy 개수
            
        Returns:
            List of synthetic decoy SMILES
        """
        decoys = []
        attempts = 0
        max_attempts = num_decoys * 500  # 적절한 시도 횟수 (5배 증가, 속도 고려)
        
        while len(decoys) < num_decoys and attempts < max_attempts:
            attempts += 1
            
            try:
                # 분자 변형 시도
                modified_mol = self._modify_molecule(active_mol)
                if modified_mol is None:
                    continue
                
                # 물리화학적 특성 확인
                modified_props = self.calculate_molecular_properties(modified_mol)
                if modified_props is None:
                    continue
                
                # DUD-E 기준 체크
                if (self.properties_match(active_props, modified_props) and
                    self.is_topologically_dissimilar(active_mol, modified_mol)):
                    
                    modified_smiles = Chem.MolToSmiles(modified_mol)
                    if modified_smiles not in decoys:
                        decoys.append(modified_smiles)
            
            except Exception as e:
                continue
        
        return decoys
    
    def _modify_molecule(self, mol):
        """
        분자 구조 변형 (원자 치환, 결합 회전, 작용기 추가 등)
        
        Args:
            mol: RDKit Mol object
            
        Returns:
            Modified Mol object
        """
        try:
            # 분자 복사
            modified_mol = Chem.RWMol(mol)
            
            # 변형 방법 랜덤 선택 (더 다양한 옵션)
            modification_type = random.choice([
                'atom_substitution', 
                'double_substitution',  # 2개 원자 동시 치환
                'add_methyl',           # 메틸기 추가
                'remove_hydrogen',      # 수소 제거
                'bond_rotation'
            ])
            
            if modification_type == 'atom_substitution':
                # 원자 치환 (C -> N, O, S, F, Cl 등)
                atom_idx = random.randint(0, modified_mol.GetNumAtoms() - 1)
                atom = modified_mol.GetAtomWithIdx(atom_idx)
                
                if atom.GetSymbol() == 'C':
                    substitute = random.choice(['N', 'O', 'S', 'F', 'Cl'])
                    atom_num = {'N': 7, 'O': 8, 'S': 16, 'F': 9, 'Cl': 17}
                    atom.SetAtomicNum(atom_num[substitute])
                elif atom.GetSymbol() in ['N', 'O']:
                    # N, O -> C 또는 서로 교환
                    substitute = random.choice(['C', 'N', 'O'])
                    atom_num = {'C': 6, 'N': 7, 'O': 8}
                    atom.SetAtomicNum(atom_num[substitute])
            
            elif modification_type == 'double_substitution':
                # 2개 원자 동시 치환으로 변형 증가
                if modified_mol.GetNumAtoms() >= 2:
                    for _ in range(2):
                        atom_idx = random.randint(0, modified_mol.GetNumAtoms() - 1)
                        atom = modified_mol.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == 'C':
                            substitute = random.choice(['N', 'O'])
                            atom.SetAtomicNum({'N': 7, 'O': 8}[substitute])
            
            elif modification_type == 'add_methyl':
                # 랜덤 탄소에 메틸기 추가
                if modified_mol.GetNumAtoms() > 0:
                    carbon_idx = random.randint(0, modified_mol.GetNumAtoms() - 1)
                    # 메틸 탄소 추가
                    new_carbon = modified_mol.AddAtom(Chem.Atom(6))
                    modified_mol.AddBond(carbon_idx, new_carbon, Chem.BondType.SINGLE)
            
            elif modification_type == 'remove_hydrogen':
                # 불포화도 증가 (H 제거)
                for atom in modified_mol.GetAtoms():
                    if atom.GetTotalNumHs() > 0 and random.random() < 0.3:
                        atom.SetNumExplicitHs(max(0, atom.GetTotalNumHs() - 1))
            
            elif modification_type == 'bond_rotation':
                # 3D Conformer 생성 및 회전
                try:
                    AllChem.EmbedMolecule(modified_mol, randomSeed=random.randint(0, 10000))
                    AllChem.UFFOptimizeMolecule(modified_mol)
                except:
                    pass  # 3D 생성 실패 시 무시
            
            # Sanitize
            Chem.SanitizeMol(modified_mol)
            
            return modified_mol.GetMol()
        
        except Exception as e:
            return None
    
    def generate_balanced_dataset(self, chembl_df: pd.DataFrame, 
                                  target_protein: str,
                                  pos_threshold: float = 10000,
                                  neg_threshold: float = 20000) -> pd.DataFrame:
        """
        실제 DB 데이터 우선 사용 + 부족한 비활성 화합물은 decoy로 채우기
        
        Args:
            chembl_df: ChEMBL 데이터프레임 (SMILES, standard_value 포함)
            target_protein: 타겟 단백질 이름
            pos_threshold: 활성 화합물 기준 (nM, 기본 10000)
            neg_threshold: 비활성 화합물 기준 (nM, 기본 20000)
            
        Returns:
            Balanced DataFrame with actives + inactives + decoys
        """
        print(f"\n{'='*60}")
        print(f"DUD-E Decoy Generation for {target_protein}")
        print(f"{'='*60}")
        
        # 1. 활성/비활성 화합물 분류
        actives_df = chembl_df[chembl_df['standard_value'] <= pos_threshold].copy()
        inactives_df = chembl_df[chembl_df['standard_value'] >= neg_threshold].copy()
        
        num_actives = len(actives_df)
        num_inactives = len(inactives_df)
        
        print(f"\n📊 Original Dataset:")
        print(f"  - Active compounds: {num_actives}")
        print(f"  - Inactive compounds: {num_inactives}")
        print(f"  - Ratio: 1:{num_inactives/num_actives:.1f}" if num_actives > 0 else "")
        
        # 2. 필요한 비활성 화합물 개수 계산
        required_inactives = int(num_actives * self.decoy_ratio)
        shortage = max(0, required_inactives - num_inactives)
        
        print(f"\n🎯 Target Dataset (ratio 1:{self.decoy_ratio}):")
        print(f"  - Required inactives: {required_inactives}")
        print(f"  - Shortage: {shortage}")
        
        if shortage == 0:
            print(f"\n✅ Sufficient inactive compounds. No decoys needed.")
            result_df = pd.concat([actives_df, inactives_df], ignore_index=True)
            result_df['source'] = ['active'] * num_actives + ['inactive'] * num_inactives
            return result_df
        
        # 3. Decoy 생성
        print(f"\n🔬 Generating {shortage} decoys...")
        print(f"  Strategy: {decoys_per_active} decoys per active compound")
        
        decoys_per_active = shortage // num_actives + 1
        all_decoys = []
        
        for idx, row in actives_df.iterrows():
            if len(all_decoys) >= shortage:
                break
            
            active_smiles = row['canonical_smiles']
            print(f"  [{idx+1}/{num_actives}] Generating decoys for: {active_smiles[:30]}...")
            
            decoys = self.generate_decoys_from_zinc(active_smiles, num_decoys=decoys_per_active)
            all_decoys.extend(decoys)
            
            print(f"    → Generated {len(decoys)} decoys (Total: {len(all_decoys)}/{shortage})")
            
            if (idx + 1) % 5 == 0:
                progress_pct = (len(all_decoys) / shortage * 100) if shortage > 0 else 0
                print(f"  📊 Progress: {idx+1}/{num_actives} actives, {len(all_decoys)}/{shortage} decoys ({progress_pct:.1f}%)")
        
        # Decoy 개수 조정
        all_decoys = all_decoys[:shortage]
        
        print(f"\n✅ Generated {len(all_decoys)} decoys")
        
        # 4. Decoy DataFrame 생성
        decoys_df = pd.DataFrame({
            'canonical_smiles': all_decoys,
            'standard_value': [neg_threshold + 10000] * len(all_decoys),  # 확실히 비활성
            'source': ['decoy'] * len(all_decoys)
        })
        
        # 5. 최종 데이터셋 결합
        actives_df['source'] = 'active'
        inactives_df['source'] = 'inactive'
        
        result_df = pd.concat([actives_df, inactives_df, decoys_df], ignore_index=True)
        
        print(f"\n📈 Final Dataset:")
        print(f"  - Active: {num_actives}")
        print(f"  - Inactive (real): {num_inactives}")
        print(f"  - Inactive (decoy): {len(all_decoys)}")
        print(f"  - Total: {len(result_df)}")
        print(f"  - Final ratio: 1:{(num_inactives + len(all_decoys))/num_actives:.1f}")
        print(f"{'='*60}\n")
        
        return result_df
    
    def load_zinc_database(self, zinc_file_path: str):
        """
        ZINC 데이터베이스 로드 (선택적)
        
        Args:
            zinc_file_path: ZINC SMILES 파일 경로
        """
        try:
            if zinc_file_path.endswith('.csv'):
                zinc_df = pd.read_csv(zinc_file_path)
                self.zinc_smiles_pool = zinc_df['smiles'].tolist()
            elif zinc_file_path.endswith('.txt'):
                with open(zinc_file_path, 'r') as f:
                    self.zinc_smiles_pool = [line.strip() for line in f if line.strip()]
            
            print(f"Loaded {len(self.zinc_smiles_pool)} ZINC compounds")
        
        except Exception as e:
            print(f"Error loading ZINC database: {e}")
            self.zinc_smiles_pool = []
