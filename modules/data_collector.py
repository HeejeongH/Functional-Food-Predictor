"""
데이터 수집 모듈 - ChEMBL과 BindingDB에서 PCI 데이터 수집
"""

import requests
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from rdkit import Chem

try:
    from chembl_webresource_client.new_client import new_client
    CHEMBL_AVAILABLE = True
except ImportError:
    CHEMBL_AVAILABLE = False
    print("Warning: chembl_webresource_client not installed")


class ChEMBLTargetCollector:
    """ChEMBL 데이터베이스에서 타겟 단백질 정보 수집"""
    
    def __init__(self):
        if not CHEMBL_AVAILABLE:
            raise ImportError("chembl_webresource_client is required for ChEMBL data collection")
        self.target = new_client.target
        self.activity = new_client.activity

    def get_target_data(self, target_list: List[str], standard_type: str = 'IC50') -> pd.DataFrame:
        """
        타겟 리스트에서 활성 데이터를 수집
        
        Args:
            target_list: 검색할 타겟 이름 리스트
            standard_type: 활성 타입 (IC50, EC50, Ki 등)
        
        Returns:
            수집된 데이터 DataFrame
        """
        all_targets = self._filter_targets(target_list)
        results = self._collect_activities(all_targets, standard_type)
        df = pd.DataFrame(results)
        
        # Uniprot 매핑 추가
        uniprot_mapping = self._get_uniprot_mapping()
        df['Standard_Name'] = df['target_pref_name'].apply(
            lambda x: self._map_with_uniprot(x, uniprot_mapping)
        )
        
        return df

    def _filter_targets(self, target_list: List[str]) -> Dict:
        """타겟 필터링"""
        all_targets = {}
        for target_i in target_list:
            filtered_target = self.target.filter(
                target_synonym__icontains=target_i, 
                organism='Homo sapiens'
            )
            for t in filtered_target:
                all_targets[t['target_chembl_id']] = t
        print(f"Found {len(all_targets)} targets")
        return all_targets

    def _collect_activities(self, all_targets: Dict, standard_type: str) -> List[Dict]:
        """활성 데이터 수집"""
        results = []
        for target_id, t in all_targets.items():
            valid_data = self.activity.filter(
                target_chembl_id=target_id,
                standard_type=standard_type
            )
            for act in valid_data:
                if act['standard_value']:
                    results.append({
                        'target_pref_name': act['target_pref_name'],
                        'target_chembl_id': act['target_chembl_id'],
                        'molecule_pref_name': act['molecule_pref_name'],
                        'molecule_chembl_id': act['molecule_chembl_id'],
                        'canonical_smiles': act['canonical_smiles'],
                        'standard_type': act['standard_type'],
                        'standard_value': act['standard_value'],
                        'standard_units': act['standard_units'],
                        'standard_relation': act['standard_relation']
                    })
        return results
    
    def _get_uniprot_mapping(self) -> Dict[str, str]:
        """Uniprot에서 유전자 심볼 매핑 가져오기"""
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': 'organism_id:9606',
            'format': 'json',
            'fields': 'accession,gene_names,protein_name',
            'size': 500
        }
        
        try:
            response = requests.get(url, params=params, timeout=10)
            data = response.json()
            
            name_to_symbol = {}
            for entry in data.get('results', []):
                gene_symbol = ""
                if 'genes' in entry and entry['genes']:
                    gene_symbol = entry['genes'][0].get('geneName', {}).get('value', '')
                
                if gene_symbol:
                    if 'proteinDescription' in entry:
                        rec_name = entry['proteinDescription'].get('recommendedName', {})
                        if rec_name:
                            full_name = rec_name.get('fullName', {}).get('value', '')
                            if full_name:
                                name_to_symbol[full_name] = gene_symbol
                        
                        alt_names = entry['proteinDescription'].get('alternativeNames', [])
                        for alt in alt_names:
                            alt_name = alt.get('fullName', {}).get('value', '')
                            if alt_name:
                                name_to_symbol[alt_name] = gene_symbol
            
            return name_to_symbol
        except Exception as e:
            print(f"Warning: Failed to get Uniprot mapping: {e}")
            return {}
    
    def _map_with_uniprot(self, target_name: str, mapping_dict: Dict[str, str]) -> str:
        """타겟 이름을 유전자 심볼로 매핑"""
        return mapping_dict.get(target_name, target_name)


def to_canonical_smiles(smiles: str) -> Optional[str]:
    """SMILES를 canonical 형태로 변환"""
    if pd.isna(smiles):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None


def collect_pci_data(
    chembl_search_list: List[str],
    standard_type: str = 'IC50',
    save_path: Optional[str] = None
) -> tuple:
    """
    ChEMBL에서 PCI 데이터 수집
    
    Args:
        chembl_search_list: ChEMBL 검색 키워드 리스트
        standard_type: 활성 타입
        save_path: 저장 경로 (선택사항)
    
    Returns:
        (chembl_df,) - ChEMBL 데이터프레임
    """
    print("Collecting ChEMBL data...")
    collector = ChEMBLTargetCollector()
    chembl_df = collector.get_target_data(chembl_search_list, standard_type).dropna()
    print(f"ChEMBL data collected: {len(chembl_df)} records")
    
    if save_path:
        chembl_df.to_csv(save_path, index=False)
        print(f"Data saved to {save_path}")
    
    return chembl_df
