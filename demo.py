"""
PCI Platform Demo Script
간단한 기능 테스트 및 데모
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
import numpy as np
from modules.feature_extractor import MolecularFeatureExtractor

def test_feature_extractor():
    """분자 특성 추출기 테스트"""
    print("=" * 60)
    print("Testing Molecular Feature Extractor")
    print("=" * 60)
    
    # 테스트 SMILES
    test_smiles = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
    ]
    
    names = ["Aspirin", "Ibuprofen", "Caffeine"]
    
    # 다양한 fingerprint 타입 테스트
    fp_types = ['ECFP4', 'MACCS', 'RDKit']
    
    for fp_type in fp_types:
        print(f"\n{fp_type} Fingerprint:")
        print("-" * 60)
        
        extractor = MolecularFeatureExtractor(fp_type=fp_type)
        
        for name, smiles in zip(names, test_smiles):
            fp = extractor.smiles_to_fingerprint(smiles)
            print(f"{name:15s} - FP Size: {len(fp)}, Non-zero: {int(fp.sum())}")
            
            # Descriptor 계산
            desc = extractor.calculate_descriptors(smiles)
            print(f"{'':15s}   MW: {desc.get('MolWt', 0):.2f}, LogP: {desc.get('LogP', 0):.2f}")
    
    print("\n✅ Feature Extractor Test Completed!")

def test_model_types():
    """사용 가능한 모델 타입 확인"""
    print("\n" + "=" * 60)
    print("Available Model Types")
    print("=" * 60)
    
    from modules.model_trainer import ModelTrainer
    
    for key, value in ModelTrainer.AVAILABLE_MODELS.items():
        print(f"  - {key:20s}: {value}")
    
    print("\n" + "=" * 60)
    print("Available Feature Selection Methods")
    print("=" * 60)
    
    for key, value in ModelTrainer.FEATURE_SELECTION_METHODS.items():
        print(f"  - {key:20s}: {value}")

def create_sample_dataset():
    """샘플 데이터셋 생성"""
    print("\n" + "=" * 60)
    print("Creating Sample Dataset with Decoys")
    print("=" * 60)
    
    # 샘플 화합물 (활성/비활성)
    active_compounds = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    ]
    
    inactive_compounds = [
        "C",  # Methane
        "CC",  # Ethane
        "O",  # Water
    ]
    
    # DataFrame 생성
    data = []
    for smiles in active_compounds:
        data.append({'canonical_smiles': smiles, 'standard_value': 100, 'standard_relation': '='})
    for smiles in inactive_compounds:
        data.append({'canonical_smiles': smiles, 'standard_value': 50000, 'standard_relation': '='})
    
    df = pd.DataFrame(data)
    
    print(f"\nCreated dataset with {len(df)} compounds")
    print(f"  - Active: {len(active_compounds)}")
    print(f"  - Inactive: {len(inactive_compounds)}")
    
    # 특성 추출 테스트
    from modules.feature_extractor import prepare_training_data
    
    try:
        prepared_df = prepare_training_data(
            chembl_df=df,
            protein_name="TestProtein",
            fp_type='ECFP4',
            fp_size=512,
            include_descriptors=True,
            pos_threshold=1000,
            neg_threshold=10000,
            dataset_type='binary'
        )
        
        print(f"\n✅ Prepared training data:")
        print(f"  - Total samples: {len(prepared_df)}")
        print(f"  - Total features: {len([c for c in prepared_df.columns if c.startswith(('FP_', 'DESC_'))])}")
        print(f"  - Active: {(prepared_df['Y'] == 1).sum()}")
        print(f"  - Inactive: {(prepared_df['Y'] == 0).sum()}")
        
        return prepared_df
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return None

def main():
    """메인 데모 함수"""
    print("\n" + "=" * 60)
    print("PCI Prediction Platform - Demo Script")
    print("=" * 60)
    
    # 1. Feature Extractor 테스트
    test_feature_extractor()
    
    # 2. 모델 타입 확인
    test_model_types()
    
    # 3. 샘플 데이터셋 생성
    sample_df = create_sample_dataset()
    
    print("\n" + "=" * 60)
    print("Demo Completed Successfully! 🎉")
    print("=" * 60)
    print("\nNext Steps:")
    print("  1. Run the web application: ./run.sh")
    print("  2. Upload your own data")
    print("  3. Train models and analyze results")
    print("=" * 60 + "\n")

if __name__ == "__main__":
    main()
