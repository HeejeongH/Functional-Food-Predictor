"""
PCI Prediction Platform - Streamlit Web Application
단백질-화합물 상호작용 예측 통합 플랫폼
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import sys
from datetime import datetime
import json
import plotly.express as px
import plotly.graph_objects as go

# 모듈 경로 추가
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.data_collector import collect_pci_data, CHEMBL_AVAILABLE
from modules.feature_extractor import MolecularFeatureExtractor, prepare_training_data
from modules.model_trainer import ModelTrainer, TABPFN_AVAILABLE
from modules.shap_analyzer import analyze_model_with_shap
from modules.foodb_predictor import (
    prepare_foodb_features, 
    batch_predict_foodb,
    load_foodb_from_csv
)

# 페이지 설정
st.set_page_config(
    page_title="PCI Prediction Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 스타일
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #ff7f0e;
        margin-top: 2rem;
        margin-bottom: 1rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)

# 세션 상태 초기화
if 'collected_data' not in st.session_state:
    st.session_state.collected_data = None
if 'prepared_data' not in st.session_state:
    st.session_state.prepared_data = None
if 'trained_model' not in st.session_state:
    st.session_state.trained_model = None
if 'shap_results' not in st.session_state:
    st.session_state.shap_results = None
if 'foodb_predictions' not in st.session_state:
    st.session_state.foodb_predictions = None

# 메인 헤더
st.markdown('<div class="main-header">🧬 PCI Prediction Platform</div>', unsafe_allow_html=True)
st.markdown("**단백질-화합물 상호작용 예측 통합 플랫폼**")

# 사이드바
with st.sidebar:
    st.image("https://via.placeholder.com/300x100.png?text=PCI+Platform", use_column_width=True)
    st.markdown("---")
    
    st.markdown("### 📋 작업 단계")
    steps = [
        "1️⃣ 데이터 수집",
        "2️⃣ 특성 변환",
        "3️⃣ 모델 학습",
        "4️⃣ SHAP 분석",
        "5️⃣ FooDB 예측"
    ]
    for step in steps:
        st.markdown(f"- {step}")
    
    st.markdown("---")
    st.markdown("### ℹ️ 정보")
    st.info("""
    **개발자**: PCI Research Team
    
    **버전**: 1.0.0
    
    **라이브러리**:
    - ChEMBL API
    - RDKit
    - Scikit-learn
    - TabPFN
    - SHAP
    """)

# 메인 탭
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🔍 1. 데이터 수집", 
    "⚙️ 2. 특성 변환", 
    "🤖 3. 모델 학습",
    "📊 4. SHAP 분석",
    "🍎 5. FooDB 예측"
])

# ========================
# 탭 1: 데이터 수집
# ========================
with tab1:
    st.markdown('<div class="sub-header">🔍 ChEMBL 데이터 수집</div>', unsafe_allow_html=True)
    
    if not CHEMBL_AVAILABLE:
        st.error("⚠️ ChEMBL Web Resource Client가 설치되지 않았습니다. `pip install chembl_webresource_client`를 실행하세요.")
    else:
        col1, col2 = st.columns([2, 1])
        
        with col1:
            target_input = st.text_area(
                "타겟 유전자/단백질 이름 (쉼표로 구분)",
                value="FTO, Alpha-ketoglutarate-dependent dioxygenase",
                help="여러 개의 타겟을 쉼표로 구분하여 입력하세요"
            )
            
            standard_type = st.selectbox(
                "활성 데이터 타입",
                options=['IC50', 'EC50', 'Ki', 'Kd'],
                help="수집할 활성 데이터의 타입을 선택하세요"
            )
        
        with col2:
            st.markdown("#### 수집 옵션")
            protein_name = st.text_input(
                "단백질 이름",
                value="FTO",
                help="저장 파일명에 사용될 이름"
            )
        
        if st.button("📥 데이터 수집 시작", type="primary"):
            with st.spinner("ChEMBL에서 데이터를 수집하는 중..."):
                try:
                    target_list = [t.strip() for t in target_input.split(',')]
                    
                    # 데이터 수집
                    chembl_df = collect_pci_data(
                        chembl_search_list=target_list,
                        standard_type=standard_type
                    )
                    
                    st.session_state.collected_data = chembl_df
                    st.session_state.protein_name = protein_name
                    
                    st.success(f"✅ {len(chembl_df)}개의 데이터를 수집했습니다!")
                    
                except Exception as e:
                    st.error(f"❌ 오류 발생: {str(e)}")
    
    # 수집된 데이터 표시
    if st.session_state.collected_data is not None:
        st.markdown("---")
        st.markdown("### 📊 수집된 데이터")
        
        df = st.session_state.collected_data
        
        col1, col2, col3 = st.columns(3)
        col1.metric("총 데이터 수", len(df))
        col2.metric("Unique SMILES", df['canonical_smiles'].nunique())
        col3.metric("타겟 수", df['target_chembl_id'].nunique())
        
        # IC50 분포 시각화
        if 'standard_value' in df.columns:
            fig = px.histogram(
                df, 
                x=np.log10(df['standard_value'].astype(float)),
                title="IC50 분포 (log10 scale)",
                labels={'x': 'log10(IC50)', 'y': 'Count'}
            )
            st.plotly_chart(fig, use_container_width=True)
        
        # 데이터 테이블
        st.dataframe(df.head(100), use_container_width=True)
        
        # 다운로드 버튼
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 CSV 다운로드",
            data=csv,
            file_name=f"{protein_name}_chembl_data_{datetime.now().strftime('%Y%m%d')}.csv",
            mime="text/csv"
        )

# ========================
# 탭 2: 특성 변환
# ========================
with tab2:
    st.markdown('<div class="sub-header">⚙️ 분자 특성 변환</div>', unsafe_allow_html=True)
    
    if st.session_state.collected_data is None:
        st.warning("⚠️ 먼저 데이터를 수집해주세요.")
    else:
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("#### Fingerprint 설정")
            fp_type = st.selectbox(
                "Fingerprint 타입",
                options=list(MolecularFeatureExtractor.FINGERPRINT_TYPES.keys()),
                format_func=lambda x: MolecularFeatureExtractor.FINGERPRINT_TYPES[x]
            )
            
            fp_size = st.select_slider(
                "Fingerprint 크기",
                options=[256, 512, 1024, 2048],
                value=1024
            ) if fp_type != 'MACCS' else 167
            
            if fp_type == 'MACCS':
                st.info("MACCS Keys는 고정된 167개의 비트를 사용합니다.")
        
        with col2:
            st.markdown("#### 학습 데이터 설정")
            pos_threshold = st.number_input(
                "활성 임계값 (nM)",
                value=10000,
                min_value=1,
                help="이 값 이하는 활성(positive)으로 분류"
            )
            
            neg_threshold = st.number_input(
                "비활성 임계값 (nM)",
                value=20000,
                min_value=1,
                help="이 값 이상은 비활성(negative)으로 분류"
            )
            
            dataset_type = st.radio(
                "데이터셋 타입",
                options=['binary', 'active_only'],
                format_func=lambda x: "이진 분류 (활성 + 비활성)" if x == 'binary' else "활성 화합물만"
            )
            
            include_descriptors = st.checkbox("Molecular Descriptors 포함", value=True)
        
        # Decoy 생성 옵션 (이진 분류일 때만)
        if dataset_type == 'binary':
            st.markdown("---")
            st.markdown("#### 🎯 Decoy 생성 설정 (DUDE-style)")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                use_decoys = st.checkbox(
                    "Decoy 사용",
                    value=False,
                    help="DUDE-style decoy를 생성하여 비활성 샘플로 사용"
                )
            
            with col2:
                decoy_ratio = st.number_input(
                    "Decoy 비율 (1:N)",
                    value=50,
                    min_value=1,
                    max_value=100,
                    step=1,
                    help="활성 화합물 1개당 생성할 Decoy 개수",
                    disabled=not use_decoys
                )
            
            with col3:
                decoy_method = st.selectbox(
                    "Decoy 생성 방법",
                    options=['dude', 'random'],
                    format_func=lambda x: "DUDE-style (물리화학적 유사)" if x == 'dude' else "Random (빠름)",
                    disabled=not use_decoys
                )
            
            if use_decoys:
                st.info("""
                💡 **Decoy란?**
                - 물리화학적 특성은 활성 화합물과 유사하지만 구조적으로 다른 화합물
                - DUDE (Database of Useful Decoys: Enhanced) 방식 사용
                - 모델이 단순한 특성이 아닌 구조를 학습하도록 유도
                - 일반적으로 1:50 비율 권장 (활성 10개 → Decoy 500개)
                """)
        else:
            use_decoys = False
            decoy_ratio = 50
            decoy_method = 'dude'
        
        if st.button("🔄 특성 변환 시작", type="primary"):
            with st.spinner("분자 특성을 변환하는 중..."):
                try:
                    prepared_df = prepare_training_data(
                        chembl_df=st.session_state.collected_data,
                        protein_name=st.session_state.protein_name,
                        fp_type=fp_type,
                        fp_size=fp_size,
                        include_descriptors=include_descriptors,
                        pos_threshold=pos_threshold,
                        neg_threshold=neg_threshold,
                        dataset_type=dataset_type,
                        use_decoys=use_decoys,
                        decoy_ratio=float(decoy_ratio),
                        decoy_method=decoy_method,
                        decoy_source=st.session_state.collected_data if use_decoys else None
                    )
                    
                    st.session_state.prepared_data = prepared_df
                    st.session_state.fp_type = fp_type
                    st.session_state.fp_size = fp_size
                    
                    st.success(f"✅ {len(prepared_df)}개의 화합물 데이터가 준비되었습니다!")
                    
                except ValueError as ve:
                    st.error(f"❌ 데이터 오류: {str(ve)}")
                    st.info("""
                    💡 **문제 해결 방법**:
                    - IC50 임계값을 조정해보세요 (활성: 10000 → 50000, 비활성: 20000 → 10000)
                    - 수집된 데이터에 충분한 IC50 값이 있는지 확인하세요
                    - "데이터 수집" 탭에서 IC50 분포를 확인하세요
                    """)
                except Exception as e:
                    st.error(f"❌ 예상치 못한 오류 발생: {str(e)}")
                    with st.expander("🔍 디버깅 정보"):
                        st.code(str(e))
                        st.write("ChEMBL 데이터를 다시 수집하거나 다른 타겟을 시도해보세요.")
    
    # 준비된 데이터 표시
    if st.session_state.prepared_data is not None:
        st.markdown("---")
        st.markdown("### 📊 준비된 학습 데이터")
        
        df = st.session_state.prepared_data
        
        col1, col2, col3 = st.columns(3)
        col1.metric("총 샘플 수", len(df))
        col2.metric("피처 수", len([c for c in df.columns if c.startswith(('FP_', 'DESC_'))]))
        
        if 'Y' in df.columns:
            active_count = (df['Y'] == 1).sum()
            inactive_count = (df['Y'] == 0).sum()
            col3.metric("활성/비활성 비율", f"{active_count}/{inactive_count}")
            
            # 클래스 분포
            fig = px.pie(
                values=[active_count, inactive_count],
                names=['Active', 'Inactive'],
                title="클래스 분포"
            )
            st.plotly_chart(fig, use_container_width=True)
        
        st.dataframe(df.head(100), use_container_width=True)
        
        # 다운로드 버튼
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 CSV 다운로드",
            data=csv,
            file_name=f"{st.session_state.protein_name}_prepared_data_{datetime.now().strftime('%Y%m%d')}.csv",
            mime="text/csv"
        )

# ========================
# 탭 3: 모델 학습
# ========================
with tab3:
    st.markdown('<div class="sub-header">🤖 머신러닝 모델 학습</div>', unsafe_allow_html=True)
    
    if st.session_state.prepared_data is None:
        st.warning("⚠️ 먼저 데이터를 준비해주세요.")
    else:
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("#### 모델 설정")
            model_type = st.selectbox(
                "모델 타입",
                options=list(ModelTrainer.AVAILABLE_MODELS.keys()),
                format_func=lambda x: ModelTrainer.AVAILABLE_MODELS[x]
            )
            
            if model_type == 'TabPFN' and not TABPFN_AVAILABLE:
                st.error("⚠️ TabPFN이 설치되지 않았습니다.")
            
            test_size = st.slider(
                "테스트 데이터 비율",
                min_value=0.1,
                max_value=0.3,
                value=0.2,
                step=0.05
            )
        
        with col2:
            st.markdown("#### 피처 선택")
            feature_selection = st.selectbox(
                "피처 선택 방법",
                options=list(ModelTrainer.FEATURE_SELECTION_METHODS.keys()),
                format_func=lambda x: ModelTrainer.FEATURE_SELECTION_METHODS[x]
            )
            
            if feature_selection != 'none':
                n_features = st.select_slider(
                    "선택할 피처 수",
                    options=[50, 100, 200, 300, 400, 500],
                    value=200
                )
            else:
                n_features = None
        
        if st.button("🚀 모델 학습 시작", type="primary"):
            with st.spinner("모델을 학습하는 중..."):
                try:
                    trainer = ModelTrainer(
                        model_type=model_type,
                        n_features=n_features,
                        feature_selection_method=feature_selection
                    )
                    
                    metrics = trainer.train(
                        df=st.session_state.prepared_data,
                        test_size=test_size
                    )
                    
                    st.session_state.trained_model = trainer
                    st.session_state.metrics = metrics
                    
                    # 모델 저장
                    model_dir = "models"
                    os.makedirs(model_dir, exist_ok=True)
                    model_path = os.path.join(
                        model_dir, 
                        f"{st.session_state.protein_name}_{model_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.joblib"
                    )
                    trainer.save_model(model_path)
                    
                    st.success(f"✅ 모델 학습 완료! (저장: {model_path})")
                    
                except Exception as e:
                    st.error(f"❌ 오류 발생: {str(e)}")
    
    # 학습 결과 표시
    if st.session_state.trained_model is not None:
        st.markdown("---")
        st.markdown("### 📊 모델 성능")
        
        metrics = st.session_state.metrics
        
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Accuracy", f"{metrics['accuracy']:.4f}")
        col2.metric("F1 Score", f"{metrics['f1']:.4f}")
        col3.metric("AUC", f"{metrics['auc']:.4f}")
        col4.metric("Precision", f"{metrics['precision']:.4f}")
        
        # Confusion Matrix
        cm = np.array(metrics['confusion_matrix'])
        fig = px.imshow(
            cm,
            text_auto=True,
            labels=dict(x="Predicted", y="Actual", color="Count"),
            x=['Inactive', 'Active'],
            y=['Inactive', 'Active'],
            title="Confusion Matrix"
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # 메트릭 상세 정보
        with st.expander("📈 상세 메트릭"):
            st.json(metrics)

# ========================
# 탭 4: SHAP 분석
# ========================
with tab4:
    st.markdown('<div class="sub-header">📊 SHAP 분석</div>', unsafe_allow_html=True)
    
    if st.session_state.trained_model is None:
        st.warning("⚠️ 먼저 모델을 학습해주세요.")
    else:
        col1, col2 = st.columns(2)
        
        with col1:
            max_display = st.slider(
                "표시할 최대 피처 수",
                min_value=10,
                max_value=50,
                value=20
            )
        
        with col2:
            max_samples = st.slider(
                "SHAP 계산 샘플 수",
                min_value=50,
                max_value=500,
                value=100,
                step=50
            )
        
        if st.button("📊 SHAP 분석 시작", type="primary"):
            with st.spinner("SHAP 값을 계산하는 중..."):
                try:
                    df = st.session_state.prepared_data
                    feature_columns = [c for c in df.columns if c.startswith(('FP_', 'DESC_'))]
                    
                    X = df[feature_columns].values
                    y = df['Y'].values
                    
                    from sklearn.model_selection import train_test_split
                    X_train, X_test, y_train, y_test = train_test_split(
                        X, y, test_size=0.2, random_state=42, stratify=y
                    )
                    
                    # SHAP 분석
                    shap_results = analyze_model_with_shap(
                        model=st.session_state.trained_model.model,
                        X_train=X_train[:max_samples],
                        X_test=X_test[:max_samples],
                        feature_names=feature_columns,
                        save_dir=None,
                        max_display=max_display
                    )
                    
                    st.session_state.shap_results = shap_results
                    
                    st.success("✅ SHAP 분석 완료!")
                    
                except Exception as e:
                    st.error(f"❌ 오류 발생: {str(e)}")
        
        # SHAP 결과 표시
        if st.session_state.shap_results is not None:
            st.markdown("---")
            st.markdown("### 📊 SHAP 분석 결과")
            
            results = st.session_state.shap_results
            
            # 피처 중요도 표
            st.markdown("#### 상위 중요 피처")
            importance_df = results['importance_df'].head(max_display)
            
            fig = px.bar(
                importance_df,
                x='mean_abs_shap',
                y='feature',
                orientation='h',
                title="피처 중요도 (SHAP)"
            )
            fig.update_layout(yaxis={'categoryorder':'total ascending'})
            st.plotly_chart(fig, use_container_width=True)
            
            st.dataframe(importance_df, use_container_width=True)
            
            # SHAP 플롯
            col1, col2 = st.columns(2)
            with col1:
                if 'summary_plot' in results:
                    st.pyplot(results['summary_plot'])
            with col2:
                if 'bar_plot' in results:
                    st.pyplot(results['bar_plot'])
            
            # 다운로드
            csv = importance_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 피처 중요도 CSV 다운로드",
                data=csv,
                file_name=f"shap_importance_{datetime.now().strftime('%Y%m%d')}.csv",
                mime="text/csv"
            )

# ========================
# 탭 5: FooDB 예측
# ========================
with tab5:
    st.markdown('<div class="sub-header">🍎 FooDB 화합물 예측</div>', unsafe_allow_html=True)
    
    if st.session_state.trained_model is None:
        st.warning("⚠️ 먼저 모델을 학습해주세요.")
    else:
        st.markdown("#### FooDB 데이터 로드")
        
        uploaded_file = st.file_uploader(
            "FooDB CSV 파일 업로드",
            type=['csv'],
            help="FooDB 화합물 데이터 CSV 파일을 업로드하세요"
        )
        
        col1, col2 = st.columns(2)
        with col1:
            threshold = st.slider(
                "예측 임계값",
                min_value=0.0,
                max_value=1.0,
                value=0.5,
                step=0.05
            )
        
        with col2:
            batch_size = st.select_slider(
                "배치 크기",
                options=[100, 500, 1000, 2000],
                value=500
            )
        
        if uploaded_file is not None:
            if st.button("🔮 예측 시작", type="primary"):
                with st.spinner("FooDB 화합물을 예측하는 중..."):
                    try:
                        # FooDB 데이터 로드
                        foodb_df = pd.read_csv(uploaded_file)
                        st.info(f"📊 {len(foodb_df)}개의 화합물을 로드했습니다.")
                        
                        # 피처 추출
                        foodb_features = prepare_foodb_features(
                            foodb_df=foodb_df,
                            fp_type=st.session_state.fp_type,
                            fp_size=st.session_state.fp_size,
                            include_descriptors=True
                        )
                        
                        # 예측
                        feature_columns = st.session_state.trained_model.feature_names
                        
                        predictions_df = batch_predict_foodb(
                            foodb_df=foodb_features,
                            trainer=st.session_state.trained_model,
                            feature_columns=feature_columns,
                            batch_size=batch_size,
                            threshold=threshold
                        )
                        
                        st.session_state.foodb_predictions = predictions_df
                        
                        st.success(f"✅ {len(predictions_df)}개의 활성 화합물을 발견했습니다!")
                        
                    except Exception as e:
                        st.error(f"❌ 오류 발생: {str(e)}")
        
        # 예측 결과 표시
        if st.session_state.foodb_predictions is not None:
            st.markdown("---")
            st.markdown("### 📊 예측 결과")
            
            pred_df = st.session_state.foodb_predictions
            
            col1, col2 = st.columns(2)
            col1.metric("활성 화합물 수", len(pred_df))
            col2.metric("평균 확률", f"{pred_df['probability_class_1'].mean():.4f}")
            
            # 확률 분포
            fig = px.histogram(
                pred_df,
                x='probability_class_1',
                nbins=50,
                title="활성 확률 분포"
            )
            st.plotly_chart(fig, use_container_width=True)
            
            # 상위 화합물 표시
            st.markdown("#### 🏆 상위 활성 예측 화합물")
            display_columns = ['id', 'name', 'canonical_SMILES', 'probability_class_1']
            display_columns = [c for c in display_columns if c in pred_df.columns]
            st.dataframe(pred_df[display_columns].head(50), use_container_width=True)
            
            # 다운로드 버튼
            csv = pred_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 전체 예측 결과 CSV 다운로드",
                data=csv,
                file_name=f"foodb_predictions_{datetime.now().strftime('%Y%m%d')}.csv",
                mime="text/csv"
            )

# 푸터
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: gray;'>
    <p>PCI Prediction Platform v1.0.0 | © 2024 PCI Research Team</p>
    <p>Powered by Streamlit, RDKit, TabPFN, and SHAP</p>
</div>
""", unsafe_allow_html=True)
