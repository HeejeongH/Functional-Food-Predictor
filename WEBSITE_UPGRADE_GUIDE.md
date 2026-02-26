# 웹사이트 업그레이드 가이드

## 🎨 현재 vs 업그레이드 옵션

### 현재 (Streamlit)
**장점**:
- ✅ 빠른 개발 (이미 완성!)
- ✅ Python만으로 구현
- ✅ 데이터 과학에 최적화

**단점**:
- ❌ 디자인 커스터마이징 제한적
- ❌ 단순한 UI
- ❌ 모바일 최적화 부족
- ❌ 전문적인 느낌 부족

---

## 🚀 업그레이드 옵션

### Option 1: **Gradio** (가장 빠른 업그레이드) ⭐⭐⭐⭐⭐

**난이도**: ⭐ (매우 쉬움)
**디자인**: ⭐⭐⭐⭐ (깔끔)
**추천**: **현재 코드 거의 그대로 사용!**

```python
# app/gradio_app.py
import gradio as gr
from modules.data_collector import collect_pci_data
from modules.feature_extractor import prepare_training_data

def predict_pci(target_gene, ic50_threshold, use_decoys, decoy_ratio):
    # 기존 함수 그대로 사용
    chembl_df = collect_pci_data([target_gene])
    prepared_df = prepare_training_data(
        chembl_df, target_gene, 
        use_decoys=use_decoys, 
        decoy_ratio=decoy_ratio
    )
    return f"Active: {len(prepared_df)}, Decoy: {(prepared_df['Y']==0).sum()}"

# Gradio 인터페이스 (매우 예쁨!)
with gr.Blocks(theme=gr.themes.Soft()) as demo:
    gr.Markdown("# 🧬 PCI Prediction Platform")
    
    with gr.Tab("Data Collection"):
        target_input = gr.Textbox(label="Target Gene", value="FTO")
        collect_btn = gr.Button("Collect Data", variant="primary")
        output = gr.Textbox(label="Result")
    
    with gr.Tab("Feature Extraction"):
        fp_type = gr.Dropdown(["ECFP4", "MACCS"], label="Fingerprint")
        decoy_check = gr.Checkbox(label="Use Decoys", value=True)
        decoy_ratio = gr.Slider(1, 100, value=50, label="Decoy Ratio")
    
    collect_btn.click(predict_pci, [target_input, ...], output)

demo.launch()
```

**장점**:
- ✅ Streamlit보다 훨씬 예쁨
- ✅ 기존 코드 거의 그대로 사용
- ✅ 모바일 반응형
- ✅ 테마 커스터마이징 쉬움
- ✅ HuggingFace에 무료 배포 가능

**단점**:
- ⚠️ Streamlit만큼 세밀한 제어는 어려움

---

### Option 2: **Flask + Vue.js** (프로페셔널) ⭐⭐⭐⭐⭐

**난이도**: ⭐⭐⭐ (중간)
**디자인**: ⭐⭐⭐⭐⭐ (완전 자유)
**추천**: **최고의 커스터마이징**

#### 백엔드 (Flask)
```python
# backend/app.py
from flask import Flask, request, jsonify
from flask_cors import CORS
from modules.data_collector import collect_pci_data

app = Flask(__name__)
CORS(app)

@app.route('/api/collect', methods=['POST'])
def collect_data():
    data = request.json
    target = data['target']
    
    chembl_df = collect_pci_data([target])
    
    return jsonify({
        'success': True,
        'count': len(chembl_df),
        'data': chembl_df.to_dict('records')[:100]
    })

@app.route('/api/train', methods=['POST'])
def train_model():
    # 모델 학습 로직
    pass

if __name__ == '__main__':
    app.run(debug=True, port=5000)
```

#### 프론트엔드 (Vue.js)
```vue
<!-- frontend/src/App.vue -->
<template>
  <div id="app">
    <nav class="navbar">
      <h1>🧬 PCI Prediction Platform</h1>
    </nav>
    
    <div class="container">
      <div class="card">
        <h2>Data Collection</h2>
        <input v-model="targetGene" placeholder="Target Gene" />
        <button @click="collectData" class="btn-primary">
          Collect Data
        </button>
      </div>
      
      <div class="results" v-if="results">
        <h3>Results: {{ results.count }} compounds</h3>
        <table>
          <tr v-for="item in results.data" :key="item.id">
            <td>{{ item.canonical_smiles }}</td>
            <td>{{ item.standard_value }}</td>
          </tr>
        </table>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  data() {
    return {
      targetGene: 'FTO',
      results: null
    }
  },
  methods: {
    async collectData() {
      const response = await fetch('http://localhost:5000/api/collect', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ target: this.targetGene })
      })
      this.results = await response.json()
    }
  }
}
</script>

<style>
.card {
  background: white;
  border-radius: 12px;
  padding: 2rem;
  box-shadow: 0 4px 6px rgba(0,0,0,0.1);
}
.btn-primary {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color: white;
  padding: 12px 24px;
  border: none;
  border-radius: 8px;
  cursor: pointer;
}
</style>
```

**장점**:
- ✅ 완전한 디자인 자유
- ✅ 프로페셔널한 느낌
- ✅ 최신 웹 기술
- ✅ 확장성 높음
- ✅ API 재사용 가능

**단점**:
- ❌ 개발 시간 많이 소요 (1-2주)
- ❌ 프론트엔드 지식 필요
- ❌ 배포 복잡

---

### Option 3: **FastAPI + React** (현대적) ⭐⭐⭐⭐⭐

**난이도**: ⭐⭐⭐ (중간)
**디자인**: ⭐⭐⭐⭐⭐ (완전 자유)
**추천**: **최신 기술 스택**

#### 백엔드 (FastAPI)
```python
# backend/main.py
from fastapi import FastAPI, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from modules.data_collector import collect_pci_data
from modules.model_trainer import ModelTrainer

app = FastAPI(title="PCI Prediction API")
app.add_middleware(CORSMiddleware, allow_origins=["*"])

class CollectRequest(BaseModel):
    target: str
    standard_type: str = 'IC50'

class TrainRequest(BaseModel):
    protein_name: str
    model_type: str
    use_decoys: bool = True
    decoy_ratio: float = 50.0

@app.post("/api/collect")
async def collect_data(request: CollectRequest):
    chembl_df = collect_pci_data([request.target], request.standard_type)
    return {
        "success": True,
        "count": len(chembl_df),
        "preview": chembl_df.head(10).to_dict('records')
    }

@app.post("/api/train")
async def train_model(request: TrainRequest, background_tasks: BackgroundTasks):
    # 백그라운드 작업으로 모델 학습
    background_tasks.add_task(train_model_task, request)
    return {"status": "training_started"}

@app.get("/")
async def root():
    return {"message": "PCI Prediction API"}
```

#### 프론트엔드 (React + TailwindCSS)
```jsx
// frontend/src/App.jsx
import React, { useState } from 'react';
import axios from 'axios';

function App() {
  const [targetGene, setTargetGene] = useState('FTO');
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);

  const collectData = async () => {
    setLoading(true);
    try {
      const response = await axios.post('http://localhost:8000/api/collect', {
        target: targetGene,
        standard_type: 'IC50'
      });
      setResults(response.data);
    } catch (error) {
      console.error(error);
    }
    setLoading(false);
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100">
      <nav className="bg-white shadow-lg">
        <div className="max-w-7xl mx-auto px-4 py-4">
          <h1 className="text-3xl font-bold text-indigo-600">
            🧬 PCI Prediction Platform
          </h1>
        </div>
      </nav>

      <div className="max-w-7xl mx-auto px-4 py-8">
        <div className="bg-white rounded-xl shadow-xl p-8">
          <h2 className="text-2xl font-semibold mb-6">Data Collection</h2>
          
          <div className="space-y-4">
            <input
              type="text"
              value={targetGene}
              onChange={(e) => setTargetGene(e.target.value)}
              placeholder="Enter target gene (e.g., FTO)"
              className="w-full px-4 py-3 border border-gray-300 rounded-lg focus:ring-2 focus:ring-indigo-500"
            />
            
            <button
              onClick={collectData}
              disabled={loading}
              className="w-full bg-gradient-to-r from-indigo-500 to-purple-600 text-white py-3 rounded-lg font-semibold hover:shadow-lg transition-all"
            >
              {loading ? 'Collecting...' : 'Collect Data'}
            </button>
          </div>

          {results && (
            <div className="mt-8">
              <h3 className="text-xl font-semibold mb-4">
                Results: {results.count} compounds
              </h3>
              <div className="overflow-x-auto">
                <table className="min-w-full divide-y divide-gray-200">
                  <thead className="bg-gray-50">
                    <tr>
                      <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">
                        SMILES
                      </th>
                      <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 uppercase">
                        IC50 (nM)
                      </th>
                    </tr>
                  </thead>
                  <tbody className="bg-white divide-y divide-gray-200">
                    {results.preview.map((item, idx) => (
                      <tr key={idx} className="hover:bg-gray-50">
                        <td className="px-6 py-4 text-sm font-mono">
                          {item.canonical_smiles}
                        </td>
                        <td className="px-6 py-4 text-sm">
                          {item.standard_value}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

export default App;
```

**장점**:
- ✅ 최신 기술 (FastAPI, React)
- ✅ 자동 API 문서 (Swagger)
- ✅ 타입 안정성 (Pydantic)
- ✅ 매우 빠른 성능
- ✅ 확장성 최고

**단점**:
- ❌ 학습 곡선
- ❌ 개발 시간 (1-2주)

---

### Option 4: **Next.js + Python Backend** (풀스택) ⭐⭐⭐⭐⭐

**난이도**: ⭐⭐⭐⭐ (어려움)
**디자인**: ⭐⭐⭐⭐⭐ (최고)
**추천**: **프로덕션 레벨**

```jsx
// frontend/pages/index.js
import { useState } from 'react';
import Head from 'next/head';

export default function Home() {
  return (
    <div className="container">
      <Head>
        <title>PCI Prediction Platform</title>
      </Head>

      <main className="main">
        <h1 className="title">
          🧬 PCI Prediction Platform
        </h1>
        
        <div className="grid">
          <div className="card">
            <h2>Data Collection &rarr;</h2>
            <p>Collect PCI data from ChEMBL and BindingDB</p>
          </div>

          <div className="card">
            <h2>Feature Extraction &rarr;</h2>
            <p>Generate fingerprints and molecular descriptors</p>
          </div>

          <div className="card">
            <h2>Model Training &rarr;</h2>
            <p>Train ML models with TabPFN, RandomForest, etc.</p>
          </div>

          <div className="card">
            <h2>SHAP Analysis &rarr;</h2>
            <p>Interpret model predictions and feature importance</p>
          </div>
        </div>
      </main>
    </div>
  );
}
```

**장점**:
- ✅ SEO 최적화
- ✅ 서버 사이드 렌더링
- ✅ 최고의 성능
- ✅ 프로덕션 준비 완료

**단점**:
- ❌ 가장 복잡
- ❌ 개발 시간 많이 소요 (2-3주)

---

## 🎯 추천 순서

### 1단계: Gradio로 빠른 업그레이드 (1-2일) ⭐⭐⭐⭐⭐

```bash
# 설치
pip install gradio

# 실행
python app/gradio_app.py
```

**이유**:
- 기존 코드 거의 그대로 사용
- Streamlit보다 훨씬 예쁨
- 빠른 개발

### 2단계: FastAPI + React (1-2주)

**더 전문적인 느낌 필요하면**

### 3단계: Next.js 풀스택 (2-3주)

**완전한 프로덕션 서비스**

---

## 💰 배포 옵션

### 무료 배포

| 플랫폼 | Gradio | Flask/FastAPI | Next.js |
|--------|--------|---------------|---------|
| **HuggingFace Spaces** | ✅ 무료 | ✅ 무료 | ❌ |
| **Vercel** | ❌ | ❌ | ✅ 무료 |
| **Railway** | ✅ 무료 | ✅ 무료 | ✅ 무료 |
| **Render** | ✅ 무료 | ✅ 무료 | ✅ 무료 |
| **Heroku** | ⚠️ 유료 | ⚠️ 유료 | ⚠️ 유료 |

### 프로페셔널 배포 (유료)

- **AWS**: EC2, Lambda, Amplify
- **Google Cloud**: Cloud Run, App Engine
- **Azure**: App Service
- **DigitalOcean**: Droplet, App Platform

---

## 🎨 디자인 예시

### Gradio Theme
```python
theme = gr.themes.Soft(
    primary_hue="indigo",
    secondary_hue="purple",
    neutral_hue="slate",
    font=[gr.themes.GoogleFont("Inter"), "sans-serif"]
)
```

### TailwindCSS (React)
```jsx
<div className="bg-gradient-to-r from-blue-500 to-purple-600 text-white p-8 rounded-2xl shadow-2xl">
  <h1 className="text-4xl font-bold mb-4">PCI Prediction</h1>
  <p className="text-lg opacity-90">Predict protein-compound interactions</p>
</div>
```

---

## 📊 비교 요약

| 항목 | Streamlit | Gradio | Flask+Vue | FastAPI+React |
|------|-----------|--------|-----------|---------------|
| **개발 시간** | 완료 ✅ | 1-2일 | 1-2주 | 1-2주 |
| **디자인 자유도** | ⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **학습 곡선** | 쉬움 | 쉬움 | 중간 | 중간 |
| **모바일 대응** | ⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **확장성** | ⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **API 제공** | ❌ | ⚠️ | ✅ | ✅ |
| **추천도** | 현재 | **즉시!** | 나중에 | 최종 |

---

## 🚀 즉시 시작하기

### Gradio로 업그레이드 (가장 추천!)

제가 바로 Gradio 버전을 만들어드릴까요?
- 기존 코드 재사용
- 1-2일이면 완성
- Streamlit보다 훨씬 예쁨
- HuggingFace에 무료 배포

**원하시면 지금 바로 Gradio 앱을 만들어드리겠습니다!** 🎨
