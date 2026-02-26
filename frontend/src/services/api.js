import axios from 'axios'

const API_BASE_URL = '/api'

const api = axios.create({
  baseURL: API_BASE_URL,
  headers: {
    'Content-Type': 'application/json',
  },
  timeout: 600000, // 10 minutes for long-running tasks
})

// Health check
export const checkHealth = async () => {
  const response = await api.get('/health')
  return response.data
}

// Step 1: Data Collection
export const collectData = async (params) => {
  const response = await api.post('/collect-data', params)
  return response.data
}

// Step 2: Feature Preparation
export const prepareFeatures = async (params) => {
  const response = await api.post('/prepare-features', params)
  return response.data
}

// Step 3: Model Training
export const trainModel = async (params) => {
  const response = await api.post('/train-model', params)
  return response.data
}

// Step 4: SHAP Analysis
export const runShapAnalysis = async (params) => {
  const response = await api.post('/shap-analysis', params)
  return response.data
}

// Step 5: FooDB Prediction
export const predictFooDB = async (params) => {
  const response = await api.post('/predict-foodb', params)
  return response.data
}

// Download results
export const downloadResult = (filename) => {
  return `${API_BASE_URL}/results/${filename}`
}

// List sessions
export const listSessions = async () => {
  const response = await api.get('/sessions')
  return response.data
}

export default api
