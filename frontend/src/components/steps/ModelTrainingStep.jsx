import React, { useState } from 'react'
import {
  Box,
  Typography,
  Button,
  Grid,
  Card,
  CardContent,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  CircularProgress,
  Alert,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
} from '@mui/material'
import { ArrowForward, ArrowBack, ModelTraining } from '@mui/icons-material'
import { useSnackbar } from 'notistack'
import { trainModel } from '../../services/api'

const modelTypes = ['TabPFN', 'RandomForest', 'GradientBoosting', 'SVM', 'LogisticRegression']
const featureSelectionMethods = ['mutual_info', 'rfe', 'pca', 'univariate']

export default function ModelTrainingStep({ onNext, onBack, sessionId, updateWorkflowData }) {
  const { enqueueSnackbar } = useSnackbar()
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  
  const [modelType, setModelType] = useState('TabPFN')
  const [featureSelection, setFeatureSelection] = useState('mutual_info')
  const [nFeatures, setNFeatures] = useState(100)

  const handleTrain = async () => {
    if (!sessionId) {
      enqueueSnackbar('먼저 데이터를 수집하고 특성을 변환해주세요', { variant: 'warning' })
      return
    }

    setLoading(true)
    try {
      const data = await trainModel({
        session_id: sessionId,
        model_type: modelType,
        feature_selection: featureSelection,
        n_features: nFeatures,
        test_size: 0.2,
      })
      
      setResult(data)
      updateWorkflowData({ modelData: data })
      enqueueSnackbar(`✅ 모델 학습 완료! F1: ${data.training_metrics.f1_score.toFixed(3)}`, { 
        variant: 'success' 
      })
    } catch (error) {
      enqueueSnackbar(`❌ 모델 학습 실패: ${error.message}`, { variant: 'error' })
    } finally {
      setLoading(false)
    }
  }

  return (
    <Box>
      <Typography variant="h4" gutterBottom sx={{ fontWeight: 700, mb: 3 }}>
        🤖 Step 3: 머신러닝 모델 학습
      </Typography>

      <Grid container spacing={3}>
        <Grid item xs={12} md={5}>
          <Card elevation={2}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                모델 설정
              </Typography>

              <FormControl fullWidth sx={{ mb: 2 }}>
                <InputLabel>ML 모델</InputLabel>
                <Select value={modelType} label="ML 모델" onChange={(e) => setModelType(e.target.value)}>
                  {modelTypes.map((type) => (
                    <MenuItem key={type} value={type}>{type}</MenuItem>
                  ))}
                </Select>
              </FormControl>

              <FormControl fullWidth sx={{ mb: 2 }}>
                <InputLabel>Feature Selection</InputLabel>
                <Select value={featureSelection} label="Feature Selection" onChange={(e) => setFeatureSelection(e.target.value)}>
                  {featureSelectionMethods.map((method) => (
                    <MenuItem key={method} value={method}>{method}</MenuItem>
                  ))}
                </Select>
              </FormControl>

              <Button
                fullWidth
                variant="contained"
                size="large"
                onClick={handleTrain}
                disabled={loading || !sessionId}
                startIcon={loading ? <CircularProgress size={20} /> : <ModelTraining />}
              >
                {loading ? '학습 중...' : '모델 학습 시작'}
              </Button>
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={7}>
          {result && (
            <Card elevation={2}>
              <CardContent>
                <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                  학습 결과
                </Typography>
                <Alert severity="success" sx={{ mb: 2 }}>
                  ✅ 모델: {result.model_type} | 특성: {result.n_selected_features}개
                </Alert>

                <TableContainer component={Paper} variant="outlined">
                  <Table size="small">
                    <TableHead>
                      <TableRow>
                        <TableCell><strong>지표</strong></TableCell>
                        <TableCell align="right"><strong>값</strong></TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      <TableRow>
                        <TableCell>Accuracy</TableCell>
                        <TableCell align="right">{result.training_metrics.accuracy.toFixed(4)}</TableCell>
                      </TableRow>
                      <TableRow>
                        <TableCell>Precision</TableCell>
                        <TableCell align="right">{result.training_metrics.precision.toFixed(4)}</TableCell>
                      </TableRow>
                      <TableRow>
                        <TableCell>Recall</TableCell>
                        <TableCell align="right">{result.training_metrics.recall.toFixed(4)}</TableCell>
                      </TableRow>
                      <TableRow>
                        <TableCell>F1 Score</TableCell>
                        <TableCell align="right"><strong>{result.training_metrics.f1_score.toFixed(4)}</strong></TableCell>
                      </TableRow>
                      <TableRow>
                        <TableCell>AUC</TableCell>
                        <TableCell align="right"><strong>{result.training_metrics.auc.toFixed(4)}</strong></TableCell>
                      </TableRow>
                    </TableBody>
                  </Table>
                </TableContainer>
              </CardContent>
            </Card>
          )}
        </Grid>
      </Grid>

      <Box sx={{ display: 'flex', justifyContent: 'space-between', mt: 4 }}>
        <Button variant="outlined" size="large" onClick={onBack} startIcon={<ArrowBack />}>
          이전
        </Button>
        <Button variant="contained" size="large" onClick={onNext} endIcon={<ArrowForward />} disabled={!result}>
          다음
        </Button>
      </Box>
    </Box>
  )
}
