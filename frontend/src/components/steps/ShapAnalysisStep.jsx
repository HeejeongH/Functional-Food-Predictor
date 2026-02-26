import React, { useState } from 'react'
import {
  Box,
  Typography,
  Button,
  Grid,
  Card,
  CardContent,
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
import { ArrowForward, ArrowBack, Analytics, Download } from '@mui/icons-material'
import { useSnackbar } from 'notistack'
import { runShapAnalysis, downloadResult } from '../../services/api'

export default function ShapAnalysisStep({ onNext, onBack, sessionId, updateWorkflowData }) {
  const { enqueueSnackbar } = useSnackbar()
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)

  const handleAnalyze = async () => {
    if (!sessionId) {
      enqueueSnackbar('먼저 모델을 학습해주세요', { variant: 'warning' })
      return
    }

    setLoading(true)
    try {
      const data = await runShapAnalysis({ session_id: sessionId })
      setResult(data)
      updateWorkflowData({ shapData: data })
      enqueueSnackbar(`✅ SHAP 분석 완료!`, { variant: 'success' })
    } catch (error) {
      enqueueSnackbar(`❌ SHAP 분석 실패: ${error.message}`, { variant: 'error' })
    } finally {
      setLoading(false)
    }
  }

  return (
    <Box>
      <Typography variant="h4" gutterBottom sx={{ fontWeight: 700, mb: 3 }}>
        📊 Step 4: SHAP 분석
      </Typography>
      <Typography variant="body1" color="text.secondary" sx={{ mb: 4 }}>
        SHAP을 사용하여 모델의 예측에 중요한 화학적 특성을 도출합니다.
      </Typography>

      <Grid container spacing={3}>
        <Grid item xs={12} md={6}>
          <Card elevation={2}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                SHAP 분석 실행
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
                SHAP (SHapley Additive exPlanations)는 머신러닝 모델의 예측을 설명하는 강력한 도구입니다.
              </Typography>

              <Button
                fullWidth
                variant="contained"
                size="large"
                onClick={handleAnalyze}
                disabled={loading || !sessionId}
                startIcon={loading ? <CircularProgress size={20} /> : <Analytics />}
              >
                {loading ? 'SHAP 분석 중...' : 'SHAP 분석 시작'}
              </Button>

              {result && (
                <Box sx={{ mt: 3 }}>
                  <Alert severity="success">✅ SHAP 분석 완료!</Alert>
                  <Box sx={{ mt: 2 }}>
                    <Button
                      fullWidth
                      variant="outlined"
                      startIcon={<Download />}
                      href={downloadResult(result.csv_download.split('/').pop())}
                      sx={{ mb: 1 }}
                    >
                      Feature Importance CSV 다운로드
                    </Button>
                  </Box>
                </Box>
              )}
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={6}>
          {result && (
            <Card elevation={2}>
              <CardContent>
                <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                  Top 10 중요 특성
                </Typography>
                <TableContainer component={Paper} variant="outlined">
                  <Table size="small">
                    <TableHead>
                      <TableRow>
                        <TableCell><strong>순위</strong></TableCell>
                        <TableCell><strong>특성</strong></TableCell>
                        <TableCell align="right"><strong>중요도</strong></TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {result.top_features.slice(0, 10).map((item, index) => (
                        <TableRow key={index}>
                          <TableCell>{index + 1}</TableCell>
                          <TableCell>{item.feature}</TableCell>
                          <TableCell align="right">{item.importance.toFixed(4)}</TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </TableContainer>

                <Box sx={{ mt: 2 }}>
                  {result.summary_plot && (
                    <img
                      src={result.summary_plot}
                      alt="SHAP Summary Plot"
                      style={{ width: '100%', borderRadius: 8 }}
                    />
                  )}
                </Box>
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
