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
  Slider,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
  Chip,
} from '@mui/material'
import { ArrowBack, Refresh, Download, Restaurant } from '@mui/icons-material'
import { useSnackbar } from 'notistack'
import { predictFooDB, downloadResult } from '../../services/api'

export default function FooDBPredictionStep({ onBack, onReset, sessionId, updateWorkflowData }) {
  const { enqueueSnackbar } = useSnackbar()
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [threshold, setThreshold] = useState(0.5)

  const handlePredict = async () => {
    if (!sessionId) {
      enqueueSnackbar('먼저 모든 이전 단계를 완료해주세요', { variant: 'warning' })
      return
    }

    setLoading(true)
    try {
      const data = await predictFooDB({
        session_id: sessionId,
        prediction_threshold: threshold,
      })
      setResult(data)
      updateWorkflowData({ foodbData: data })
      enqueueSnackbar(`✅ FooDB 예측 완료! ${data.active_compounds}개 활성 화합물 발견`, { 
        variant: 'success' 
      })
    } catch (error) {
      enqueueSnackbar(`❌ FooDB 예측 실패: ${error.message}`, { variant: 'error' })
    } finally {
      setLoading(false)
    }
  }

  return (
    <Box>
      <Typography variant="h4" gutterBottom sx={{ fontWeight: 700, mb: 3 }}>
        🍎 Step 5: FooDB 실제 데이터 예측
      </Typography>
      <Typography variant="body1" color="text.secondary" sx={{ mb: 4 }}>
        70,000개 이상의 식품 화합물 데이터베이스에서 활성 화합물을 예측합니다.
      </Typography>

      <Grid container spacing={3}>
        <Grid item xs={12} md={5}>
          <Card elevation={2}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                예측 설정
              </Typography>

              <Typography variant="body2" gutterBottom>
                활성 임계값: {threshold.toFixed(2)}
              </Typography>
              <Slider
                value={threshold}
                onChange={(e, newValue) => setThreshold(newValue)}
                min={0.1}
                max={0.9}
                step={0.05}
                marks
                valueLabelDisplay="auto"
                sx={{ mb: 3 }}
              />

              <Button
                fullWidth
                variant="contained"
                size="large"
                onClick={handlePredict}
                disabled={loading || !sessionId}
                startIcon={loading ? <CircularProgress size={20} /> : <Restaurant />}
              >
                {loading ? 'FooDB 예측 중...' : 'FooDB 예측 시작'}
              </Button>

              {result && (
                <Box sx={{ mt: 3 }}>
                  <Alert severity="success" sx={{ mb: 2 }}>
                    ✅ 예측 완료!
                  </Alert>
                  <Grid container spacing={1}>
                    <Grid item xs={12}>
                      <Chip
                        label={`전체: ${result.total_compounds}개`}
                        color="primary"
                        sx={{ width: '100%', fontSize: '1rem', py: 2 }}
                      />
                    </Grid>
                    <Grid item xs={12}>
                      <Chip
                        label={`활성: ${result.active_compounds}개 (${result.active_percentage.toFixed(2)}%)`}
                        color="success"
                        sx={{ width: '100%', fontSize: '1rem', py: 2 }}
                      />
                    </Grid>
                  </Grid>

                  <Box sx={{ mt: 2 }}>
                    <Button
                      fullWidth
                      variant="outlined"
                      startIcon={<Download />}
                      href={downloadResult(result.active_csv_download.split('/').pop())}
                      sx={{ mb: 1 }}
                    >
                      활성 화합물 CSV 다운로드
                    </Button>
                    <Button
                      fullWidth
                      variant="outlined"
                      startIcon={<Download />}
                      href={downloadResult(result.csv_download.split('/').pop())}
                    >
                      전체 예측 결과 CSV 다운로드
                    </Button>
                  </Box>
                </Box>
              )}
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={7}>
          {result && (
            <Card elevation={2}>
              <CardContent>
                <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                  Top 10 활성 화합물
                </Typography>
                <TableContainer component={Paper} variant="outlined">
                  <Table size="small">
                    <TableHead>
                      <TableRow>
                        <TableCell><strong>순위</strong></TableCell>
                        <TableCell><strong>화합물명</strong></TableCell>
                        <TableCell align="right"><strong>확률</strong></TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {result.top_predictions.slice(0, 10).map((item, index) => (
                        <TableRow key={index}>
                          <TableCell>{index + 1}</TableCell>
                          <TableCell>{item.name || item.public_id}</TableCell>
                          <TableCell align="right">
                            <Chip
                              label={item.prediction_proba.toFixed(3)}
                              color="success"
                              size="small"
                            />
                          </TableCell>
                        </TableRow>
                      ))}
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
        <Button variant="contained" size="large" onClick={onReset} startIcon={<Refresh />} color="secondary">
          새로 시작
        </Button>
      </Box>
    </Box>
  )
}
