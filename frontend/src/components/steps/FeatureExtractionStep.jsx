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
  TextField,
  Switch,
  FormControlLabel,
  Slider,
  CircularProgress,
  Alert,
  Chip,
} from '@mui/material'
import { ArrowForward, ArrowBack } from '@mui/icons-material'
import { useSnackbar } from 'notistack'
import { prepareFeatures } from '../../services/api'

const fingerprintTypes = ['ECFP4', 'ECFP6', 'MACCS', 'AtomPair', 'TopologicalTorsion', 'RDKit']
const fingerprintSizes = [256, 512, 1024, 2048]

export default function FeatureExtractionStep({ onNext, onBack, sessionId, updateWorkflowData }) {
  const { enqueueSnackbar } = useSnackbar()
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  
  const [fingerprintType, setFingerprintType] = useState('ECFP4')
  const [fingerprintSize, setFingerprintSize] = useState(1024)
  const [activeThreshold, setActiveThreshold] = useState(10000)
  const [inactiveThreshold, setInactiveThreshold] = useState(20000)
  const [useDecoys, setUseDecoys] = useState(true)
  const [decoyRatio, setDecoyRatio] = useState(50)
  const [decoyMethod, setDecoyMethod] = useState('dude')
  const [includeDescriptors, setIncludeDescriptors] = useState(true)

  const handlePrepare = async () => {
    if (!sessionId) {
      enqueueSnackbar('먼저 데이터를 수집해주세요', { variant: 'warning' })
      return
    }

    setLoading(true)
    try {
      const data = await prepareFeatures({
        session_id: sessionId,
        fingerprint_type: fingerprintType,
        fingerprint_size: fingerprintSize,
        active_threshold: activeThreshold,
        inactive_threshold: inactiveThreshold,
        use_decoys: useDecoys,
        decoy_ratio: decoyRatio,
        decoy_method: decoyMethod,
        include_descriptors: includeDescriptors,
      })
      
      setResult(data)
      updateWorkflowData({ featureData: data })
      enqueueSnackbar(`✅ 특성 변환 완료! ${data.total_samples}개 샘플, ${data.n_features}개 특성`, { 
        variant: 'success' 
      })
    } catch (error) {
      enqueueSnackbar(`❌ 특성 변환 실패: ${error.message}`, { variant: 'error' })
    } finally {
      setLoading(false)
    }
  }

  return (
    <Box>
      <Typography variant="h4" gutterBottom sx={{ fontWeight: 700, mb: 3 }}>
        🧬 Step 2: 분자 특성 변환
      </Typography>
      <Typography variant="body1" color="text.secondary" sx={{ mb: 4 }}>
        Fingerprint와 Molecular Descriptor로 화합물의 화학적 특성을 수치화합니다.
      </Typography>

      <Grid container spacing={3}>
        <Grid item xs={12} md={6}>
          <Card elevation={2}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                Fingerprint 설정
              </Typography>

              <FormControl fullWidth sx={{ mb: 2 }}>
                <InputLabel>Fingerprint 유형</InputLabel>
                <Select
                  value={fingerprintType}
                  label="Fingerprint 유형"
                  onChange={(e) => setFingerprintType(e.target.value)}
                >
                  {fingerprintTypes.map((type) => (
                    <MenuItem key={type} value={type}>{type}</MenuItem>
                  ))}
                </Select>
              </FormControl>

              <FormControl fullWidth sx={{ mb: 2 }}>
                <InputLabel>Fingerprint 크기</InputLabel>
                <Select
                  value={fingerprintSize}
                  label="Fingerprint 크기"
                  onChange={(e) => setFingerprintSize(e.target.value)}
                  disabled={fingerprintType === 'MACCS'}
                >
                  {fingerprintSizes.map((size) => (
                    <MenuItem key={size} value={size}>{size}</MenuItem>
                  ))}
                </Select>
              </FormControl>

              <FormControlLabel
                control={
                  <Switch
                    checked={includeDescriptors}
                    onChange={(e) => setIncludeDescriptors(e.target.checked)}
                  />
                }
                label="Molecular Descriptors 포함 (MW, LogP, TPSA 등)"
              />
            </CardContent>
          </Card>

          <Card elevation={2} sx={{ mt: 2 }}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                Decoy 생성 (DUDE 방식)
              </Typography>

              <FormControlLabel
                control={
                  <Switch
                    checked={useDecoys}
                    onChange={(e) => setUseDecoys(e.target.checked)}
                  />
                }
                label="Decoy (Negative Sample) 생성"
              />

              {useDecoys && (
                <>
                  <Typography variant="body2" gutterBottom sx={{ mt: 2 }}>
                    Decoy 비율: 1:{decoyRatio}
                  </Typography>
                  <Slider
                    value={decoyRatio}
                    onChange={(e, newValue) => setDecoyRatio(newValue)}
                    min={10}
                    max={100}
                    step={10}
                    marks
                    valueLabelDisplay="auto"
                    sx={{ mb: 2 }}
                  />

                  <FormControl fullWidth>
                    <InputLabel>Decoy 생성 방법</InputLabel>
                    <Select
                      value={decoyMethod}
                      label="Decoy 생성 방법"
                      onChange={(e) => setDecoyMethod(e.target.value)}
                    >
                      <MenuItem value="dude">DUDE-style (물리화학적 유사 + 구조적 상이)</MenuItem>
                      <MenuItem value="random">Random (빠른 생성)</MenuItem>
                    </Select>
                  </FormControl>
                </>
              )}
            </CardContent>
          </Card>
        </Grid>

        <Grid item xs={12} md={6}>
          <Card elevation={2}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                Activity 임계값 설정
              </Typography>

              <TextField
                fullWidth
                type="number"
                label="Active 임계값 (nM)"
                value={activeThreshold}
                onChange={(e) => setActiveThreshold(Number(e.target.value))}
                helperText="이 값 이하 = Active (Y=1)"
                sx={{ mb: 2 }}
              />

              <TextField
                fullWidth
                type="number"
                label="Inactive 임계값 (nM)"
                value={inactiveThreshold}
                onChange={(e) => setInactiveThreshold(Number(e.target.value))}
                helperText="이 값 이상 = Inactive (Y=0)"
                sx={{ mb: 3 }}
              />

              <Button
                fullWidth
                variant="contained"
                size="large"
                onClick={handlePrepare}
                disabled={loading || !sessionId}
                startIcon={loading ? <CircularProgress size={20} /> : null}
              >
                {loading ? '특성 변환 중...' : '특성 변환 시작'}
              </Button>

              {result && (
                <Box sx={{ mt: 3 }}>
                  <Alert severity="success" sx={{ mb: 2 }}>
                    ✅ 특성 변환 완료!
                  </Alert>
                  
                  <Grid container spacing={1}>
                    <Grid item xs={6}>
                      <Chip label={`활성: ${result.n_actives}`} color="success" sx={{ width: '100%' }} />
                    </Grid>
                    <Grid item xs={6}>
                      <Chip label={`비활성: ${result.n_inactives}`} color="error" sx={{ width: '100%' }} />
                    </Grid>
                    <Grid item xs={12}>
                      <Chip label={`총 특성: ${result.n_features}개`} color="primary" sx={{ width: '100%' }} />
                    </Grid>
                  </Grid>
                </Box>
              )}
            </CardContent>
          </Card>
        </Grid>
      </Grid>

      <Box sx={{ display: 'flex', justifyContent: 'space-between', mt: 4 }}>
        <Button
          variant="outlined"
          size="large"
          onClick={onBack}
          startIcon={<ArrowBack />}
        >
          이전
        </Button>
        <Button
          variant="contained"
          size="large"
          onClick={onNext}
          endIcon={<ArrowForward />}
          disabled={!result}
        >
          다음
        </Button>
      </Box>
    </Box>
  )
}
