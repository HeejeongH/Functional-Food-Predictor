import React, { useState } from 'react'
import {
  Box,
  Typography,
  TextField,
  Button,
  Grid,
  Chip,
  CircularProgress,
  Alert,
  Card,
  CardContent,
  Divider,
  LinearProgress,
} from '@mui/material'
import {
  Add as AddIcon,
  Delete as DeleteIcon,
  CloudDownload as CloudDownloadIcon,
  ArrowForward as ArrowForwardIcon,
} from '@mui/icons-material'
import { useSnackbar } from 'notistack'
import { collectData } from '../../services/api'

export default function DataCollectionStep({ onNext, sessionId, setSessionId, updateWorkflowData }) {
  const { enqueueSnackbar } = useSnackbar()
  const [geneInput, setGeneInput] = useState('')
  const [genes, setGenes] = useState(['FTO'])
  const [ic50Threshold, setIc50Threshold] = useState(10000)
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)

  const handleAddGene = () => {
    if (geneInput.trim() && !genes.includes(geneInput.trim().toUpperCase())) {
      setGenes([...genes, geneInput.trim().toUpperCase()])
      setGeneInput('')
    }
  }

  const handleDeleteGene = (geneToDelete) => {
    setGenes(genes.filter((gene) => gene !== geneToDelete))
  }

  const handleCollect = async () => {
    if (genes.length === 0) {
      enqueueSnackbar('최소 1개의 유전자를 입력해주세요', { variant: 'warning' })
      return
    }

    setLoading(true)
    try {
      const data = await collectData({
        gene_names: genes,
        ic50_threshold: ic50Threshold,
      })
      
      setResult(data)
      setSessionId(data.session_id)
      updateWorkflowData({ collectionData: data })
      enqueueSnackbar(`✅ 데이터 수집 완료! 총 ${data.total_compounds}개 화합물`, { 
        variant: 'success' 
      })
    } catch (error) {
      console.error('Data collection error:', error)
      enqueueSnackbar(`❌ 데이터 수집 실패: ${error.message}`, { variant: 'error' })
    } finally {
      setLoading(false)
    }
  }

  return (
    <Box>
      <Typography variant="h4" gutterBottom sx={{ fontWeight: 700, mb: 3 }}>
        🧬 Step 1: 데이터 수집
      </Typography>
      <Typography variant="body1" color="text.secondary" sx={{ mb: 4 }}>
        ChEMBL과 BindingDB에서 타겟 유전자의 IC50 데이터를 자동으로 수집합니다.
      </Typography>

      <Grid container spacing={3}>
        {/* Input Section */}
        <Grid item xs={12} md={6}>
          <Card elevation={2} sx={{ height: '100%' }}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                타겟 유전자 설정
              </Typography>
              
              <Box sx={{ mb: 3 }}>
                <TextField
                  fullWidth
                  label="유전자 이름 입력"
                  placeholder="예: FTO, PPARG, EGFR"
                  value={geneInput}
                  onChange={(e) => setGeneInput(e.target.value.toUpperCase())}
                  onKeyPress={(e) => {
                    if (e.key === 'Enter') {
                      handleAddGene()
                    }
                  }}
                  InputProps={{
                    endAdornment: (
                      <Button 
                        onClick={handleAddGene}
                        startIcon={<AddIcon />}
                        disabled={!geneInput.trim()}
                      >
                        추가
                      </Button>
                    ),
                  }}
                />
              </Box>

              <Box sx={{ mb: 3 }}>
                <Typography variant="body2" color="text.secondary" gutterBottom>
                  선택된 유전자:
                </Typography>
                <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 1 }}>
                  {genes.map((gene) => (
                    <Chip
                      key={gene}
                      label={gene}
                      onDelete={() => handleDeleteGene(gene)}
                      color="primary"
                      deleteIcon={<DeleteIcon />}
                    />
                  ))}
                </Box>
              </Box>

              <Divider sx={{ my: 2 }} />

              <TextField
                fullWidth
                type="number"
                label="IC50 임계값 (nM)"
                value={ic50Threshold}
                onChange={(e) => setIc50Threshold(Number(e.target.value))}
                helperText="이 값 이하를 활성 화합물로 간주합니다"
                sx={{ mb: 3 }}
              />

              <Button
                fullWidth
                variant="contained"
                size="large"
                onClick={handleCollect}
                disabled={loading || genes.length === 0}
                startIcon={loading ? <CircularProgress size={20} /> : <CloudDownloadIcon />}
              >
                {loading ? '데이터 수집 중...' : '데이터 수집 시작'}
              </Button>
            </CardContent>
          </Card>
        </Grid>

        {/* Result Section */}
        <Grid item xs={12} md={6}>
          <Card elevation={2} sx={{ height: '100%' }}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                수집 결과
              </Typography>

              {loading && (
                <Box sx={{ textAlign: 'center', py: 4 }}>
                  <CircularProgress size={60} />
                  <Typography variant="body2" color="text.secondary" sx={{ mt: 2 }}>
                    ChEMBL 및 BindingDB에서 데이터를 가져오는 중...
                  </Typography>
                  <LinearProgress sx={{ mt: 2 }} />
                </Box>
              )}

              {!loading && !result && (
                <Alert severity="info" sx={{ mt: 2 }}>
                  타겟 유전자를 설정하고 '데이터 수집 시작' 버튼을 클릭하세요.
                </Alert>
              )}

              {result && (
                <Box>
                  <Grid container spacing={2} sx={{ mb: 3 }}>
                    <Grid item xs={6}>
                      <Card sx={{ bgcolor: 'primary.main', color: 'white', textAlign: 'center', p: 2 }}>
                        <Typography variant="h4" sx={{ fontWeight: 700 }}>
                          {result.total_compounds}
                        </Typography>
                        <Typography variant="body2">총 화합물</Typography>
                      </Card>
                    </Grid>
                    <Grid item xs={6}>
                      <Card sx={{ bgcolor: 'secondary.main', color: 'white', textAlign: 'center', p: 2 }}>
                        <Typography variant="h4" sx={{ fontWeight: 700 }}>
                          {result.unique_smiles}
                        </Typography>
                        <Typography variant="body2">고유 SMILES</Typography>
                      </Card>
                    </Grid>
                    <Grid item xs={6}>
                      <Card sx={{ bgcolor: 'success.main', color: 'white', textAlign: 'center', p: 2 }}>
                        <Typography variant="h4" sx={{ fontWeight: 700 }}>
                          {result.chembl_count}
                        </Typography>
                        <Typography variant="body2">ChEMBL</Typography>
                      </Card>
                    </Grid>
                    <Grid item xs={6}>
                      <Card sx={{ bgcolor: 'info.main', color: 'white', textAlign: 'center', p: 2 }}>
                        <Typography variant="h4" sx={{ fontWeight: 700 }}>
                          {result.bindingdb_count}
                        </Typography>
                        <Typography variant="body2">BindingDB</Typography>
                      </Card>
                    </Grid>
                  </Grid>

                  <Alert severity="success" sx={{ mb: 2 }}>
                    ✅ 데이터 수집이 완료되었습니다!
                  </Alert>

                  <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
                    <strong>IC50 통계:</strong>
                  </Typography>
                  <Typography variant="caption" display="block">
                    • 최소: {result.ic50_stats.min.toFixed(2)} nM
                  </Typography>
                  <Typography variant="caption" display="block">
                    • 최대: {result.ic50_stats.max.toFixed(2)} nM
                  </Typography>
                  <Typography variant="caption" display="block">
                    • 평균: {result.ic50_stats.mean.toFixed(2)} nM
                  </Typography>
                  <Typography variant="caption" display="block" sx={{ mb: 2 }}>
                    • 중앙값: {result.ic50_stats.median.toFixed(2)} nM
                  </Typography>

                  <Button
                    fullWidth
                    variant="contained"
                    color="primary"
                    size="large"
                    onClick={onNext}
                    endIcon={<ArrowForwardIcon />}
                  >
                    다음 단계로
                  </Button>
                </Box>
              )}
            </CardContent>
          </Card>
        </Grid>
      </Grid>
    </Box>
  )
}
