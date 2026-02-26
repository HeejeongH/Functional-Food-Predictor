import React, { useState } from 'react'
import {
  Box,
  Container,
  AppBar,
  Toolbar,
  Typography,
  Stepper,
  Step,
  StepLabel,
  Paper,
  Button,
  Grid,
  useTheme,
} from '@mui/material'
import {
  Science as ScienceIcon,
  BiotechOutlined as BiotechIcon,
} from '@mui/icons-material'
import { useSnackbar } from 'notistack'

// Step Components
import DataCollectionStep from './components/steps/DataCollectionStep'
import FeatureExtractionStep from './components/steps/FeatureExtractionStep'
import ModelTrainingStep from './components/steps/ModelTrainingStep'
import ShapAnalysisStep from './components/steps/ShapAnalysisStep'
import FooDBPredictionStep from './components/steps/FooDBPredictionStep'

const steps = [
  '데이터 수집',
  '특성 변환',
  '모델 학습',
  'SHAP 분석',
  'FooDB 예측'
]

function App() {
  const theme = useTheme()
  const { enqueueSnackbar } = useSnackbar()
  const [activeStep, setActiveStep] = useState(0)
  const [sessionId, setSessionId] = useState(null)
  const [workflowData, setWorkflowData] = useState({})

  const handleNext = () => {
    setActiveStep((prevActiveStep) => prevActiveStep + 1)
  }

  const handleBack = () => {
    setActiveStep((prevActiveStep) => prevActiveStep - 1)
  }

  const handleReset = () => {
    setActiveStep(0)
    setSessionId(null)
    setWorkflowData({})
    enqueueSnackbar('워크플로우가 초기화되었습니다', { variant: 'info' })
  }

  const updateWorkflowData = (stepData) => {
    setWorkflowData((prev) => ({ ...prev, ...stepData }))
  }

  const getStepContent = (step) => {
    switch (step) {
      case 0:
        return (
          <DataCollectionStep
            onNext={handleNext}
            sessionId={sessionId}
            setSessionId={setSessionId}
            updateWorkflowData={updateWorkflowData}
          />
        )
      case 1:
        return (
          <FeatureExtractionStep
            onNext={handleNext}
            onBack={handleBack}
            sessionId={sessionId}
            workflowData={workflowData}
            updateWorkflowData={updateWorkflowData}
          />
        )
      case 2:
        return (
          <ModelTrainingStep
            onNext={handleNext}
            onBack={handleBack}
            sessionId={sessionId}
            workflowData={workflowData}
            updateWorkflowData={updateWorkflowData}
          />
        )
      case 3:
        return (
          <ShapAnalysisStep
            onNext={handleNext}
            onBack={handleBack}
            sessionId={sessionId}
            workflowData={workflowData}
            updateWorkflowData={updateWorkflowData}
          />
        )
      case 4:
        return (
          <FooDBPredictionStep
            onBack={handleBack}
            onReset={handleReset}
            sessionId={sessionId}
            workflowData={workflowData}
            updateWorkflowData={updateWorkflowData}
          />
        )
      default:
        return 'Unknown step'
    }
  }

  return (
    <Box sx={{ minHeight: '100vh', background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)' }}>
      {/* Header */}
      <AppBar 
        position="static" 
        elevation={0}
        sx={{ 
          background: 'rgba(255, 255, 255, 0.95)',
          backdropFilter: 'blur(10px)',
          borderBottom: '1px solid rgba(0, 0, 0, 0.05)'
        }}
      >
        <Toolbar sx={{ py: 2 }}>
          <Box sx={{ display: 'flex', alignItems: 'center', flexGrow: 1 }}>
            <BiotechIcon sx={{ fontSize: 40, mr: 2, color: theme.palette.primary.main }} />
            <Box>
              <Typography 
                variant="h5" 
                component="div" 
                sx={{ 
                  fontWeight: 700,
                  background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                  WebkitBackgroundClip: 'text',
                  WebkitTextFillColor: 'transparent',
                  backgroundClip: 'text',
                }}
              >
                PCI Prediction Platform
              </Typography>
              <Typography variant="caption" sx={{ color: 'text.secondary' }}>
                단백질-화합물 상호작용 예측 플랫폼
              </Typography>
            </Box>
          </Box>
          <ScienceIcon sx={{ fontSize: 40, color: theme.palette.secondary.main, opacity: 0.8 }} />
        </Toolbar>
      </AppBar>

      {/* Main Content */}
      <Container maxWidth="xl" sx={{ py: 6 }}>
        {/* Stepper */}
        <Paper 
          elevation={3}
          sx={{ 
            p: 4, 
            mb: 4,
            borderRadius: 4,
            background: 'rgba(255, 255, 255, 0.98)',
          }}
        >
          <Stepper activeStep={activeStep} alternativeLabel>
            {steps.map((label) => (
              <Step key={label}>
                <StepLabel
                  sx={{
                    '& .MuiStepLabel-label': {
                      fontWeight: 500,
                      fontSize: '0.95rem',
                    },
                  }}
                >
                  {label}
                </StepLabel>
              </Step>
            ))}
          </Stepper>
        </Paper>

        {/* Step Content */}
        <Paper 
          elevation={3}
          className="fade-in"
          sx={{ 
            p: 4,
            borderRadius: 4,
            background: 'rgba(255, 255, 255, 0.98)',
            minHeight: 500,
          }}
        >
          {getStepContent(activeStep)}
        </Paper>

        {/* Info Footer */}
        <Box sx={{ mt: 4, textAlign: 'center' }}>
          <Typography variant="body2" sx={{ color: 'rgba(255, 255, 255, 0.9)' }}>
            Powered by ChEMBL, BindingDB, FooDB | Built with React + Flask + TabPFN
          </Typography>
          <Typography variant="caption" sx={{ color: 'rgba(255, 255, 255, 0.7)' }}>
            Version 2.0.0 - Professional Edition
          </Typography>
        </Box>
      </Container>
    </Box>
  )
}

export default App
