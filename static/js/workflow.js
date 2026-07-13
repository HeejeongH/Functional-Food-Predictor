/* =============================================================
   PCI Prediction — Workflow Controller v4
   5-Page SPA: data-prep / train / results / shap / predict
   ============================================================= */

const API = window.location.origin;

/* ──────────────────────────────────────────────────────────────
   STATE
────────────────────────────────────────────────────────────── */
let collectedProteins   = [];   // proteins gathered in Page 1
let lastTrainResults    = [];   // [{model_id, model_type, protein, metrics…}]
let foodbUploadedPath   = null; // server path after upload
let foodbUploadedSmiles = [];   // SMILES extracted from CSV client-side
let predResults         = [];   // last prediction rows for CSV export
let aucChart            = null;
let metricsChart        = null;
let shapChart           = null;

const PAGE_META = {
  'data-prep': ['01', '데이터 준비'],
  'train':     ['02', '모델 학습'],
  'results':   ['03', '결과 비교'],
  'shap':      ['04', 'SHAP 분석'],
  'predict':   ['05', '예측'],
};

/* ──────────────────────────────────────────────────────────────
   BOOT
────────────────────────────────────────────────────────────── */
document.addEventListener('DOMContentLoaded', () => {
  initNav();
  initModelCards();
  initDragDrop();
  checkHealth();
  setInterval(checkHealth, 30_000);
});

/* ──────────────────────────────────────────────────────────────
   NAVIGATION
────────────────────────────────────────────────────────────── */
function initNav() {
  document.querySelectorAll('.sb-item[data-page]').forEach(el => {
    el.addEventListener('click', e => {
      e.preventDefault();
      navigateTo(el.dataset.page);
    });
  });

  document.getElementById('burger')?.addEventListener('click', () => {
    document.getElementById('sidebar').classList.toggle('open');
  });

  // close sidebar on outside click (mobile)
  document.addEventListener('click', e => {
    const sb = document.getElementById('sidebar');
    if (sb?.classList.contains('open') &&
        !sb.contains(e.target) &&
        !e.target.closest('#burger')) {
      sb.classList.remove('open');
    }
  });
}

function navigateTo(page) {
  // pages
  document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
  const target = document.getElementById('page-' + page);
  if (target) target.classList.add('active');

  // sidebar items
  document.querySelectorAll('.sb-item[data-page]').forEach(el => {
    el.classList.toggle('active', el.dataset.page === page);
  });

  // topbar
  const meta = PAGE_META[page] || ['--', page];
  const stepEl  = document.getElementById('tb-step-num');
  const titleEl = document.getElementById('tb-title');
  if (stepEl)  stepEl.textContent  = meta[0];
  if (titleEl) titleEl.textContent = meta[1];

  // close sidebar on mobile
  document.getElementById('sidebar')?.classList.remove('open');

  // lazy-load data on page enter
  if (page === 'results') loadCompareModels();
  if (page === 'shap')    loadShapModels();
  if (page === 'predict') loadPredictModels();
}

/* ──────────────────────────────────────────────────────────────
   HEALTH CHECK
────────────────────────────────────────────────────────────── */
async function checkHealth() {
  const dot   = document.getElementById('sys-dot');
  const label = document.getElementById('sys-label');
  try {
    const { data } = await axios.get(`${API}/health`);
    dot.className = 'sys-dot ok';
    if (label) label.textContent = '시스템 정상';
    setEl('st-ds',  data.data?.collected_datasets ?? 0);
    setEl('st-mdl', data.data?.trained_models     ?? 0);
  } catch {
    dot.className = 'sys-dot err';
    if (label) label.textContent = '연결 오류';
  }
}

/* ──────────────────────────────────────────────────────────────
   ████  PAGE 1 — 데이터 준비
────────────────────────────────────────────────────────────── */

/* ── A: 데이터 수집 ── */
async function runCollect() {
  const raw     = getVal('inp-proteins').trim();
  const stdtype = getVal('inp-stdtype');
  if (!raw) { toast('타겟 단백질 이름을 입력하세요', 'warn'); return; }

  const btn = byId('btn-collect');
  setBtnLoading(btn, '수집 중…');
  setStatus('st-collect', 'running', '수집 중');

  const res = byId('res-collect');
  showBox(res, 'info', '<i class="fas fa-spinner fa-spin"></i> ChEMBL + BindingDB 에서 데이터를 가져오는 중입니다…');

  try {
    const targets = raw.split(',').map(s => s.trim()).filter(Boolean);
    const { data } = await axios.post(`${API}/api/data/collect`, {
      target_list:   targets,
      standard_type: stdtype,
    });

    collectedProteins = targets;
    populateProteinDropdowns(targets);

    showBox(res, 'ok',
      `<strong><i class="fas fa-check-circle"></i> 수집 완료</strong><br>
       단백질: <em>${targets.join(', ')}</em><br>
       총 화합물: <strong>${fmt(data.total_compounds)}</strong>개
       &nbsp;<span class="muted">(ChEMBL: ${fmt(data.chembl_count)} / BindingDB: ${fmt(data.bindingdb_count)})</span>`
    );
    setStatus('st-collect', 'done', '완료');
    checkHealth();
    toast('데이터 수집 완료!', 'ok');

  } catch (err) {
    showBox(res, 'err', '수집 실패: ' + extractErr(err));
    setStatus('st-collect', 'err', '실패');
    toast('수집 실패', 'err');
  } finally {
    setBtnReady(btn, '<i class="fas fa-cloud-download-alt"></i> 수집 시작');
  }
}

/* ── B: Decoy 생성 ── */
async function runDecoy() {
  const protein = getVal('inp-decoy-protein');
  if (!protein) { toast('단백질을 선택하세요 (A 단계 수집 필요)', 'warn'); return; }

  const btn = byId('btn-decoy');
  setBtnLoading(btn, 'Decoy 생성 중…');
  setStatus('st-decoy', 'running', '생성 중');

  const res = byId('res-decoy');
  showBox(res, 'info', '<i class="fas fa-spinner fa-spin"></i> DUD-E 방식으로 Decoy 화합물을 생성하는 중입니다…');

  try {
    const { data } = await axios.post(`${API}/api/decoy/generate`, {
      protein_name:  protein,
      decoy_ratio:   parseFloat(getVal('inp-decoy-ratio')),
      pos_threshold: parseFloat(getVal('inp-pos-thr')),
      neg_threshold: parseFloat(getVal('inp-neg-thr')),
    });

    showBox(res, 'ok',
      `<strong><i class="fas fa-check-circle"></i> Decoy 생성 완료</strong><br>
       <div class="mini-stat-row">
         ${mStat(data.num_actives, '활성', '#10b981')}
         ${mStat(data.num_inactives_real, '실 비활성', '#6366f1')}
         ${mStat(data.num_decoys_generated, 'Decoy', '#f59e0b')}
         ${mStat(data.total_compounds, '합계', '#0ea5e9')}
       </div>
       <span class="muted">최종 비율: ${data.final_ratio ?? '–'}</span>`
    );
    setStatus('st-decoy', 'done', '완료');
    toast('Decoy 생성 완료!', 'ok');

  } catch (err) {
    showBox(res, 'err', 'Decoy 생성 실패: ' + extractErr(err));
    setStatus('st-decoy', 'err', '실패');
    toast('Decoy 생성 실패', 'err');
  } finally {
    setBtnReady(btn, '<i class="fas fa-balance-scale"></i> Decoy 생성');
  }
}

/* ── C: 특성 변환 ── */
async function runTransform() {
  const protein = getVal('inp-trans-protein');
  if (!protein) { toast('단백질을 선택하세요 (A 단계 수집 필요)', 'warn'); return; }

  const btn = byId('btn-transform');
  setBtnLoading(btn, '변환 중…');
  setStatus('st-transform', 'running', '변환 중');

  const res = byId('res-transform');
  showBox(res, 'info', '<i class="fas fa-spinner fa-spin"></i> SMILES → Fingerprint + Descriptor 변환 중입니다… (수 분 소요될 수 있음)');

  try {
    const { data } = await axios.post(`${API}/api/features/transform`, {
      protein_name:     protein,
      fingerprint_type: getVal('inp-fptype'),
      dataset_type:     'dataset',
      dataset_ratio:    getVal('inp-dsr'),
      ignore3D:         getVal('inp-3d') === 'true',
    });

    const fpSize  = data.fingerprint_size  ?? '–';
    const descCnt = data.descriptor_count  ?? '–';
    const total   = data.feature_count     ?? (fpSize !== '–' && descCnt !== '–' ? fpSize + descCnt : '–');

    showBox(res, 'ok',
      `<strong><i class="fas fa-check-circle"></i> 특성 변환 완료</strong><br>
       <div class="mini-stat-row">
         ${mStat(fpSize,             'Fingerprint',   '#8b5cf6')}
         ${mStat(descCnt,            'Descriptor',    '#06b6d4')}
         ${mStat(total,              '총 특성',       '#0ea5e9')}
         ${mStat(data.total_compounds ?? '–', '화합물', '#10b981')}
       </div>
       <span class="muted">타입: ${data.fingerprint_type ?? getVal('inp-fptype')} · ${getVal('inp-3d')==='true'?'2D만':'2D+3D'} · ${getVal('inp-dsr')}</span>`
    );
    setStatus('st-transform', 'done', '완료');
    checkHealth();
    toast('특성 변환 완료!', 'ok');

  } catch (err) {
    showBox(res, 'err', '특성 변환 실패: ' + extractErr(err));
    setStatus('st-transform', 'err', '실패');
    toast('특성 변환 실패', 'err');
  } finally {
    setBtnReady(btn, '<i class="fas fa-cogs"></i> 변환 실행');
  }
}

/* populate all protein <select> elements across pages */
function populateProteinDropdowns(proteins) {
  const opts = proteins.map(p => `<option value="${p}">${p}</option>`).join('');
  ['inp-decoy-protein', 'inp-trans-protein', 'inp-train-protein', 'inp-shap-protein-filter'].forEach(id => {
    const el = byId(id);
    if (!el) return;
    el.disabled = false;
    el.innerHTML = opts;
  });
}

/* ──────────────────────────────────────────────────────────────
   ████  PAGE 2 — 모델 학습
────────────────────────────────────────────────────────────── */

function initModelCards() {
  document.querySelectorAll('.model-card').forEach(card => {
    const cb = card.querySelector('input[type=checkbox]');
    if (!cb) return;
    // reflect initial checked state
    if (cb.checked) card.classList.add('checked');
    card.addEventListener('click', e => {
      // let the native checkbox toggle first
      setTimeout(() => {
        card.classList.toggle('checked', cb.checked);
      }, 0);
    });
  });
}

async function runTrainAll() {
  const protein  = getVal('inp-train-protein');
  const testSize = parseFloat(getVal('inp-test-size'))  || 0.2;
  const seed     = parseInt(getVal('inp-seed'))         || 42;

  if (!protein) { toast('단백질을 선택하세요 (데이터 준비 페이지에서 수집 필요)', 'warn'); return; }

  // gather checked models
  const checked = [...document.querySelectorAll('.model-card input[type=checkbox]:checked')]
    .map(cb => cb.value);

  if (!checked.length) { toast('학습할 모델을 하나 이상 선택하세요', 'warn'); return; }

  // show progress card
  const progCard = byId('train-progress-card');
  const logEl    = byId('train-log');
  const spinEl   = byId('train-spin');
  progCard?.classList.remove('hidden');
  if (logEl) logEl.innerHTML = '';

  const btn = byId('btn-train');
  setBtnLoading(btn, '학습 중…');

  lastTrainResults = [];

  for (const modelType of checked) {
    const logItem = appendLogItem(logEl, modelType, 'running');

    try {
      const { data } = await axios.post(`${API}/api/models/train`, {
        protein_name: protein,
        model_type:   modelType,
        feature_type: 'fingerprint',
        test_size:    testSize,
        random_state: seed,
      });

      // store result
      lastTrainResults.push({
        model_id:    data.model_id,
        model_type:  modelType,
        protein:     protein,
        accuracy:    data.accuracy,
        precision:   data.precision,
        recall:      data.recall,
        f1_score:    data.f1_score,
        roc_auc:     data.roc_auc,
      });

      updateLogItem(logItem, 'done',
        `AUC ${pct(data.roc_auc)} · F1 ${pct(data.f1_score)} · Acc ${pct(data.accuracy)}`
      );
      toast(`${modelType} 학습 완료 (AUC ${pct(data.roc_auc)})`, 'ok');

    } catch (err) {
      updateLogItem(logItem, 'err', extractErr(err));
      toast(`${modelType} 학습 실패`, 'err');
    }
  }

  // finish
  if (spinEl) spinEl.className = 'fas fa-check-circle';
  setBtnReady(btn, '<i class="fas fa-play-circle"></i><span>선택 모델 학습 시작</span>');

  checkHealth();
  populateShapModels();
  populatePredModels();

  if (lastTrainResults.length) {
    toast('모든 학습 완료! 결과 비교 페이지를 확인하세요.', 'ok');
  }
}

function appendLogItem(container, name, state) {
  const el = document.createElement('div');
  el.className = `log-item ${state}`;
  el.innerHTML = `
    <span class="log-icon"><i class="fas fa-spinner fa-spin"></i></span>
    <span class="log-name">${name}</span>
    <span class="log-detail">학습 중…</span>`;
  container?.appendChild(el);
  container?.scrollTo(0, container.scrollHeight);
  return el;
}

function updateLogItem(el, state, detail) {
  el.className = `log-item ${state}`;
  const iconMap = { done: 'fa-check-circle', err: 'fa-times-circle', running: 'fa-spinner fa-spin' };
  el.querySelector('.log-icon').innerHTML = `<i class="fas ${iconMap[state] || 'fa-circle'}"></i>`;
  el.querySelector('.log-detail').textContent = detail;
}

/* ──────────────────────────────────────────────────────────────
   ████  PAGE 3 — 결과 비교
────────────────────────────────────────────────────────────── */

async function loadCompareModels() {
  const listEl = byId('compare-model-list');
  if (listEl) listEl.innerHTML = '<span class="muted"><i class="fas fa-spinner fa-spin"></i> 불러오는 중…</span>';

  try {
    const { data } = await axios.get(`${API}/api/models/list`);
    const models = data.models ?? [];

    if (!models.length) {
      if (listEl) listEl.innerHTML = '<span class="muted">저장된 모델이 없습니다. 모델 학습 페이지에서 먼저 학습하세요.</span>';
      clearCompareCharts();
      return;
    }

    // render checkboxes
    if (listEl) {
      listEl.innerHTML = models.map(m =>
        `<label class="cmp-check">
           <input type="checkbox" value="${m.model_id}" checked onchange="renderCompareCharts()"/>
           <span>${m.model_type ?? m.model_id}</span>
           <span class="cmp-protein">${m.protein_name ?? ''}</span>
         </label>`
      ).join('');
    }

    // merge with in-memory lastTrainResults
    models.forEach(m => {
      if (!lastTrainResults.find(r => r.model_id === m.model_id)) {
        lastTrainResults.push({
          model_id:   m.model_id,
          model_type: m.model_type,
          protein:    m.protein_name,
          accuracy:   m.metrics?.accuracy,
          precision:  m.metrics?.precision,
          recall:     m.metrics?.recall,
          f1_score:   m.metrics?.f1_score,
          roc_auc:    m.metrics?.roc_auc,
          train_size: m.train_size,
          test_size:  m.test_size,
        });
      }
    });

    renderCompareCharts();

  } catch (err) {
    if (listEl) listEl.innerHTML = `<span class="err-text">불러오기 실패: ${extractErr(err)}</span>`;
  }
}

function getSelectedModelIds() {
  return [...document.querySelectorAll('#compare-model-list input[type=checkbox]:checked')]
    .map(cb => cb.value);
}

function renderCompareCharts() {
  const selectedIds = getSelectedModelIds();
  const rows = lastTrainResults.filter(r => selectedIds.includes(r.model_id));

  if (!rows.length) { clearCompareCharts(); return; }

  const labels = rows.map(r => r.model_type || r.model_id.split('_')[1] || r.model_id);

  // ── AUC chart ──
  const aucCtx = byId('chart-auc')?.getContext('2d');
  if (aucCtx) {
    aucChart?.destroy();
    aucChart = new Chart(aucCtx, {
      type: 'bar',
      data: {
        labels,
        datasets: [{
          label: 'ROC-AUC',
          data: rows.map(r => r.roc_auc ?? 0),
          backgroundColor: rows.map((_, i) => CHART_COLORS[i % CHART_COLORS.length] + 'cc'),
          borderColor:     rows.map((_, i) => CHART_COLORS[i % CHART_COLORS.length]),
          borderWidth: 2,
          borderRadius: 6,
        }],
      },
      options: chartOpts('ROC-AUC', 0, 1),
    });
  }

  // ── F1 / Precision / Recall chart ──
  const metCtx = byId('chart-metrics')?.getContext('2d');
  if (metCtx) {
    metricsChart?.destroy();
    metricsChart = new Chart(metCtx, {
      type: 'bar',
      data: {
        labels,
        datasets: [
          { label: 'F1',        data: rows.map(r => r.f1_score  ?? 0), backgroundColor: '#8b5cf6cc', borderColor: '#8b5cf6', borderWidth: 2, borderRadius: 6 },
          { label: 'Precision', data: rows.map(r => r.precision ?? 0), backgroundColor: '#0ea5e9cc', borderColor: '#0ea5e9', borderWidth: 2, borderRadius: 6 },
          { label: 'Recall',    data: rows.map(r => r.recall    ?? 0), backgroundColor: '#10b981cc', borderColor: '#10b981', borderWidth: 2, borderRadius: 6 },
        ],
      },
      options: chartOpts('Metrics', 0, 1),
    });
  }

  // ── table ──
  const tbody = byId('metrics-tbody');
  if (tbody) {
    tbody.innerHTML = rows.map(r => `
      <tr>
        <td style="font-size:.78rem;font-family:monospace">${r.model_id}</td>
        <td>${r.model_type ?? '–'}</td>
        <td>${r.protein ?? '–'}</td>
        <td>${fmtPct(r.accuracy)}</td>
        <td>${fmtPct(r.precision)}</td>
        <td>${fmtPct(r.recall)}</td>
        <td><strong>${fmtPct(r.f1_score)}</strong></td>
        <td><strong style="color:var(--p3)">${fmtPct(r.roc_auc)}</strong></td>
        <td>${fmt(r.train_size)}</td>
        <td>${fmt(r.test_size)}</td>
      </tr>`).join('');
  }
}

function clearCompareCharts() {
  aucChart?.destroy(); aucChart = null;
  metricsChart?.destroy(); metricsChart = null;
  const tbody = byId('metrics-tbody');
  if (tbody) tbody.innerHTML = '<tr><td colspan="10" class="empty-cell">모델을 불러오세요</td></tr>';
}

/* ──────────────────────────────────────────────────────────────
   ████  PAGE 4 — SHAP 분석
────────────────────────────────────────────────────────────── */

async function loadShapModels() {
  populateShapModels();
}

async function populateShapModels() {
  const sel = byId('inp-shap-model');
  if (!sel) return;

  try {
    const { data } = await axios.get(`${API}/api/models/list`);
    const models = data.models ?? [];
    sel.innerHTML = models.length
      ? models.map(m => `<option value="${m.model_id}">${m.model_type ?? m.model_id} — ${m.protein_name ?? ''}</option>`).join('')
      : '<option value="">학습된 모델 없음</option>';
  } catch {
    sel.innerHTML = '<option value="">모델 목록 로드 실패</option>';
  }
}

async function runShap() {
  const modelId = getVal('inp-shap-model');
  const topN    = parseInt(getVal('inp-shap-n')) || 20;

  if (!modelId) { toast('분석할 모델을 선택하세요', 'warn'); return; }

  const btn = byId('btn-shap');
  setBtnLoading(btn, 'SHAP 분석 중…');

  const statusEl = byId('res-shap-status');
  showBox(statusEl, 'info', '<i class="fas fa-spinner fa-spin"></i> SHAP 분석 실행 중… (모델 크기에 따라 수 분 소요)');

  byId('shap-result-card')?.classList.add('hidden');

  try {
    const { data } = await axios.post(`${API}/api/shap/analyze`, {
      model_id:     modelId,
      feature_type: 'fingerprint',
      top_n:        topN,
    });

    showBox(statusEl, 'ok', `<i class="fas fa-check-circle"></i> SHAP 분석 완료 (샘플 ${data.shap_values_summary?.samples_analyzed ?? '?'}개 분석)`);

    // label
    const labelEl = byId('shap-model-label');
    if (labelEl) labelEl.textContent = `모델: ${modelId}`;

    renderShapChart(data.top_features, topN);
    renderShapTable(data.top_features);
    byId('shap-result-card')?.classList.remove('hidden');
    toast('SHAP 분석 완료!', 'ok');

  } catch (err) {
    showBox(statusEl, 'err', 'SHAP 분석 실패: ' + extractErr(err));
    toast('SHAP 분석 실패', 'err');
  } finally {
    setBtnReady(btn, '<i class="fas fa-microscope"></i> SHAP 분석 실행');
  }
}

function renderShapChart(features, topN) {
  const ctx = byId('chart-shap')?.getContext('2d');
  if (!ctx || !features?.length) return;

  const top = features.slice(0, Math.min(topN, 20));
  const labels = top.map(f => {
    const name = f.feature ?? f.name ?? f.feature_name ?? `Feature_${f.feature_index ?? '?'}`;
    return name.length > 22 ? name.slice(0, 22) + '…' : name;
  });
  const vals = top.map(f => parseFloat(f.importance ?? f.shap_value ?? f.mean_shap ?? 0));

  shapChart?.destroy();
  shapChart = new Chart(ctx, {
    type: 'bar',
    data: {
      labels,
      datasets: [{
        label: 'SHAP 값 (절댓값)',
        data: vals,
        backgroundColor: vals.map((_, i) => `hsl(${195 + i * 8}, 80%, 55%)99`),
        borderColor: '#06b6d4',
        borderWidth: 1,
        borderRadius: 4,
      }],
    },
    options: {
      indexAxis: 'y',
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        tooltip: { callbacks: { label: ctx => ` ${ctx.parsed.x.toFixed(5)}` } },
      },
      scales: {
        x: { grid: { color: '#e2e8f033' }, ticks: { color: '#94a3b8', font: { size: 11 } } },
        y: { grid: { display: false },     ticks: { color: '#cbd5e1', font: { size: 11 } } },
      },
    },
  });
}

function renderShapTable(features) {
  const tbody = byId('shap-tbody');
  if (!tbody || !features?.length) return;

  const maxV = Math.max(...features.map(f => Math.abs(f.importance ?? f.shap_value ?? f.mean_shap ?? 0)));

  tbody.innerHTML = features.slice(0, 30).map((f, i) => {
    const name = f.feature ?? f.name ?? f.feature_name ?? `Feature_${f.feature_index ?? i}`;
    const val  = parseFloat(f.importance ?? f.shap_value ?? f.mean_shap ?? 0);
    const rel  = maxV ? (Math.abs(val) / maxV * 100).toFixed(1) : '0';
    return `<tr>
      <td class="muted">${i + 1}</td>
      <td style="font-size:.8rem;max-width:180px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="${name}">${name}</td>
      <td style="font-family:monospace;font-size:.8rem">${Math.abs(val).toFixed(5)}</td>
      <td>
        <div style="display:flex;align-items:center;gap:.5rem">
          <div style="flex:1;height:6px;background:#1e293b;border-radius:3px">
            <div style="width:${rel}%;height:100%;background:var(--p4);border-radius:3px"></div>
          </div>
          <span style="font-size:.75rem;color:#94a3b8;width:40px">${rel}%</span>
        </div>
      </td>
    </tr>`;
  }).join('');
}

/* ──────────────────────────────────────────────────────────────
   ████  PAGE 5 — 예측
────────────────────────────────────────────────────────────── */

function switchTab(tab, btn) {
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
  btn.classList.add('active');
  byId('tab-' + tab)?.classList.add('active');
}

async function loadPredictModels() {
  populatePredModels();
}

async function populatePredModels() {
  const sel = byId('inp-pred-model');
  if (!sel) return;

  try {
    const { data } = await axios.get(`${API}/api/models/list`);
    const models = data.models ?? [];
    sel.innerHTML = models.length
      ? models.map(m => `<option value="${m.model_id}">${m.model_type ?? m.model_id} — ${m.protein_name ?? ''}</option>`).join('')
      : '<option value="">새로고침 버튼을 눌러 모델을 불러오세요</option>';
  } catch {
    sel.innerHTML = '<option value="">모델 목록 로드 실패</option>';
  }
}

function loadSampleSmiles() {
  byId('inp-smiles').value =
`CCO
CC(=O)O
c1ccccc1C(=O)O
CN1C=NC2=C1C(=O)N(C(=O)N2C)C
CC(C)Cc1ccc(cc1)C(C)C(=O)O
OC(=O)c1ccc(cc1)N
c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34
CC(=O)Nc1ccc(O)cc1
OC(=O)C1CCCN1
c1ccc(cc1)-c1ccncc1`;
  toast('샘플 SMILES 10개 불러옴', 'info');
}

async function runSmilesPred() {
  const modelId = getVal('inp-pred-model');
  if (!modelId) { toast('예측에 사용할 모델을 선택하세요', 'warn'); return; }

  const raw = getVal('inp-smiles').trim();
  if (!raw)  { toast('SMILES를 입력하세요', 'warn'); return; }

  const smilesList = raw.split('\n').map(s => s.trim()).filter(Boolean);
  if (!smilesList.length) { toast('유효한 SMILES가 없습니다', 'warn'); return; }

  const btn = byId('btn-smiles-pred');
  setBtnLoading(btn, '예측 중…');
  setStatus('st-smiles', 'running', '예측 중');
  byId('pred-result-card')?.classList.add('hidden');

  try {
    const { data } = await axios.post(`${API}/api/models/predict`, {
      smiles_list:  smilesList,
      model_id:     modelId,
      feature_type: 'fingerprint',
    });

    predResults = data.predictions ?? [];
    renderPredResults(predResults, `SMILES ${smilesList.length}개 예측 완료`);
    setStatus('st-smiles', 'done', '완료');
    toast(`예측 완료: ${smilesList.length}개`, 'ok');

  } catch (err) {
    setStatus('st-smiles', 'err', '실패');
    toast('예측 실패: ' + extractErr(err), 'err');
  } finally {
    setBtnReady(btn, '<i class="fas fa-play-circle"></i> 예측 실행');
  }
}

/* ── FooDB CSV upload ── */
function onCsvSelect(e) {
  const file = e.target.files?.[0];
  if (!file) return;

  const infoEl = byId('csv-info');
  showBox(infoEl, 'info', `<i class="fas fa-spinner fa-spin"></i> "${file.name}" 읽는 중…`);

  const reader = new FileReader();
  reader.onload = ev => {
    const text  = ev.target.result;
    const lines = text.split('\n').filter(Boolean);
    const header = lines[0]?.split(',').map(c => c.trim().replace(/"/g, '').toLowerCase());

    // detect SMILES column
    const smilesIdx = header?.findIndex(c =>
      ['moldb_smiles','canonical_smiles','smiles'].includes(c)
    );

    if (smilesIdx < 0) {
      showBox(infoEl, 'err', `SMILES 컬럼을 찾을 수 없습니다.<br>사용 가능한 컬럼: <em>${header?.join(', ')}</em>`);
      return;
    }

    foodbUploadedSmiles = lines.slice(1).map(line => {
      const cols = line.split(',');
      return cols[smilesIdx]?.trim().replace(/"/g, '') ?? '';
    }).filter(Boolean);

    showBox(infoEl, 'ok',
      `<i class="fas fa-check-circle"></i> <strong>${file.name}</strong> 로드 완료<br>
       화합물: <strong>${foodbUploadedSmiles.length.toLocaleString()}</strong>개
       &nbsp;(SMILES 컬럼: <em>${header[smilesIdx]}</em>)`
    );

    // also upload to server for server-side prediction
    uploadCsvToServer(file);
  };
  reader.readAsText(file);
}

async function uploadCsvToServer(file) {
  try {
    const form = new FormData();
    form.append('file', file);
    const { data } = await axios.post(`${API}/api/foodb/upload`, form);
    foodbUploadedPath = data.file_path;
  } catch {
    // non-fatal: we can still use client-side extracted SMILES
    foodbUploadedPath = null;
  }
}

async function runFoodbPred() {
  const modelId = getVal('inp-pred-model');
  if (!modelId) { toast('예측에 사용할 모델을 선택하세요', 'warn'); return; }

  const topN   = parseInt(getVal('inp-topn'))   || 100;
  const thresh = parseFloat(getVal('inp-thresh')) ?? 0.5;
  const mode   = getVal('inp-pred-mode');

  const btn = byId('btn-foodb-pred');

  // prefer server-side path, fallback to client-extracted SMILES list
  if (foodbUploadedPath) {
    // server-side batch predict
    setBtnLoading(btn, '예측 중…');
    setStatus('st-foodb', 'running', '예측 중');
    byId('pred-result-card')?.classList.add('hidden');

    try {
      const params = new URLSearchParams({ model_id: modelId, top_n: topN });
      const { data } = await axios.post(
        `${API}/api/foodb/predict?${params.toString()}`
      );

      predResults = data.top_active_compounds ?? [];
      const summary = `총 ${fmt(data.total_compounds)}개 중 활성 ${fmt(data.active_compounds)}개 (임계값 ${thresh})`;
      renderPredResults(predResults, summary, data);
      setStatus('st-foodb', 'done', '완료');
      toast(`FooDB 예측 완료: ${fmt(data.active_compounds)}개 활성`, 'ok');

    } catch (err) {
      setStatus('st-foodb', 'err', '실패');
      toast('FooDB 예측 실패: ' + extractErr(err), 'err');
    } finally {
      setBtnReady(btn, '<i class="fas fa-play-circle"></i> FooDB 예측');
    }

  } else if (foodbUploadedSmiles.length) {
    // client-side SMILES → /api/models/predict
    setBtnLoading(btn, '예측 중…');
    setStatus('st-foodb', 'running', '예측 중');
    byId('pred-result-card')?.classList.add('hidden');

    const smilesList = mode === 'top'
      ? foodbUploadedSmiles.slice(0, topN * 5)   // sample more than needed
      : foodbUploadedSmiles;

    try {
      const CHUNK = 200;
      let allPreds = [];

      for (let i = 0; i < smilesList.length; i += CHUNK) {
        const chunk = smilesList.slice(i, i + CHUNK);
        const { data } = await axios.post(`${API}/api/models/predict`, {
          smiles_list:  chunk,
          model_id:     modelId,
          feature_type: 'fingerprint',
        });
        allPreds = allPreds.concat(data.predictions ?? []);
      }

      // sort by active probability, filter by threshold
      allPreds.sort((a, b) => (b.probability_active ?? 0) - (a.probability_active ?? 0));
      const active   = allPreds.filter(p => (p.probability_active ?? 0) >= thresh);
      const topSlice = mode === 'top' ? active.slice(0, topN) : active;

      predResults = topSlice;
      const summary = `총 ${fmt(allPreds.length)}개 예측 · 활성(≥${thresh}): ${fmt(active.length)}개`;
      renderPredResults(topSlice, summary);
      setStatus('st-foodb', 'done', '완료');
      toast(`FooDB 예측 완료: 활성 ${fmt(active.length)}개`, 'ok');

    } catch (err) {
      setStatus('st-foodb', 'err', '실패');
      toast('FooDB 예측 실패: ' + extractErr(err), 'err');
    } finally {
      setBtnReady(btn, '<i class="fas fa-play-circle"></i> FooDB 예측');
    }

  } else {
    toast('먼저 FooDB CSV 파일을 업로드하세요', 'warn');
  }
}

/* render pred result card */
function renderPredResults(rows, summary, extra) {
  const card     = byId('pred-result-card');
  const summaryEl = byId('pred-result-summary');
  const statsEl  = byId('pred-stats');
  const tbody    = byId('pred-tbody');

  if (summaryEl) summaryEl.textContent = summary;
  card?.classList.remove('hidden');

  // stats boxes
  const total   = extra?.total_compounds ?? rows.length;
  const active  = extra?.active_compounds  ?? rows.filter(r => (r.probability_active ?? r.probability ?? 0) >= 0.5).length;
  const inactive = total - active;

  if (statsEl) {
    statsEl.innerHTML = `
      ${psBox(total,   '전체',   '#0ea5e9')}
      ${psBox(active,  '활성',   '#10b981')}
      ${psBox(inactive,'비활성', '#ef4444')}
      ${extra ? psBox(fmtPct(extra.statistics?.mean_probability_active), '평균 활성확률', '#8b5cf6') : ''}`;
  }

  // table
  if (tbody) {
    tbody.innerHTML = rows.slice(0, 200).map((r, i) => {
      const prob = r.probability_active ?? r.probability ?? 0;
      const pred = r.prediction ?? (prob >= 0.5 ? 1 : 0);
      const smiles = r.smiles ?? r.canonical_smiles ?? '–';
      const name   = r.compound_name ?? '';
      const probPct = (prob * 100).toFixed(1);
      const barColor = prob >= 0.7 ? '#10b981' : prob >= 0.5 ? '#f59e0b' : '#ef4444';

      return `<tr>
        <td class="muted">${i + 1}</td>
        <td style="font-size:.75rem;font-family:monospace;max-width:160px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="${smiles}">
          ${smiles.length > 28 ? smiles.slice(0, 28) + '…' : smiles}
          ${name ? `<br><span style="color:#94a3b8;font-family:inherit">${name}</span>` : ''}
        </td>
        <td>
          <div class="prob-bar-wrap">
            <div style="flex:1;height:6px;background:#1e293b;border-radius:3px">
              <div style="width:${probPct}%;height:100%;background:${barColor};border-radius:3px"></div>
            </div>
            <span style="font-size:.8rem;font-weight:700;color:${barColor};width:44px;text-align:right">${probPct}%</span>
          </div>
        </td>
        <td>
          <span class="pred-badge ${pred === 1 ? 'active' : 'inactive'}">${pred === 1 ? '활성' : '비활성'}</span>
        </td>
      </tr>`;
    }).join('');

    if (rows.length > 200) {
      const tr = document.createElement('tr');
      tr.innerHTML = `<td colspan="4" class="empty-cell">… +${rows.length - 200}개 (CSV 다운로드로 전체 확인)</td>`;
      tbody.appendChild(tr);
    }
  }
}

function downloadPredResult() {
  if (!predResults.length) { toast('다운로드할 결과가 없습니다', 'warn'); return; }

  const header = 'smiles,compound_name,probability_active,probability_inactive,prediction\n';
  const rows = predResults.map(r =>
    [
      `"${r.smiles ?? r.canonical_smiles ?? ''}"`,
      `"${r.compound_name ?? ''}"`,
      (r.probability_active  ?? '').toString(),
      (r.probability_inactive ?? '').toString(),
      (r.prediction ?? '').toString(),
    ].join(',')
  );
  const csv  = header + rows.join('\n');
  const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement('a');
  a.href = url;
  a.download = `pci_predictions_${Date.now()}.csv`;
  a.click();
  URL.revokeObjectURL(url);
  toast('CSV 다운로드 시작!', 'ok');
}

/* ──────────────────────────────────────────────────────────────
   DRAG-AND-DROP for FooDB upload zone
────────────────────────────────────────────────────────────── */
function initDragDrop() {
  const zone = byId('upload-zone');
  if (!zone) return;

  zone.addEventListener('dragover', e => {
    e.preventDefault();
    zone.classList.add('drag-over');
  });
  zone.addEventListener('dragleave', () => zone.classList.remove('drag-over'));
  zone.addEventListener('drop', e => {
    e.preventDefault();
    zone.classList.remove('drag-over');
    const file = e.dataTransfer.files?.[0];
    if (file && file.name.endsWith('.csv')) {
      const dt = new DataTransfer();
      dt.items.add(file);
      const inp = byId('inp-csv');
      if (inp) {
        inp.files = dt.files;
        onCsvSelect({ target: { files: [file] } });
      }
    } else {
      toast('CSV 파일만 업로드 가능합니다', 'warn');
    }
  });
}

/* ──────────────────────────────────────────────────────────────
   TOAST
────────────────────────────────────────────────────────────── */
function toast(msg, type = 'info') {
  const stack = byId('toast-stack');
  if (!stack) return;
  const el = document.createElement('div');
  const iconMap = { ok:'fa-check-circle', err:'fa-times-circle', warn:'fa-exclamation-triangle', info:'fa-info-circle' };
  el.className = `toast ${type}`;
  el.innerHTML = `<i class="fas ${iconMap[type] || 'fa-info-circle'}"></i> ${msg}`;
  stack.appendChild(el);
  // stagger entrance
  requestAnimationFrame(() => el.classList.add('show'));
  setTimeout(() => {
    el.classList.remove('show');
    setTimeout(() => el.remove(), 400);
  }, 3800);
}

/* ──────────────────────────────────────────────────────────────
   CHART HELPERS
────────────────────────────────────────────────────────────── */
const CHART_COLORS = ['#f59e0b','#10b981','#8b5cf6','#0ea5e9','#ef4444','#06b6d4','#ec4899','#84cc16'];

function chartOpts(title, minY = 0, maxY = 1) {
  return {
    responsive: true,
    maintainAspectRatio: false,
    plugins: {
      legend: { display: title.includes('/') },
      tooltip: { callbacks: { label: ctx => ` ${(ctx.parsed.y * 100).toFixed(1)}%` } },
    },
    scales: {
      x: { grid: { color: '#e2e8f022' }, ticks: { color: '#94a3b8', font: { size: 11 } } },
      y: {
        min: minY, max: maxY,
        grid: { color: '#e2e8f022' },
        ticks: { color: '#94a3b8', font: { size: 11 }, callback: v => (v * 100).toFixed(0) + '%' },
      },
    },
  };
}

/* ──────────────────────────────────────────────────────────────
   DOM HELPERS
────────────────────────────────────────────────────────────── */
const byId  = id  => document.getElementById(id);
const getVal = id => byId(id)?.value ?? '';
const setEl  = (id, v) => { const el = byId(id); if (el) el.textContent = v; };
const fmt    = v  => (v != null && v !== '–') ? Number(v).toLocaleString() : '–';
const pct    = v  => v != null ? (v * 100).toFixed(1) + '%' : '–';
const fmtPct = v  => v != null ? (v * 100).toFixed(1) + '%' : '–';

function setBtnLoading(btn, text) {
  if (!btn) return;
  btn.disabled = true;
  btn._origHTML = btn.innerHTML;
  btn.innerHTML = `<i class="fas fa-spinner fa-spin"></i> ${text}`;
}
function setBtnReady(btn, html) {
  if (!btn) return;
  btn.disabled = false;
  btn.innerHTML = html || btn._origHTML || html;
}

function showBox(el, type, html) {
  if (!el) return;
  el.className = `result-box ${type}`;
  el.innerHTML = html;
  el.classList.remove('hidden');
}

function setStatus(id, state, text) {
  const el = byId(id);
  if (!el) return;
  const iconMap = {
    idle:    '<span class="dot-idle"></span>',
    running: '<i class="fas fa-spinner fa-spin" style="color:#f59e0b"></i>',
    done:    '<i class="fas fa-check-circle" style="color:#10b981"></i>',
    err:     '<i class="fas fa-times-circle" style="color:#ef4444"></i>',
  };
  el.innerHTML = `${iconMap[state] || iconMap.idle} ${text}`;
}

/* mini stat for result boxes */
function mStat(val, label, color) {
  return `<span class="mini-stat" style="border-color:${color}22">
    <strong style="color:${color}">${fmt(val)}</strong>
    <em>${label}</em>
  </span>`;
}

/* pred stats box */
function psBox(val, label, color) {
  return `<div class="ps-box" style="border-top:3px solid ${color}">
    <div class="ps-val" style="color:${color}">${typeof val === 'string' ? val : fmt(val)}</div>
    <div class="ps-lbl">${label}</div>
  </div>`;
}

function extractErr(err) {
  if (!err.response) return err.message ?? '알 수 없는 오류';
  const d = err.response.data;
  if (!d) return err.message;
  if (typeof d.detail === 'string') return d.detail;
  if (typeof d.detail === 'object') return JSON.stringify(d.detail).slice(0, 200);
  return JSON.stringify(d).slice(0, 200);
}
