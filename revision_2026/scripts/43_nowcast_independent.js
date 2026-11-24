// Third, independent implementation of the Figure 4 nowcasts (JavaScript / Node.js).
// Shares no code with the Python (scripts/41_figure4.py) or R (scripts/42_check_figure4.R) versions:
//   - reads the preprint .xlsx by unzipping it (jszip) and parsing the worksheet XML directly;
//   - parses the three Google Trends CSVs by hand;
//   - fits OLS through the normal equations (X'X b = X'y) solved by Gaussian elimination with partial pivoting.
// Usage (from the project folder): node scripts/43_nowcast_independent.js
// Output: out/fig4_predictions_js.csv (schedule, model, month, pred) and out/fig4_mae_js.csv
const fs = require('fs');
const JSZip = require('jszip');

function colIndex(ref) {                       // "AB12" -> zero-based column index
  const letters = ref.match(/^[A-Z]+/)[0]; let n = 0;
  for (const ch of letters) n = n * 26 + (ch.charCodeAt(0) - 64);
  return n - 1;
}
function unescapeXml(s) { return s.replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/&apos;/g, "'").replace(/&amp;/g, '&'); }

async function readPreprint(path) {
  const zip = await JSZip.loadAsync(fs.readFileSync(path));
  const shared = [];
  const ss = zip.file('xl/sharedStrings.xml');
  if (ss) {
    const xml = await ss.async('string');
    for (const si of xml.match(/<si>[\s\S]*?<\/si>/g) || []) {
      shared.push(unescapeXml((si.match(/<t[^>]*>([\s\S]*?)<\/t>/g) || []).map(t => t.replace(/<[^>]+>/g, '')).join('')));
    }
  }
  const wbRels = await zip.file('xl/_rels/workbook.xml.rels').async('string');
  const wb = await zip.file('xl/workbook.xml').async('string');
  const rid = wb.match(/<sheet [^>]*r:id="([^"]+)"/)[1];                         // first sheet
  const target = wbRels.match(new RegExp(`Id="${rid}"[^>]*Target="([^"]+)"`))?.[1] || wbRels.match(new RegExp(`Target="([^"]+)"[^>]*Id="${rid}"`))[1];
  const sheet = await zip.file('xl/' + target.replace(/^\/?xl\//, '')).async('string');
  const rows = {};
  for (const row of sheet.match(/<row [\s\S]*?<\/row>/g) || []) {
    const r = Number(row.match(/ r="(\d+)"/)[1]); rows[r] = {};
    for (const c of row.match(/<c [^>]*?(?:\/>|>[\s\S]*?<\/c>)/g) || []) {
      const ref = c.match(/ r="([A-Z]+\d+)"/)[1]; const t = (c.match(/ t="(\w+)"/) || [])[1];
      const v = (c.match(/<v>([\s\S]*?)<\/v>/) || [])[1];
      if (v === undefined) continue;
      rows[r][colIndex(ref)] = t === 's' ? shared[Number(v)] : (t === 'str' ? unescapeXml(v) : Number(v));
    }
  }
  const rnums = Object.keys(rows).map(Number).sort((a, b) => a - b);
  const hr = rnums.find(r => Object.values(rows[r]).includes('Confirmed Cases'));
  const col = Number(Object.entries(rows[hr]).find(([, v]) => v === 'Confirmed Cases')[0]);
  const out = [];
  for (const r of rnums) if (r > hr && typeof rows[r][col] === 'number') out.push(rows[r][col]);
  if (out.length !== 72) throw new Error(`expected 72 case values, got ${out.length}`);
  return out;
}
function readTrends(path) {
  const lines = fs.readFileSync(path, 'utf8').split(/\r?\n/).filter(l => /^\d{4}-\d{2},/.test(l));
  if (lines.length !== 72) throw new Error(`${path}: ${lines.length} rows`);
  return lines.map(l => { const v = l.split(',')[1].trim(); return v === '<1' ? 0.5 : Number(v); });
}
function solve(A, b) {                         // Gaussian elimination, partial pivoting
  const n = b.length, M = A.map((row, i) => [...row, b[i]]);
  for (let k = 0; k < n; k++) {
    let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
    [M[k], M[p]] = [M[p], M[k]];
    for (let i = k + 1; i < n; i++) { const f = M[i][k] / M[k][k]; for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j]; }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) { let s = M[i][n]; for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j]; x[i] = s / M[i][i]; }
  return x;
}

(async () => {
  const cases = await readPreprint(process.env.PREPRINT_XLSX || 'data_raw/paper_finalized_data_RP.xlsx');
  const tr = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv', 'data_raw/gt_dengue_term_day3_2026-11-23.csv'].map(readTrends);
  const search = cases.map((_, i) => (tr[0][i] + tr[1][i] + tr[2][i]) / 3);
  const LY = cases.map(c => Math.log(c + 1)), LX = search.map(s => Math.log(s + 0.5));
  const month = i => (i % 12) + 1, year = i => 2019 + Math.floor(i / 12);
  const label = i => `${year(i)}-${String(month(i)).padStart(2, '0')}`;
  const row = (i, model) => { const r = new Array(12).fill(0); r[month(i) - 1] = 1; r.push(LY[i - 1]); if (model === 'M2') r.push(LX[i]); return r; };
  const fit = (train, model) => {
    const k = row(train[0], model).length, A = Array.from({ length: k }, () => new Array(k).fill(0)), b = new Array(k).fill(0);
    for (const i of train) { const x = row(i, model); for (let p = 0; p < k; p++) { b[p] += x[p] * LY[i]; for (let q = 0; q < k; q++) A[p][q] += x[p] * x[q]; } }
    return solve(A, b);
  };
  const pred = (beta, i, model) => Math.max(Math.expm1(row(i, model).reduce((s, v, j) => s + v * beta[j], 0)), 0);
  const START = 6, TEST = Array.from({ length: 36 }, (_, j) => 36 + j);          // July 2019; Jan 2022 .. Dec 2024
  const range = (a, b) => Array.from({ length: b - a }, (_, j) => a + j);
  const lines = ['schedule,model,month,pred'], mae = ['schedule,model,mae'];
  for (const sched of ['season', 'monthly']) {
    for (const model of ['M0', 'M1', 'M2']) {
      let err = 0;
      for (const t of TEST) {
        let p;
        if (model === 'M0') p = cases[t - 1];
        else if (sched === 'season') p = pred(fit(range(START, 36 + 12 * (year(t) - 2022)), model), t, model);
        else p = pred(fit(range(START, t), model), t, model);
        lines.push(`${sched},${model},${label(t)},${p.toPrecision(17)}`); err += Math.abs(p - cases[t]);
      }
      mae.push(`${sched},${model},${(err / TEST.length).toPrecision(17)}`);
    }
  }
  fs.writeFileSync('out/fig4_predictions_js.csv', lines.join('\n') + '\n');
  fs.writeFileSync('out/fig4_mae_js.csv', mae.join('\n') + '\n');
  console.log(mae.join('\n'));
})().catch(e => { console.error(e); process.exit(1); });
