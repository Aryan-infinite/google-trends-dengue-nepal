// Third, independent implementation of Supplementary Figure S1 values (JavaScript). Reads the preprint .xlsx by
// unzipping it (jszip) and parsing the worksheet XML; ACF with denominator n; PACF as the lag-k coefficient of an OLS
// regression on a constant and lags 1..k (normal equations, Gaussian elimination); Ljung-Box Q at lags 6 and 12.
// Usage: node scripts/52_acf_independent.js   -> out/figS1_values_js.csv
const fs = require('fs'); const JSZip = require('jszip');
const col = r => { let n = 0; for (const ch of r.match(/^[A-Z]+/)[0]) n = n * 26 + ch.charCodeAt(0) - 64; return n - 1; };
async function cases(path) {
  const z = await JSZip.loadAsync(fs.readFileSync(path)), sh = [];
  const ss = z.file('xl/sharedStrings.xml'); if (ss) for (const si of (await ss.async('string')).match(/<si>[\s\S]*?<\/si>/g) || [])
    sh.push((si.match(/<t[^>]*>([\s\S]*?)<\/t>/g) || []).map(t => t.replace(/<[^>]+>/g, '')).join(''));
  const wb = await z.file('xl/workbook.xml').async('string'), rels = await z.file('xl/_rels/workbook.xml.rels').async('string');
  const rid = wb.match(/<sheet [^>]*r:id="([^"]+)"/)[1];
  const tgt = (rels.match(new RegExp(`Id="${rid}"[^>]*Target="([^"]+)"`)) || rels.match(new RegExp(`Target="([^"]+)"[^>]*Id="${rid}"`)))[1];
  const xml = await z.file('xl/' + tgt.replace(/^\/?xl\//, '')).async('string'), rows = {};
  for (const row of xml.match(/<row [\s\S]*?<\/row>/g)) { const r = +row.match(/ r="(\d+)"/)[1]; rows[r] = {};
    for (const c of row.match(/<c [^>]*?(?:\/>|>[\s\S]*?<\/c>)/g) || []) { const v = (c.match(/<v>([\s\S]*?)<\/v>/) || [])[1]; if (v === undefined) continue;
      const t = (c.match(/ t="(\w+)"/) || [])[1]; rows[r][col(c.match(/ r="([A-Z]+\d+)"/)[1])] = t === 's' ? sh[+v] : +v; } }
  const rn = Object.keys(rows).map(Number).sort((a, b) => a - b), hr = rn.find(r => Object.values(rows[r]).includes('Confirmed Cases'));
  const k = +Object.entries(rows[hr]).find(([, v]) => v === 'Confirmed Cases')[0];
  const out = rn.filter(r => r > hr && typeof rows[r][k] === 'number').map(r => rows[r][k]); if (out.length !== 72) throw new Error(out.length); return out;
}
const trend = p => fs.readFileSync(p, 'utf8').split(/\r?\n/).filter(l => /^\d{4}-\d{2},/.test(l)).map(l => { const v = l.split(',')[1].trim(); return v === '<1' ? 0.5 : +v; });
function solve(A, b) { const n = b.length, M = A.map((r, i) => [...r, b[i]]);
  for (let k = 0; k < n; k++) { let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i; [M[k], M[p]] = [M[p], M[k]];
    for (let i = k + 1; i < n; i++) { const f = M[i][k] / M[k][k]; for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j]; } }
  const x = Array(n).fill(0); for (let i = n - 1; i >= 0; i--) { let s = M[i][n]; for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j]; x[i] = s / M[i][i]; } return x; }
(async () => {
  const c = (await cases(process.env.PREPRINT_XLSX || 'data_raw/paper_finalized_data_RP.xlsx')).slice(12);
  const t = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv', 'data_raw/gt_dengue_term_day3_2026-11-23.csv'].map(trend);
  const s = t[0].map((_, i) => (t[0][i] + t[1][i] + t[2][i]) / 3).slice(12), n = 60, K = 20;
  const mon = i => i % 12, adj = v => { const m = Array(12).fill(0); v.forEach((x, i) => m[mon(i)] += x / 5); return v.map((x, i) => x - m[mon(i)]); };
  const ly = c.map(x => Math.log(x + 1)), lx = s.map(x => Math.log(x + 0.5));
  const S = { ly_raw: ly, ly_adj: adj(ly), lx_raw: lx, lx_adj: adj(lx) }, out = ['series,stat,value'];
  for (const [key, v] of Object.entries(S)) {
    const mu = v.reduce((a, b) => a + b) / n, d = v.map(x => x - mu), den = d.reduce((a, x) => a + x * x, 0);
    const acf = []; for (let k = 1; k <= K; k++) { let s_ = 0; for (let i = 0; i + k < n; i++) s_ += d[i] * d[i + k]; acf.push(s_ / den); }
    const pacf = [];
    for (let k = 1; k <= K; k++) { const A = Array.from({ length: k + 1 }, () => Array(k + 1).fill(0)), b = Array(k + 1).fill(0);
      for (let i = k; i < n; i++) { const x = [1]; for (let j = 1; j <= k; j++) x.push(v[i - j]);
        for (let p = 0; p <= k; p++) { b[p] += x[p] * v[i]; for (let q = 0; q <= k; q++) A[p][q] += x[p] * x[q]; } }
      pacf.push(solve(A, b)[k]); }
    for (const L of [6, 12]) { let q = 0; for (let k = 1; k <= L; k++) q += acf[k - 1] ** 2 / (n - k); out.push(`${key},lbq${L},${(n * (n + 2) * q).toPrecision(17)}`); }
    acf.forEach((a, k) => out.push(`${key},acf${k + 1},${a.toPrecision(17)}`)); pacf.forEach((a, k) => out.push(`${key},pacf${k + 1},${a.toPrecision(17)}`));
  }
  fs.writeFileSync('out/figS1_values_js.csv', out.join('\n') + '\n'); console.log(`wrote ${out.length - 1} values`);
})().catch(e => { console.error(e); process.exit(1); });
