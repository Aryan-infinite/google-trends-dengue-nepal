// Third, independent implementation of Supplementary Figure S3 values (JavaScript). Reads the RAW EDCD district workbook
// (all yearly sheets, via jszip + XML parsing), assigns each of the 77 districts to its pre-2015 development region from
// its own district list, reads the regional Google Trends files, and computes per region (2020-2023): seasonally adjusted
// correlation, effective sample size (lagged Pearson autocorrelations, lags 1-12, floor 4), Fisher 95% CI, zero months and
// 2023:2022 peak ratios.  Usage: node scripts/57_regional_independent.js
const fs = require('fs'); const JSZip = require('jszip');
const REGION = {
  Eastern: ['TAPLEJUNG', 'PANCHTHAR', 'ILAM', 'JHAPA', 'MORANG', 'SUNSARI', 'DHANKUTA', 'TERHATHUM', 'SANKHUWASABHA', 'BHOJPUR',
            'SOLUKHUMBU', 'OKHALDHUNGA', 'KHOTANG', 'UDAYAPUR', 'SAPTARI', 'SIRAHA'],
  Central: ['DHANUSA', 'DHANUSHA', 'MAHOTTARI', 'SARLAHI', 'SINDHULI', 'RAMECHHAP', 'DOLAKHA', 'SINDHUPALCHOK', 'KAVREPALANCHOK', 'KAVRE',
            'LALITPUR', 'BHAKTAPUR', 'KATHMANDU', 'NUWAKOT', 'RASUWA', 'DHADING', 'MAKWANPUR', 'MAKAWANPUR', 'RAUTAHAT', 'BARA', 'PARSA', 'CHITWAN', 'CHITAWAN'],
  Western: ['GORKHA', 'LAMJUNG', 'TANAHU', 'TANAHUN', 'SYANGJA', 'KASKI', 'MANANG', 'MUSTANG', 'MYAGDI', 'PARBAT', 'BAGLUNG', 'GULMI', 'PALPA',
            'NAWALPARASI', 'RUPANDEHI', 'KAPILBASTU', 'KAPILVASTU', 'ARGHAKHANCHI'],
  'Mid-Western': ['DANG', 'PYUTHAN', 'ROLPA', 'RUKUM', 'SALYAN', 'BANKE', 'BARDIYA', 'SURKHET', 'DAILEKH', 'JAJARKOT', 'DOLPA', 'JUMLA', 'KALIKOT', 'MUGU', 'HUMLA'],
  'Far-Western': ['BAJURA', 'BAJHANG', 'ACHHAM', 'DOTI', 'KAILALI', 'KANCHANPUR', 'DADELDHURA', 'BAITADI', 'DARCHULA'] };
const regionOf = name => { const words = name.toUpperCase().replace(/[^A-Z ]/g, ' ').split(/\s+/).filter(Boolean);
  const hit = Object.entries(REGION).filter(([, l]) => l.some(d => words.includes(d))).map(([r]) => r); return hit; };
const col = r => { let n = 0; for (const ch of r.match(/^[A-Z]+/)[0]) n = n * 26 + ch.charCodeAt(0) - 64; return n - 1; };
async function sheets(path) {
  const z = await JSZip.loadAsync(fs.readFileSync(path)), sh = [];
  const ss = z.file('xl/sharedStrings.xml'); if (ss) for (const si of (await ss.async('string')).match(/<si>[\s\S]*?<\/si>/g) || [])
    sh.push((si.match(/<t[^>]*>([\s\S]*?)<\/t>/g) || []).map(t => t.replace(/<[^>]+>/g, '')).join('').replace(/&amp;/g, '&'));
  const wb = await z.file('xl/workbook.xml').async('string'), rels = await z.file('xl/_rels/workbook.xml.rels').async('string'), out = {};
  for (const m of wb.matchAll(/<sheet [^>]*name="([^"]+)"[^>]*r:id="([^"]+)"/g)) {
    const tgt = (rels.match(new RegExp(`Id="${m[2]}"[^>]*Target="([^"]+)"`)) || rels.match(new RegExp(`Target="([^"]+)"[^>]*Id="${m[2]}"`)))[1];
    const xml = await z.file('xl/' + tgt.replace(/^\/?xl\//, '')).async('string'), rows = [];
    for (const row of xml.match(/<row [\s\S]*?<\/row>/g) || []) { const cells = {};
      for (const c of row.match(/<c [^>]*?(?:\/>|>[\s\S]*?<\/c>)/g) || []) { const v = (c.match(/<v>([\s\S]*?)<\/v>/) || [])[1]; if (v === undefined) continue;
        const t = (c.match(/ t="(\w+)"/) || [])[1]; cells[col(c.match(/ r="([A-Z]+\d+)"/)[1])] = t === 's' ? sh[+v] : +v; }
      rows.push(cells); }
    out[m[1]] = rows; }
  return out;
}
const trends = p => fs.readFileSync(p, 'utf8').split(/\r?\n/).filter(l => /^202[0-3]-\d{2},/.test(l)).map(l => { const v = l.split(',')[1].trim(); return v === '<1' ? 0.5 : +v; });
const mean = a => a.reduce((x, y) => x + y, 0) / a.length;
const pear = (a, b) => { const ma = mean(a), mb = mean(b); let s = 0, sa = 0, sb = 0; for (let i = 0; i < a.length; i++) { s += (a[i] - ma) * (b[i] - mb); sa += (a[i] - ma) ** 2; sb += (b[i] - mb) ** 2; } return s / Math.sqrt(sa * sb); };
const adj = v => { const m = Array(12).fill(0); v.forEach((x, i) => m[i % 12] += x / (v.length / 12)); return v.map((x, i) => x - m[i % 12]); };
(async () => {
  const S = await sheets('data_raw/edcd_district_monthly_2019_2023.xlsx'), RC = {}; let assigned = 0, problems = [];
  for (const r of Object.keys(REGION)) RC[r] = [];
  for (const y of [2020, 2021, 2022, 2023]) {
    const rows = S[`${y} Monthly`], hdr = rows.findIndex(r => Object.values(r).includes('District with Code'));
    const tot = {}; for (const r of Object.keys(REGION)) tot[r] = Array(12).fill(0); let nd = 0;
    for (const r of rows.slice(hdr + 1)) { const nm = r[1]; if (typeof nm !== 'string' || !/^\d{3} /.test(nm)) continue;
      const hit = regionOf(nm); if (hit.length !== 1) { problems.push(`${y} ${nm}: ${hit}`); continue; } nd++;
      for (let k = 0; k < 12; k++) tot[hit[0]][k] += typeof r[3 + k] === 'number' ? r[3 + k] : 0; }
    if (nd !== 77) problems.push(`${y}: ${nd} districts assigned`); assigned += nd;
    for (const r of Object.keys(REGION)) RC[r].push(...tot[r]); }
  if (problems.length) { console.error('district assignment problems:', problems); process.exit(1); }
  const out = ['region,r,n_eff,ci_lo,ci_hi,zeros,ratio_cases,ratio_search'];
  for (const rg of ['Central', 'Eastern', 'Western']) {
    const s = trends(`data_raw/gt_dengue_region_${rg.toLowerCase()}.csv`), c = RC[rg]; if (s.length !== 48 || c.length !== 48) throw new Error(rg);
    const a = adj(s.map(v => Math.log(v + 0.5))), b = adj(c.map(v => Math.log(v + 1))), r = pear(a, b);
    let sum = 0; for (let k = 1; k <= 12; k++) sum += pear(a.slice(0, -k), a.slice(k)) * pear(b.slice(0, -k), b.slice(k));
    const ne = Math.max(48 / (1 + 2 * sum), 4), z = Math.atanh(r), se = 1 / Math.sqrt(ne - 3);
    const mx = (v, i, j) => Math.max(...v.slice(i, j));
    out.push([rg, r, ne, Math.tanh(z - 1.96 * se), Math.tanh(z + 1.96 * se), s.filter(v => v === 0).length, mx(c, 36, 48) / mx(c, 24, 36), mx(s, 36, 48) / mx(s, 24, 36)]
      .map(v => typeof v === 'number' ? v.toPrecision(17) : v).join(','));
  }
  fs.writeFileSync('out/figS3_values_js.csv', out.join('\n') + '\n'); console.log(`${assigned / 4} districts per year assigned to exactly one region`); console.log(out.join('\n'));
})().catch(e => { console.error(e); process.exit(1); });
