// Third, independent implementation of the Supplementary Figure S2 values (JavaScript / Node.js).
//   - reads the preprint .xlsx by unzipping it (jszip) and parsing the worksheet XML directly;
//   - parses the three Google Trends CSVs by hand;
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
// Shares no code with scripts/62_figureS2.py or scripts/63_check_figureS2.R. Event categories are assigned here by
// keyword rules on event_type (Python uses an explicit lookup table; R uses its own table).
// Usage: node scripts/64_s3_independent.js
// Output: out/figS2_values_js.csv, out/figS2_events_js.csv
function parseCsv(text) {                      // RFC 4180: quoted fields, doubled quotes, commas inside quotes
  const rows = []; let row = [], f = '', q = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (q) { if (ch === '"') { if (text[i + 1] === '"') { f += '"'; i++; } else q = false; } else f += ch; }
    else if (ch === '"') q = true; else if (ch === ',') { row.push(f); f = ''; }
    else if (ch === '\n' || ch === '\r') { if (ch === '\r' && text[i + 1] === '\n') i++; row.push(f); rows.push(row); row = []; f = ''; }
    else f += ch;
  }
  if (f !== '' || row.length) { row.push(f); rows.push(row); }
  const hdr = rows.shift(); return rows.filter(r => r.length > 1).map(r => Object.fromEntries(hdr.map((h, i) => [h, r[i]])));
}
function category(t) {
  if (/^no national alert/.test(t) || /\(5 dates\)/.test(t)) return null;
  if (/^official count/.test(t)) return 'count';          // before 'alert': one type reads "(no formal alert)"
  if (/alert|response/.test(t)) return 'alert';
  if (/reporting begins|international/.test(t)) return 'report';
  if (/count/.test(t)) return 'count';
  throw new Error('unknown event type: ' + t);
}
function xOf(d) {
  const [y, m, dd] = d.split('-').map(Number); const nd = new Date(Date.UTC(y, m, 0)).getUTCDate();
  return [y, (m - 1) + (dd - 0.5) / nd];
}
(async () => {
  const cases = await readPreprint(process.env.PREPRINT_XLSX || 'data_raw/paper_finalized_data_RP.xlsx');
  const tr = ['data_raw/gt_dengue_term.csv', 'data_raw/gt_dengue_term_day2_2026-11-22.csv', 'data_raw/gt_dengue_term_day3_2026-11-23.csv'].map(readTrends);
  const s = cases.map((_, i) => (tr[0][i] + tr[1][i] + tr[2][i]) / 3);
  const out = ['year,month,cases_rel,search_rel'];
  for (let y = 2019; y <= 2024; y++) {
    const i0 = (y - 2019) * 12;
    const c = cases.slice(i0, i0 + 12).map((v, m) => (y === 2019 && m < 5) ? NaN : v), ss = s.slice(i0, i0 + 12);
    const cm = Math.max(...c.filter(v => !Number.isNaN(v))), sm = Math.max(...ss);
    for (let m = 0; m < 12; m++) out.push([y, m + 1, Number.isNaN(c[m]) ? 'NA' : (c[m] / cm).toPrecision(17), (ss[m] / sm).toPrecision(17)].join(','));
  }
  fs.writeFileSync('out/figS2_values_js.csv', out.join('\n') + '\n');
  const ev = [];
  for (const r of parseCsv(fs.readFileSync('governance/announcement_events.csv', 'utf8'))) { const c = category(r.event_type); if (c) ev.push([r.date, c]); }
  for (const r of parseCsv(fs.readFileSync('out/06_2024_reporting_lag.csv', 'utf8'))) if (!/-12-31$/.test(r.as_of)) ev.push([r.as_of, 'count']);
  const eo = ['date,year,category,x'];
  for (const [d, c] of ev) { const [y, x] = xOf(d); eo.push([d, y, c, x.toPrecision(17)].join(',')); }
  fs.writeFileSync('out/figS2_events_js.csv', eo.join('\n') + '\n');
  console.log(`JS: ${out.length - 1} monthly rows, ${ev.length} events`);
})();
