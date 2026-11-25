// Build manuscript/Tables_4-5.docx (editable Word tables) from out/tables45_values.json (scripts/72).
// JavaScript (docx); no number is computed here: every cell string comes from the JSON, so the R check (scripts/74)
// compares the document with its own recomputation. A4 landscape, Arial 10 pt in tables, three-line style
// (rules above and below the header and at the foot), as in Tables 1-3.
// Usage: node scripts/73_tables45_docx.js
const fs = require('fs');
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, BorderStyle, AlignmentType,
        PageOrientation, VerticalAlign } = require('docx');
const V = JSON.parse(fs.readFileSync('out/tables45_values.json', 'utf8'));
const FONT = 'Arial', SZ = 20, SZN = 18;                       // half-points: 10 pt cells, 9 pt notes
const NONE = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' }, RULE = { style: BorderStyle.SINGLE, size: 6, color: '000000' };
const TW = 14678;                                              // text width (DXA): A4 landscape, 0.75 in margins
function runs(text, opt = {}) {                                // supports _x_ for italics (n_eff, P) via a tiny marker syntax
  const out = []; const parts = text.split(/(\{i:[^}]*\}|\{sub:[^}]*\})/);
  for (const p of parts) {
    if (!p) continue;
    let m;
    if ((m = p.match(/^\{i:(.*)\}$/))) out.push(new TextRun({ text: m[1], italics: true, font: FONT, size: opt.size || SZ, bold: opt.bold }));
    else if ((m = p.match(/^\{sub:(.*)\}$/))) out.push(new TextRun({ text: m[1], subScript: true, font: FONT, size: opt.size || SZ, bold: opt.bold, italics: true }));
    else out.push(new TextRun({ text: p, font: FONT, size: opt.size || SZ, bold: opt.bold }));
  }
  return out;
}
function cell(text, w, o = {}) {
  return new TableCell({
    width: { size: w, type: WidthType.DXA }, columnSpan: o.span, rowSpan: o.rowSpan, verticalAlign: o.vAlign || VerticalAlign.TOP,
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    borders: { top: o.top ? RULE : NONE, bottom: o.bottom ? RULE : NONE, left: NONE, right: NONE },
    children: [new Paragraph({ alignment: o.align || AlignmentType.LEFT, children: runs(text, { bold: o.bold }) })],
  });
}
let SUPPRESS_LN = false;                                       // set by build(..., { noLineNumbers: true }) for scripts/84
function title(t) { return new Paragraph({ spacing: { before: 0, after: 120 }, keepNext: true, suppressLineNumbers: SUPPRESS_LN || undefined, children: runs(t, { bold: true, size: 22 }) }); }
function note(t) { return new Paragraph({ spacing: { before: 40, after: 0 }, suppressLineNumbers: SUPPRESS_LN || undefined, children: runs(t, { size: SZN }) }); }

// Column widths are scaled to the text width TWX (scripts/84 places the tables within 25 mm margins); at the default
// width the widths are exactly those of manuscript/Tables_4-5.docx.
function fit(w, tw) { const s = w.reduce((a, b) => a + b, 0); if (s === tw) return w.slice(); const o = w.map(x => Math.floor(x * tw / s)); o[0] += tw - o.reduce((a, b) => a + b, 0); return o; }
function build(TWX = TW, opts = {}) {
  SUPPRESS_LN = !!opts.noLineNumbers;
// ---------------- Table 4 ----------------
const w4 = fit([3500, 950, 1950, 750, 1450, 1450, 1450, 1800, 1378], TWX);
const R = AlignmentType.RIGHT, C = AlignmentType.CENTER;
// Two header rows without row spans (Word does not draw the bottom rule of a vertically merged cell reliably):
// row 1 carries only the spanning label, with its own rule beneath; row 2 carries every column label, ruled beneath.
const h4a = new TableRow({ tableHeader: true, children: [
  cell('', w4[0] + w4[1] + w4[2] + w4[3], { top: true, span: 4 }),
  cell('Prewhitened cross-correlation, {i:r} (Holm {i:P})ᵇ', w4[4] + w4[5] + w4[6], { top: true, bottom: true, span: 3, bold: true, align: C }),
  cell('', w4[7] + w4[8], { top: true, span: 2 }) ] });
const h4b = new TableRow({ tableHeader: true, children: [
  ['Analysis', AlignmentType.LEFT], ['Months', R], ['{i:r} (95% CI)ᵃ', R], ['{i:n}{sub:eff}', R], ['Lag −1', R], ['Lag 0', R], ['Lag +1', R],
  ['Lagged search → cases, coefficient ({i:P})ᶜ', R], ['Lead rule metᵈ', R]].map(([t, al], j) => cell(t, w4[j], { bottom: true, bold: true, align: al, vAlign: VerticalAlign.BOTTOM })) });
const rows4 = V.t4.map((r, i) => {
  const last = i === V.t4.length - 1; const sens = i >= 1 && i <= 7;
  const c = r.cells.slice();
  return new TableRow({ children: c.map((t, j) => cell(j === 0 && sens ? '   ' + t : t, w4[j], { bottom: last, align: j === 0 ? AlignmentType.LEFT : R })) });
});
const sensHead = new TableRow({ children: [cell('Sensitivity analyses', TWX, { span: 9, bold: false })] });
const t4 = new Table({ width: { size: TWX, type: WidthType.DXA }, columnWidths: w4, rows: [h4a, h4b, rows4[0], sensHead, ...rows4.slice(1, 8), rows4[8]] });

// ---------------- Table 5 ----------------
const w5 = opts.w5 || fit([3800, 900, 1050, 1050, 1050, 1050, 2778, 1600, 1400], TWX);
const h5a = new TableRow({ tableHeader: true, children: [
  cell('', w5[0] + w5[1], { top: true, span: 2 }),
  cell('Mean absolute error, cases per monthᵃ', w5[2] + w5[3] + w5[4] + w5[5], { top: true, bottom: true, span: 4, bold: true, align: C }),
  cell('', w5[6] + w5[7] + w5[8], { top: true, span: 3 }) ] });
const h5b = new TableRow({ tableHeader: true, children: [
  ['Schedule and search series', AlignmentType.LEFT], ['Test months', R], ['M0', R], ['M1', R], ['M2', R], ['M3', R], ['M2 − M1, cases (95% CI)ᵇ', R],
  ['Seasons with lower error for M2', R], ['Rule criteriaᶜ', R]].map(([t, al], j) => cell(t, w5[j], { bottom: true, bold: true, align: al, vAlign: VerticalAlign.BOTTOM })) });
const rows5 = []; let g = null;
V.t5.forEach((r, i) => {
  if (r.group !== g) { g = r.group; rows5.push(new TableRow({ children: [cell(g + (r.sched === 'mo' ? ' †' : ''), TWX, { span: 9 })] })); }
  const last = i === V.t5.length - 1;
  rows5.push(new TableRow({ children: r.cells.map((t, j) => cell(j === 0 ? '   ' + t : t, w5[j], { bottom: last, align: j === 0 ? AlignmentType.LEFT : R })) }));
});
const t5 = new Table({ width: { size: TWX, type: WidthType.DXA }, columnWidths: w5, rows: [h5a, h5b, ...rows5] });

const T4 = V.t4[0], D = V.t4.find(r => r.key === 'D');
// The two blocks are exported so that the Results document (scripts/84) places the identical tables in the text.
  const block4 = [
      title('Table 1. Association between search interest and reported cases, and their temporal order: primary and sensitivity analyses'),
      t4,
      note('ᵃ Pearson correlation of the seasonally adjusted series (log scale; calendar-month means removed). The 95% CI spans a Fisher interval based on the effective sample size ({i:n}{sub:eff}) and a 12-month block-bootstrap interval. The bootstrap interval lay inside the Fisher interval in every analysis except the "Fever" control, whose upper limit (0.58) is the bootstrap limit.'),
      note('ᵇ Correlation of the seasonally adjusted series after both were filtered with the autoregressive model fitted to search interest. Positive lags denote search interest preceding cases. Holm adjustment across the seven lags −3 to +3.'),
      note('ᶜ Coefficient of log search interest in month {i:t}−1 in an ordinary least squares model of log cases in month {i:t} with monthly indicators and log cases in month {i:t}−1; Newey–West standard errors (three lags); Wald {i:z} test.'),
      note('ᵈ Decision rule (Supplementary Table S3): lag +1 correlation significant after Holm adjustment and larger than the lag −1 correlation, and lagged search interest a significant predictor of cases. The rule is applied to every row for comparison; only the primary analysis decides Question 2.'),
      note(`Analysis (iv) joins December 2021 to January 2023, so its lagged terms span a 12-month gap; its lagged-search coefficient is negative (${D.fwd_text.replace('P =', '{i:P} =')}). Analysis (vii) is identical to the primary analysis, because the only value reported as "<1" in the primary series is in April 2019.`),
      note('CI, confidence interval.'),
];
const block5 = [
      title('Table 2. Accuracy of current-month case estimates for 2022–2024 and the difference between models with and without current-month search interest'),
      t5,
      note('ᵃ M0, previous month\'s count; M1, monthly indicators and the previous month\'s log count; M2, M1 plus current-month log search interest; M3, M1 plus previous-month log search interest. M0 and M1 do not use search interest, so their errors are the same for every series.'),
      note(`ᵇ Mean over test months of the paired difference in absolute error, M2 minus M1; negative values favour M2. 95% percentile interval from 4,000 bootstrap resamples of test months. With overlapping 3-month blocks (post hoc), the intervals for the "dengue" series were ${V.blocks.ps.text} (once per season) and ${V.blocks.mo.text} (every month).`),
      note('ᶜ Decision rule (Supplementary Table S3): lower mean absolute error for M2 than for M1 in at least two of three seasons, and a 95% CI for the difference that excludes zero. The rule is applied to every row for comparison; only the first row (primary schedule, "dengue") decides Question 3. The leave-one-season-out evaluation (February 2020–December 2024) has no rule.'),
      note('† Analysis added after the primary results were known; not used for any decision.'),
      note('All estimates use final rather than real-time counts, and search interest retrieved in November 2026. CI, confidence interval.'),
];
  SUPPRESS_LN = false;
  return { block4, block5 };
}
const { block4, block5 } = build();
module.exports = { build, block4, block5, runs, cell, title, note, RULE, NONE, FONT, SZ, SZN, TW };
if (require.main === module) {
  const doc = new Document({
    styles: { default: { document: { run: { font: FONT, size: SZ } } } },
    sections: [{
      properties: { page: { size: { width: 11906, height: 16838, orientation: PageOrientation.LANDSCAPE }, margin: { top: 1080, bottom: 1080, left: 1080, right: 1080 } } },
      children: [...block4, new Paragraph({ pageBreakBefore: true, children: [] }), ...block5],
    }],
  });
  Packer.toBuffer(doc).then(b => { fs.writeFileSync('manuscript/Tables_4-5.docx', b); console.log('wrote manuscript/Tables_4-5.docx'); });
}
