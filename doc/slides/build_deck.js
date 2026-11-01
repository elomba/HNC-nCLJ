// Builds component_entropy_results.pptx (6 slides) from the study results.
// Run: cd doc/slides && node build_deck.js
const path = require("path");
const pptxgen = require("pptxgenjs");
const SKILL = process.env.PPTX_SKILL;
const { applyTheme } = require(path.join(SKILL, "scripts/apply_theme.js"));

const ROOT = path.resolve(__dirname, "../..");
const FIG = {
  scaling: path.join(ROOT, "doc/figures/component_entropy_diffusion.png"),
  twopt: path.join(ROOT, "doc/figures/component_2pt_validation_delta_fit.png"),
};

const THEME = {
  name: "Ion Mobility",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1B2433", lt1: "FFFFFF", dk2: "0B4F5C", lt2: "EDF3F5",
    accent1: "0E7C86", accent2: "E3A033", accent3: "C2553D",
    accent4: "4F6D9A", accent5: "7D8B96", accent6: "2F9E73",
    hlink: "0E7C86", folHlink: "4F6D9A",
  },
};
const HEX = THEME.colors;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5 in
pres.author = "E. Lomba";
pres.title = "Excess entropy and divalent-ion mobility in charged LJ mixtures";
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
const C = pres.SchemeColor;

const FOOT = "E. Lomba · IQF-CSIC · EntropicFluids · November 2026";

pres.defineSlideMaster({
  title: "TITLE_DARK",
  background: { color: C.text2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.7, y: 1.2, w: 11.9, h: 2.25,
      fontFace: THEME.headFontFace, fontSize: 40, bold: true, color: C.background1, valign: "top", align: "left", margin: 0 },
      text: "Title" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.7, y: 3.6, w: 11.9, h: 0.7,
      fontSize: 20, color: C.accent2, valign: "top", margin: 0 }, text: "Subtitle" } },
  ],
});

pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: C.background1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.35, w: 12.1, h: 0.85,
      fontFace: THEME.headFontFace, fontSize: 30, bold: true, color: C.text1, valign: "middle", align: "left", margin: 0 },
      text: "Title" } },
    { text: { text: FOOT, options: { x: 0.6, y: 7.05, w: 8, h: 0.3, fontSize: 10, color: C.accent5, margin: 0 } } },
  ],
  slideNumber: { x: 12.3, y: 7.05, w: 0.5, h: 0.3, fontSize: 10, color: C.accent5, align: "right" },
});

// ---------- helpers ----------
function card(slide, { x, y, w, h, fill = C.background2, name }) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.08,
    fill: { color: fill }, line: { color: fill }, objectName: name });
}
function stat(slide, { x, y, w, value, label, valueColor = C.accent1, labelColor = C.text1, name, size = 40 }) {
  slide.addText(value, { x, y, w, h: 0.75, fontFace: THEME.headFontFace, fontSize: size, bold: true,
    color: valueColor, margin: 0, isTextBox: true, objectName: name + " value" });
  slide.addText(label, { x, y: y + 0.78, w, h: 0.75, fontSize: 14, color: labelColor, margin: 0,
    valign: "top", isTextBox: true, objectName: name + " label" });
}
function bullets(items, opts = {}) {
  return items.map((t, i) => {
    const runs = Array.isArray(t) ? t : [{ text: t }];
    return runs.map((r, j) => ({ text: r.text, options: {
      bold: !!r.bold, color: r.color, bullet: j === 0 ? (opts.bullet === false ? false : true) : undefined,
      breakLine: j === runs.length - 1 && i < items.length - 1, paraSpaceAfter: opts.gap || 6, ...(r.options || {}) } }));
  }).flat();
}

// ---------- 1. Title ----------
pres.addSection({ title: "Results" });
let s = pres.addSlide({ masterName: "TITLE_DARK", sectionTitle: "Results" });
s.addText("Excess entropy and divalent-ion mobility in multicomponent charged LJ mixtures", { placeholder: "title" });
s.addText("HNC screening · replicated MD · thermodynamic integration, pair entropy and 2PT", { placeholder: "body" });
const chips = [
  ["+0.30 kB", "excess entropy gained by HNC-selected solvent designs"],
  ["+18–20%", "divalent-ion self-diffusion, both signs"],
  ["0%", "effect of adding labels without changing interactions"],
];
chips.forEach(([v, l], i) => {
  const x = 0.7 + i * 4.05;
  card(s, { x, y: 4.65, w: 3.75, h: 1.95, fill: "0E6170", name: `title chip ${i + 1}` });
  s.addText(v, { x: x + 0.3, y: 4.85, w: 3.2, h: 0.7, fontFace: THEME.headFontFace, fontSize: 34, bold: true,
    color: C.accent2, margin: 0, isTextBox: true, objectName: `title chip ${i + 1} value` });
  s.addText(l, { x: x + 0.3, y: 5.6, w: 3.2, h: 0.85, fontSize: 15, color: C.background1, margin: 0,
    valign: "top", isTextBox: true, objectName: `title chip ${i + 1} label` });
});
s.addNotes("Main message: in these model electrolytes, interaction diversity, not the number of components, raises the excess entropy and the mobility of divalent ions. HNC was used to select candidate mixtures; replicated MD with thermodynamic integration (TI) provides the reference entropy and transport. Adding mechanically identical labels raises only the ideal mixing entropy and has exactly zero effect on diffusion.");

// ---------- 2. Design ----------
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Results" });
s.addText("Design: label count vs interaction diversity", { placeholder: "title" });
const steps = [
  ["1", "HNC screen", "35 integer-charge, composition and dispersion candidates at T* = 2, ρ* = 0.7; 8 selected"],
  ["2", "Replicated MD", "8 mixtures × 3 matched replicas, N = 4000, x(±2) = 0.03; probes matched across mixtures"],
  ["3", "Entropy and transport", "TI over 19 densities (fixed LJ cutoff, Debye–Hückel endpoint); MSD in one common COM frame"],
];
steps.forEach(([n, head, body], i) => {
  const y = 1.45 + i * 1.75;
  s.addShape(pres.shapes.OVAL, { x: 0.6, y: y + 0.05, w: 0.6, h: 0.6, fill: { color: C.accent1 },
    line: { color: C.accent1 }, objectName: `step ${n} circle` });
  s.addText(n, { x: 0.6, y: y + 0.05, w: 0.6, h: 0.6, align: "center", valign: "middle", fontSize: 20, bold: true,
    color: C.background1, margin: 0, isTextBox: true, objectName: `step ${n} number` });
  s.addText(head, { x: 1.4, y, w: 4.3, h: 0.45, fontSize: 20, bold: true, color: C.text1, margin: 0,
    isTextBox: true, objectName: `step ${n} head` });
  s.addText(body, { x: 1.4, y: y + 0.5, w: 4.3, h: 1.05, fontSize: 15, color: C.text1, margin: 0, valign: "top",
    isTextBox: true, objectName: `step ${n} body` });
});
const hdr = (t) => ({ text: t, options: { bold: true, color: HEX.lt1, fill: { color: HEX.dk2 } } });
const rows = [
  [hdr("Code"), hdr("Mixture (all with 3% of each ±2 probe)"), hdr("Role")],
  ["B3", "94% LJ solvent", "Baseline"],
  ["B5", "B3 + 5% of each ±1 ion", "Salt addition"],
  ["U4", "64% solvent + 30% weak solvent, ε = 0.25", "Reference"],
  ["L5", "U4 with 2 identical weak-solvent labels", "Label null"],
  ["M5", "U4 split into ε = 0.0625 / 0.5625", "Moderate split"],
  ["H5", "U4 split into ε = 0.01 / 0.81", "Strong split"],
  ["U6 / H7", "U4 / H5 plus 5% of each ±1 ion", "With salt"],
];
s.addTable(rows, { x: 6.0, y: 1.45, w: 6.75, colW: [1.0, 4.15, 1.6], fontSize: 13, color: HEX.dk1,
  border: { type: "solid", pt: 0.75, color: "C9D6DB" }, fill: { color: HEX.lt1 }, rowH: 0.48,
  valign: "middle", objectName: "composition table" });
s.addText("U4 → M5 → H5 keep charged composition, ionic strength and the mean bare LJ amplitude Σᵢⱼ xᵢxⱼεᵢⱼ = 0.7225 fixed",
  { x: 6.0, y: 5.55, w: 6.75, h: 0.8, fontSize: 14, italic: true, color: C.accent1, margin: 0, valign: "top",
    isTextBox: true, objectName: "matched amplitude note" });
s.addNotes("The controls are the core of the design. L5 is mechanically identical to U4 and only adds a label, so it tests whether the number of components matters by itself. M5 and H5 change the distribution of solvent attractions at fixed mean bare amplitude, fixed ionic strength and fixed charged composition. B3/B5 and U6/H7 repeat the comparison with and without monovalent background. Excess entropy is always relative to the ideal mixture of the same composition, so ideal mixing is excluded.");

// ---------- 3. Transport ----------
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Results" });
s.addText("Entropy speeds up ions, less than the solvent", { placeholder: "title" });
const cats = ["H5/B3", "H7/B5", "H5/U4", "H7/U6", "M5/U4", "B5/B3", "L5/U4"];
s.addChart(pres.charts.BAR, [
  { name: "Mean self-diffusion D_net", labels: cats, values: [27.6, 26.6, 9.5, 8.4, 3.0, 0.1, 0.0] },
  { name: "Divalent probes D_probe", labels: cats, values: [20.1, 17.8, 4.8, 5.2, 1.5, 1.8, 0.0] },
], {
  x: 0.5, y: 1.35, w: 7.6, h: 5.45, barDir: "col", barGrouping: "clustered", barGapWidthPct: 60,
  chartColors: [HEX.accent1, HEX.accent2],
  showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 11, dataLabelColor: HEX.dk1, dataLabelFormatCode: "0.0",
  dataLabelFontFace: "+mn-lt",
  showLegend: true, legendPos: "t", legendFontSize: 13, legendFontFace: "+mn-lt", legendColor: HEX.dk1,
  catAxisLabelColor: HEX.dk1, catAxisLabelFontSize: 13, catAxisLabelFontFace: "+mn-lt",
  valAxisLabelColor: HEX.dk1, valAxisLabelFontSize: 12, valAxisLabelFontFace: "+mn-lt",
  valAxisTitle: "Change in self-diffusion (%)", showValAxisTitle: true, valAxisTitleFontSize: 13,
  valAxisTitleColor: HEX.dk1, valAxisTitleFontFace: "+mn-lt", valAxisMinVal: 0, valAxisMaxVal: 32,
  valGridLine: { color: "DCE4E8", size: 0.75 }, catGridLine: { style: "none" },
  showTitle: true, title: "Paired changes (3 replicas)", titleFontSize: 15, titleColor: HEX.dk1, titleFontFace: "+mn-lt",
  objectName: "diffusion change chart",
});
card(s, { x: 8.5, y: 1.45, w: 4.25, h: 5.25, name: "transport findings card" });
s.addText(bullets([
  [{ text: "Full designs: ", bold: true }, { text: "Δs_ex ≈ +0.30 kB, probes +18–20%; each sign accelerates separately" }],
  [{ text: "Controlled diversity: ", bold: true }, { text: "Δs_ex = +0.103 kB, mean +8–9%, probes ≈ +5% (window-sensitive)" }],
  [{ text: "Labels only (L5/U4): ", bold: true }, { text: "identical diffusion in every replica" }],
  [{ text: "Salt only (B5/B3): ", bold: true }, { text: "Δs_ex slightly negative, no resolved speedup" }],
], { gap: 18 }), { x: 8.8, y: 1.75, w: 3.7, h: 4.75, fontSize: 17, color: C.text1, valign: "top", margin: 0,
  isTextBox: true, objectName: "transport findings" });
s.addNotes("Percent changes are exponentiated mean paired log ratios over three matched replicas. The large design changes are robust. The controlled diversification at matched mean amplitude gives a probe response of about 5 percent, roughly half the mean-fluid response; its significance depends on the MSD fitting window and is unresolved with monovalent salt. Six longer replicas (5000 tau) and an N = 8000 set are running to settle this and the residual B3 charge-sign asymmetry.");

// ---------- 4. Entropy coordinate ----------
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Results" });
s.addText("TI entropy tracks mobility; HNC and s₂ follow", { placeholder: "title" });
s.addImage({ path: FIG.scaling, x: 0.5, y: 1.55, w: 8.4, h: 3.36, altText: "Rosenfeld-reduced mean and probe self-diffusion versus MD excess entropy for the seven distinct mixtures", objectName: "rosenfeld figure" });
s.addText("Descriptive Rosenfeld fits (3–4 mixtures, R² > 0.99): B_net ≈ 0.8, B_probe ≈ 0.6. Controlled local probe slope ≈ 0.5, about half the mean-fluid value. Absolute HNC entropies are 0.06–0.09 kB too low; increments are much better.",
  { x: 0.5, y: 5.2, w: 8.4, h: 1.2, fontSize: 15, color: C.text1, margin: 0, valign: "top", isTextBox: true,
    objectName: "rosenfeld caption" });
const st = [
  ["≤ 7%", "HNC error on the controlled entropy increments"],
  ["≈ 80%", "share of s_ex carried by the MD pair entropy s₂, and of its increments"],
  ["−0.33 kB", "beyond-pair residual s_ex − s₂: composition dependent, not negligible"],
];
st.forEach(([v, l], i) => {
  const y = 1.45 + i * 1.75;
  card(s, { x: 9.3, y, w: 3.45, h: 1.55, name: `entropy stat ${i + 1} card` });
  stat(s, { x: 9.5, y: y + 0.12, w: 3.1, value: v, label: l, name: `entropy stat ${i + 1}`, size: 30 });
});
s.addNotes("Excess entropy from thermodynamic integration is the reference coordinate. Both mean and probe diffusion increase with it, but the probe slope is weaker. HNC is an efficient screening tool: it orders the mixtures correctly and predicts the increments within about 7 percent, although absolute values are too negative. The MD pair entropy reproduces the ordering, but leaves a beyond-pair residual of about a third of a kB per particle.");

// ---------- 5. 2PT ----------
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Results" });
s.addText("2PT when TI is unfeasible: calibrate, then rank", { placeholder: "title" });
s.addImage({ path: FIG.twopt, x: 0.85, y: 1.3, w: 11.6, h: 3.48, altText: "Revised 2PT excess entropy versus TI for δ = 1, 1.25 and the fitted 1.20", objectName: "2pt validation figure" });
const tp = [
  ["Calibrate one δ", "Fit over 7 mixtures: δ* = 1.20 [1.17, 1.23]; δ = 1.25 cuts RMS bias 0.067 → 0.029 kB; held-out bias 0.030 kB", C.accent6],
  ["Good within a family", "M5/U4 and H5/U4 increments reproduced within 0.01 kB; label split exactly zero", C.accent1],
  ["Limits", "H7/U6 +30%, baseline → design +17–23%; D enters the fluidicity, so 2PT cannot test an entropy–diffusion link on its own", C.accent3],
];
tp.forEach(([h, b, col], i) => {
  const x = 0.6 + i * 4.1;
  card(s, { x, y: 4.95, w: 3.85, h: 1.75, name: `2pt card ${i + 1}` });
  s.addText(h, { x: x + 0.25, y: 5.05, w: 3.4, h: 0.42, fontSize: 17, bold: true, color: col, margin: 0,
    isTextBox: true, objectName: `2pt card ${i + 1} head` });
  s.addText(b, { x: x + 0.25, y: 5.5, w: 3.4, h: 1.3, fontSize: 14, color: C.text1, margin: 0, valign: "top",
    isTextBox: true, objectName: `2pt card ${i + 1} body` });
});
s.addNotes("2PT needs only one short dense trajectory instead of a 19-density isotherm. The partition exponent delta is phenomenological; a single value calibrated over all mixtures fixes the overall level, giving absolute excess entropies within about 0.03 kB. The bias curves are nearly parallel in delta, so the composition-dependent part of the error does not go away: comparisons across solvent families or with salt are overestimated by 17 to 30 percent. Recommended use: calibrate delta on two or three TI anchor states, then rank candidates. Because diffusion enters the fluidicity, a 2PT entropy-diffusion correlation is not independent evidence.");

// ---------- 6. Outlook ----------
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Results" });
s.addText("Next: fast entropy estimators for electrolytes", { placeholder: "title" });
const cols = [
  ["2PT, velocity based", "Ready", C.accent6, [
    "Calibrated δ ≈ 1.2–1.25 for this model family",
    "One dense trajectory per candidate",
    "Use for ranking; recalibrate on TI anchors",
    "Not independent of diffusion",
  ]],
  ["Voronoi cluster entropy", "Open", C.accent2, [
    "24 trajectories packaged with TI, s₂ and D references",
    "Must reproduce TI order, the 0.03 and 0.10 kB increments and the label null",
    "Structure based, so independent of transport",
  ]],
];
cols.forEach(([h, tag, col, items], i) => {
  const x = 0.6 + i * 6.15;
  card(s, { x, y: 1.4, w: 5.95, h: 2.65, name: `outlook ${i + 1} card` });
  s.addText(h, { x: x + 0.3, y: 1.55, w: 4.3, h: 0.5, fontSize: 18, bold: true, color: C.text1, margin: 0,
    isTextBox: true, objectName: `outlook ${i + 1} head` });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: x + 4.65, y: 1.58, w: 1.0, h: 0.42, rectRadius: 0.1,
    fill: { color: col }, line: { color: col }, objectName: `outlook ${i + 1} tag` });
  s.addText(tag, { x: x + 4.65, y: 1.58, w: 1.0, h: 0.42, align: "center", valign: "middle", fontSize: 13, bold: true,
    color: C.background1, margin: 0, isTextBox: true, objectName: `outlook ${i + 1} tag text` });
  s.addText(bullets(items), { x: x + 0.3, y: 2.2, w: 5.35, h: 1.75, fontSize: 15, color: C.text1, margin: 0,
    valign: "top", isTextBox: true, objectName: `outlook ${i + 1} bullets` });
});
const road = [
  ["HNC screen", "candidate mixtures"],
  ["TI anchors", "2–3 states to calibrate and validate"],
  ["2PT / Voronoi", "rank many candidates"],
  ["Conductivity", "Green–Kubo current, not only self-diffusion"],
];
road.forEach(([h, b], i) => {
  const x = 0.6 + i * 3.1;
  s.addShape(pres.shapes.OVAL, { x, y: 4.5, w: 0.5, h: 0.5, fill: { color: C.accent1 }, line: { color: C.accent1 },
    objectName: `road ${i + 1} dot` });
  s.addText(String(i + 1), { x, y: 4.5, w: 0.5, h: 0.5, align: "center", valign: "middle", fontSize: 15, bold: true,
    color: C.background1, margin: 0, isTextBox: true, objectName: `road ${i + 1} number` });
  s.addText(h, { x: x + 0.65, y: 4.48, w: 2.3, h: 0.4, fontSize: 16, bold: true, color: C.text1, margin: 0,
    isTextBox: true, objectName: `road ${i + 1} head` });
  s.addText(b, { x: x + 0.65, y: 4.9, w: 2.3, h: 0.9, fontSize: 14, color: C.text1, margin: 0, valign: "top",
    isTextBox: true, objectName: `road ${i + 1} body` });
});
s.addText("In progress: 6 longer replicas (5000 τ) and N = 8000 runs of B3, U4, H5, U6, H7 to settle the ≈5% controlled ion effect",
  { x: 0.6, y: 6.15, w: 12.1, h: 0.45, fontSize: 15, italic: true, color: C.accent1, margin: 0, isTextBox: true,
    objectName: "in progress note" });
s.addNotes("For complex electrolytes, TI over a full isotherm for every candidate is too expensive. A practical route: screen with HNC, compute TI only for a few anchor states, and use a cheaper estimator to rank the rest. 2PT is usable now once delta is calibrated, but it borrows information from diffusion. The Voronoi cluster entropy is open: the data package is ready, and it must reproduce the TI ordering, the controlled increments and the label null. Because it uses structure only, it would give an entropy coordinate independent of transport. Finally, conductivity enhancement must be measured from the charge-current autocorrelation; self-diffusion and Nernst-Einstein estimates miss ion correlations.");

(async () => {
  const out = path.join(__dirname, "component_entropy_results.pptx");
  await pres.writeFile({ fileName: out });
  await applyTheme(out, THEME);
  console.log("Wrote", out);
})();
