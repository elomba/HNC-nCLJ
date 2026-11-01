// Builds component_entropy_results.pptx (6 slides) from the study results.
// Run: cd doc/slides && node build_deck.js
const path = require("path");
const fs = require("fs");
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

// ---------- 2b. Particle species ----------
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Results" });
s.addText("Seven kinds of particles, one size", { placeholder: "title" });
const SP = {
  A: { fill: "8FA3AD", line: "8FA3AD", txt: "" },
  B: { fill: "C7DCE1", line: "9DB7BE", txt: "" },
  B1: { fill: "FFFFFF", line: "9DB7BE", txt: "", dash: "dash" },
  B2: { fill: "3F7480", line: "3F7480", txt: "" },
  P2: { fill: HEX.accent3, line: HEX.accent3, txt: "2+" },
  M2: { fill: HEX.accent4, line: HEX.accent4, txt: "2−" },
  P1: { fill: "E8A08E", line: "E8A08E", txt: "+" },
  M1: { fill: "A9BCDB", line: "A9BCDB", txt: "−" },
};
function particle(slide, key, x, y, d, name) {
  const p = SP[key];
  slide.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: p.fill },
    line: { color: p.line, width: 1.25, dashType: p.dash || "solid" }, objectName: name });
  if (p.txt) slide.addText(p.txt, { x, y, w: d, h: d, align: "center", valign: "middle", margin: 0,
    fontSize: d > 0.4 ? 14 : 10, bold: true, color: HEX.lt1, isTextBox: true, objectName: name + " label" });
}
const legend = [
  [["A"], "Solvent A", "q = 0, ε = 1", "all mixtures, 54–94%"],
  [["B"], "Weak solvent B", "q = 0, ε = 0.25", "U4, L5, U6: 30%"],
  [["B1", "B2"], "Split solvent B₁ / B₂", "q = 0, ε = 0.01 / 0.81 (H)\nor 0.0625 / 0.5625 (M)", "M5, H5, H7: 15% each"],
  [["P2", "M2"], "Divalent probe ions", "q = ±2, ε = 1", "all mixtures: 3% each"],
  [["P1", "M1"], "Monovalent background", "q = ±1, ε = 1", "B5, U6, H7: 5% each"],
];
legend.forEach(([keys, name, prop, where], i) => {
  const y = 1.45 + i * 0.98;
  keys.forEach((k, j) => particle(s, k, 0.6 + j * 0.62, y + 0.12, 0.52, `legend ${i + 1} particle ${j + 1}`));
  s.addText(name, { x: 1.95, y: y, w: 2.9, h: 0.4, fontSize: 16, bold: true, color: C.text1, margin: 0,
    isTextBox: true, objectName: `legend ${i + 1} name` });
  s.addText(prop, { x: 1.95, y: y + 0.38, w: 2.9, h: 0.5, fontSize: 13, color: C.text1, margin: 0, valign: "top",
    isTextBox: true, objectName: `legend ${i + 1} properties` });
  s.addText(where, { x: 4.85, y: y + 0.05, w: 2.3, h: 0.7, fontSize: 13, italic: true, color: C.accent1, margin: 0,
    valign: "top", isTextBox: true, objectName: `legend ${i + 1} where` });
});
s.addText("All particles: σ = 1, m = 1;  εᵢⱼ = √(εᵢεⱼ);  Coulomb with Bjerrum length 0.5σ at T* = 2",
  { x: 0.6, y: 6.4, w: 6.6, h: 0.5, fontSize: 13, color: C.text1, margin: 0, valign: "top", isTextBox: true,
    objectName: "species footnote" });
// cartoon snapshot of H7 (2D, schematic)
const BX = 7.65, BY = 1.45, BW = 5.1, BH = 4.85, D = 0.4;
s.addShape(pres.shapes.RECTANGLE, { x: BX, y: BY, w: BW, h: BH, fill: { color: HEX.lt2 },
  line: { color: "B9C9CF", width: 1 }, objectName: "snapshot box" });
let seed = 7;
const rand = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
const counts = [["P2", 3], ["M2", 3], ["P1", 4], ["M1", 4], ["B1", 12], ["B2", 12], ["A", 43]];
const placed = [];
counts.forEach(([k, n]) => {
  for (let c = 0; c < n; c++) {
    for (let tries = 0; tries < 5000; tries++) {
      const x = BX + 0.08 + rand() * (BW - D - 0.16), y = BY + 0.08 + rand() * (BH - D - 0.16);
      if (placed.every(([px, py]) => (px - x) ** 2 + (py - y) ** 2 > (D * 1.04) ** 2)) {
        placed.push([x, y]); particle(s, k, x, y, D, `snapshot ${k} ${c + 1}`); break;
      }
    }
  }
});
s.addText("Cartoon of H7 (2D sketch with roughly its proportions; not simulation data)", { x: BX, y: 6.4, w: BW, h: 0.5,
  fontSize: 12, italic: true, color: C.accent5, margin: 0, valign: "top", isTextBox: true, objectName: "snapshot caption" });
s.addNotes("All species have the same size and mass, so the mixtures differ only in charge and in the strength of the LJ attraction epsilon. The divalent probes are the particles whose mobility we follow. The weak solvent B lowers the average attraction; splitting it into a very weak B1 and a strong B2 keeps the mean amplitude but makes the interactions more diverse. The open dashed circle is the almost non-attractive B1. The picture on the right is a 2D cartoon with roughly H7 proportions, not a simulation snapshot.");

// ---------- 2c. HNC primer ----------
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Results" });
s.addText("HNC for charged mixtures in four steps", { placeholder: "title" });
const hnc = [
  ["Pair structure", "gᵢⱼ(r),  hᵢⱼ = gᵢⱼ − 1",
    "Probability of finding a j particle at distance r from an i particle, relative to a random mixture"],
  ["Ornstein–Zernike equation", "hᵢⱼ = cᵢⱼ + Σₖ ρₖ cᵢₖ ⊗ hₖⱼ",
    "Total correlation = direct part + indirect part relayed through every species k (a p × p matrix equation in k space)"],
  ["HNC closure (Morita–Hiroike)", "gᵢⱼ = exp(−βuᵢⱼ + hᵢⱼ − cᵢⱼ)",
    "Second relation that closes the system; it neglects the bridge diagrams and is accurate for ions"],
  ["Charges (Allnatt; Lomba–Høye)", [{ text: "cᵢⱼ → −λ" }, { text: "B", options: { subscript: true } }, { text: " qᵢqⱼ / r" }],
    "Coulomb tails decay too slowly to Fourier transform: subtract the screened (Debye) tail analytically and solve for short-ranged parts"],
];
hnc.forEach(([h, eq, txt], i) => {
  const y = 1.4 + i * 1.33;
  s.addShape(pres.shapes.OVAL, { x: 0.6, y: y + 0.03, w: 0.5, h: 0.5, fill: { color: C.accent1 }, line: { color: C.accent1 },
    objectName: `hnc step ${i + 1} dot` });
  s.addText(String(i + 1), { x: 0.6, y: y + 0.03, w: 0.5, h: 0.5, align: "center", valign: "middle", fontSize: 15, bold: true,
    color: C.background1, margin: 0, isTextBox: true, objectName: `hnc step ${i + 1} number` });
  s.addText(h, { x: 1.3, y, w: 2.9, h: 0.55, fontSize: 16, bold: true, color: C.text1, margin: 0, valign: "middle",
    isTextBox: true, objectName: `hnc step ${i + 1} head` });
  s.addText(eq, { x: 4.25, y, w: 3.3, h: 0.55, fontFace: "Cambria", fontSize: 17, color: C.accent3, margin: 0,
    valign: "middle", isTextBox: true, objectName: `hnc step ${i + 1} equation` });
  s.addText(txt, { x: 1.3, y: y + 0.58, w: 6.25, h: 0.65, fontSize: 14, color: C.text1, margin: 0, valign: "top",
    isTextBox: true, objectName: `hnc step ${i + 1} text` });
});
// real HNC g(r) for B3 from the solver output
const rdf = fs.readFileSync(path.join(ROOT, "ComponentEntropyStudy_20261030/hnc/fine/baseline_3C/RDFhncr.dat"), "utf8")
  .trim().split("\n").map((l) => l.trim().split(/\s+/).map(Number)).filter((r) => r[0] >= 0.8 && r[0] <= 3.5);
const sub = rdf.filter((_, k) => k % 4 === 0);
const xs = sub.map((r) => r[0]);
s.addChart(pres.charts.SCATTER, [
  { name: "r", values: xs },
  { name: "2+ … solvent", values: sub.map((r) => r[2]) },
  { name: "2+ … 2+", values: sub.map((r) => r[3]) },
  { name: "2+ … 2−", values: sub.map((r) => r[5]) },
], {
  x: 7.85, y: 1.3, w: 4.95, h: 3.9, lineSize: 2, lineDataSymbol: "none",
  chartColors: [HEX.accent1, HEX.accent3, HEX.accent4],
  showLegend: true, legendPos: "t", legendFontSize: 12, legendFontFace: "+mn-lt", legendColor: HEX.dk1,
  showTitle: true, title: "HNC g(r) in B3 (T* = 2, ρ* = 0.7)", titleFontSize: 14, titleColor: HEX.dk1, titleFontFace: "+mn-lt",
  valAxisMinVal: 0, valAxisMaxVal: 5.5, catAxisMinVal: 0.5, catAxisMaxVal: 3.5, catAxisMajorUnit: 0.5, valAxisMajorUnit: 1,
  catAxisLabelColor: HEX.dk1, catAxisLabelFontSize: 11, catAxisLabelFontFace: "+mn-lt",
  valAxisLabelColor: HEX.dk1, valAxisLabelFontSize: 11, valAxisLabelFontFace: "+mn-lt",
  showCatAxisTitle: true, catAxisTitle: "r / σ", catAxisTitleFontSize: 12, catAxisTitleColor: HEX.dk1, catAxisTitleFontFace: "+mn-lt",
  showValAxisTitle: true, valAxisTitle: "g(r)", valAxisTitleFontSize: 12, valAxisTitleColor: HEX.dk1, valAxisTitleFontFace: "+mn-lt",
  valGridLine: { color: "DCE4E8", size: 0.75 }, catGridLine: { style: "none" }, objectName: "hnc rdf chart",
});
card(s, { x: 7.85, y: 5.35, w: 4.95, h: 1.45, name: "hnc thermo card" });
s.addText([
  { text: "Why it is useful here: ", options: { bold: true } },
  { text: "in HNC the free energy follows in closed form from g and c, so sₑₓ = uₑₓ − aₑₓ costs one solution per mixture, enough to screen dozens of compositions" },
], { x: 8.05, y: 5.45, w: 4.55, h: 1.25, fontSize: 14, color: C.text1, margin: 0, valign: "top", isTextBox: true,
  objectName: "hnc thermo text" });
s.addNotes("For students: g(r) is the quantity measured by scattering and computed by simulation. The Ornstein-Zernike equation is exact but introduces a second unknown, the direct correlation function c; a closure relation is needed. The hypernetted-chain closure of Morita and Hiroike sums a large class of diagrams and works well for ionic systems. For charged particles, c behaves like the Coulomb potential at large distances, which cannot be Fourier transformed numerically; following Allnatt, the screened long-range part is treated analytically and only short-ranged functions are discretized, as implemented in the Lomba-Hoye algorithm. The chart is a real HNC result for the B3 mixture: oppositely charged probes form strong contact pairs (peak near 5), like-charged probes avoid each other, and the solvent structure is that of a dense LJ liquid. Finally, HNC gives the excess free energy in closed form, so the excess entropy of each candidate mixture comes from a single solution.");

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
