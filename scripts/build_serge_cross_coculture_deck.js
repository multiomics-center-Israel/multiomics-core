// Build the Serge Ankri cross-co-culture multiomics summary deck via pptxgenjs.
// Plan B (guided tour): framing -> two slides per co-culture -> shared-vs-unique
// synthesis -> close. Consumes the verified aggregation outputs from
//   scripts/serge_cross_coculture_aggregate.R
// plus each run's existing per-co-culture figures.
//
//   node scripts/build_serge_cross_coculture_deck.js
//
// Tone (house rule): understated and neutral; observations not proofs; every data
// slide cites its source file; no bold/colour used for emphasis mid-sentence.

const path = require("path");
const fs = require("fs");
const pptxgen = require("pptxgenjs");

const ROOT = "/home/ozsol/multiomics-core";
const MULTI = path.join(ROOT, "outputs/Serge_multiomics");
const SUMMARY = path.join(MULTI, "cross_coculture_summary");
const SUMFIG = path.join(SUMMARY, "figures");
const SUMTAB = path.join(SUMMARY, "tables");
// Optional CLI arg overrides the output filename (kept in the same summary dir).
const OUT_NAME = process.argv[2] || "Serge_cross_coculture_summary_deck.pptx";
const OUT = path.isAbsolute(OUT_NAME) ? OUT_NAME : path.join(SUMMARY, OUT_NAME);

// --- Neutral academic palette (no emphasis reds) ---
const INK = "1B2430";        // titles
const STEEL = "3E6D9C";      // structural accents
const MUTED = "5D6580";      // secondary text
const BODY = "24303B";       // body text
const PANEL = "F2F5F9";      // light panel fill
const RULE = "C7D0DB";       // hairlines
const WHITE = "FFFFFF";

const HEADER_FONT = "Georgia";
const BODY_FONT = "Calibri";
const SLIDE_W = 13.3;
const SLIDE_H = 7.5;

// --- Small IO helpers ---
function readCSV(file) {
  const txt = fs.readFileSync(file, "utf8").replace(/\r\n?/g, "\n").trimEnd();
  const lines = txt.split("\n");
  const header = splitCSV(lines[0]);
  return lines.slice(1).map((l) => {
    const cells = splitCSV(l);
    const row = {};
    header.forEach((h, i) => (row[h] = cells[i] ?? ""));
    return row;
  });
}
function splitCSV(line) {
  // Minimal CSV split honouring double-quoted fields (pathway names may hold commas).
  const out = [];
  let cur = "", q = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (q) {
      if (c === '"' && line[i + 1] === '"') { cur += '"'; i++; }
      else if (c === '"') q = false;
      else cur += c;
    } else if (c === '"') q = true;
    else if (c === ",") { out.push(cur); cur = ""; }
    else cur += c;
  }
  out.push(cur);
  return out.map((s) => s.replace(/^"|"$/g, ""));
}
function findMultiomicsDir(token) {
  const base = path.join(MULTI, token);
  const results = fs.readdirSync(base)
    .filter((d) => d.startsWith("Results"))
    .map((d) => path.join(base, d))
    .filter((d) => fs.existsSync(path.join(d, "multiomics")))
    .sort((a, b) => fs.statSync(b).mtimeMs - fs.statSync(a).mtimeMs);
  if (!results.length) throw new Error("no multiomics dir for " + token);
  return path.join(results[0], "multiomics");
}

// --- Data ---
const summary = readCSV(path.join(SUMTAB, "per_coculture_summary.csv"));
const summaryBy = Object.fromEntries(summary.map((r) => [r.coculture, r]));
const robustTop = readCSV(path.join(SUMTAB, "robust_top_by_coculture.csv"));
const robustSumm = readCSV(path.join(SUMTAB, "robust_feature_summary.csv"));

// Resolve a KEGG pathview PNG for one co-culture + ehi id (maps live per run).
// Prefer the multi-omics data overlay (.multi_ora.multi.png, differs per run);
// the bare <ehi>.png is the data-less KEGG reference map (identical across runs).
function pathviewMap(token, ehi) {
  const dir = path.join(findMultiomicsDir(token), "multigsea/multi_ora/pathview");
  for (const suffix of [".multi_ora.multi.png", ".png"]) {
    const p = path.join(dir, ehi + suffix);
    if (fs.existsSync(p)) return p;
  }
  return null;
}

// Co-culture order + honesty notes. n_joint / pct_agree come from the summary
// table so the numbers stay correct if the aggregation is re-run.
const COCULTURES = [
  { token: "E.coli", label: "E. coli",
    note: "Modest cross-omics signal; a set of pathways reaches joint significance, roughly half agreeing on direction." },
  { token: "L.rhamnosus", label: "L. rhamnosus",
    note: "The largest set of jointly significant pathways, though fewer than half agree on direction between RNA and protein." },
  { token: "MixSpp", label: "Mixed species",
    note: "No pathway reaches joint significance in RNA and protein; shown as a cross-omics null for this contrast." },
  { token: "E.faecalis", label: "E. faecalis",
    note: "A small set of jointly significant pathways; the proteomics differential signal in this contrast is near-null." },
  { token: "B.subtilis", label: "B. subtilis",
    note: "The proteomics run for this contrast is coverage-driven (fewer proteins identified, heavily imputed, samples separate on missingness), so the protein side is not a reliable co-culture signal." },
];

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.title = "E. histolytica co-culture multiomics — cross-co-culture summary";
pres.author = "MultiOmics Center, Technion";

pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: WHITE },
  objects: [
    { rect: { x: 0, y: 0, w: SLIDE_W, h: 0.16, fill: { color: INK } } },
    { text: {
        text: "E. histolytica co-culture multiomics  ·  cross-co-culture summary  ·  2026-07-22",
        options: { x: 0.4, y: SLIDE_H - 0.34, w: 9.5, h: 0.3,
          fontFace: BODY_FONT, fontSize: 9, color: MUTED, margin: 0 } } },
  ],
  slideNumber: { x: SLIDE_W - 0.7, y: SLIDE_H - 0.34, w: 0.3, h: 0.3,
    fontFace: BODY_FONT, fontSize: 9, color: MUTED, align: "right" },
});

// --- Slide helpers ---
function addTitle(slide, text, sub, kicker) {
  if (kicker) {
    slide.addText(kicker.toUpperCase(), { x: 0.4, y: 0.28, w: SLIDE_W - 0.8, h: 0.25,
      fontFace: BODY_FONT, fontSize: 11, color: STEEL, charSpacing: 2, margin: 0 });
  }
  slide.addText(text, { x: 0.4, y: kicker ? 0.5 : 0.34, w: SLIDE_W - 0.8, h: 0.55,
    fontFace: HEADER_FONT, fontSize: 25, bold: true, color: INK, margin: 0 });
  if (sub) {
    slide.addText(sub, { x: 0.4, y: kicker ? 1.05 : 0.9, w: SLIDE_W - 0.8, h: 0.32,
      fontFace: BODY_FONT, fontSize: 12.5, italic: true, color: MUTED, margin: 0 });
  }
}
function sourceLine(slide, text, y) {
  slide.addText("Source: " + text, { x: 0.4, y: y ?? SLIDE_H - 0.62, w: SLIDE_W - 0.8, h: 0.26,
    fontFace: BODY_FONT, fontSize: 9.5, italic: true, color: MUTED, margin: 0 });
}
// Three-line reader template: question / what we see / caveat.
function readerBlock(slide, x, y, w, q, see, caveat) {
  const runs = [
    { text: "What we asked  ", options: { fontFace: BODY_FONT, fontSize: 12, color: STEEL } },
    { text: q + "\n", options: { fontFace: BODY_FONT, fontSize: 12, color: BODY } },
    { text: "What we see  ", options: { fontFace: BODY_FONT, fontSize: 12, color: STEEL } },
    { text: see + "\n", options: { fontFace: BODY_FONT, fontSize: 12, color: BODY } },
    { text: "Caveat  ", options: { fontFace: BODY_FONT, fontSize: 12, color: STEEL } },
    { text: caveat, options: { fontFace: BODY_FONT, fontSize: 12, color: BODY } },
  ];
  slide.addText(runs, { x, y, w, h: 2.6, fontFace: BODY_FONT, lineSpacingMultiple: 1.15,
    paraSpaceAfter: 6, valign: "top", margin: 0 });
}
function img(slide, file, o) {
  slide.addImage({ path: file, x: o.x, y: o.y, w: o.w, h: o.h,
    sizing: { type: "contain", w: o.w, h: o.h } });
  if (o.cap) {
    slide.addText(o.cap, { x: o.x, y: o.y + o.h + 0.02, w: o.w, h: 0.3,
      fontFace: BODY_FONT, fontSize: 9.5, italic: true, color: MUTED, align: "center", margin: 0 });
  }
}

// =====================================================================
// A. FRAMING
// =====================================================================
// 1. Title
{
  const s = pres.addSlide();
  s.background = { color: INK };
  s.addShape(pres.shapes.RECTANGLE, { x: 0, y: 0, w: 0.22, h: SLIDE_H, fill: { color: STEEL }, line: { width: 0 } });
  s.addText("Entamoeba histolytica in co-culture", { x: 0.9, y: 1.5, w: 11.5, h: 0.9,
    fontFace: HEADER_FONT, fontSize: 40, bold: true, color: WHITE, margin: 0 });
  s.addText("RNA-seq + proteomics, integrated per co-culture", { x: 0.9, y: 2.5, w: 11.5, h: 0.7,
    fontFace: HEADER_FONT, fontSize: 26, color: "CFD8E3", margin: 0 });
  s.addShape(pres.shapes.RECTANGLE, { x: 0.9, y: 3.35, w: 2.3, h: 0.045, fill: { color: STEEL }, line: { width: 0 } });
  s.addText("Shared themes across five co-cultures, and what is specific to each", {
    x: 0.9, y: 3.6, w: 11.5, h: 0.5, fontFace: BODY_FONT, fontSize: 16, color: "CFD8E3", margin: 0 });
  s.addText("MultiOmics Center, Technion  ·  2026-07-22", { x: 0.9, y: 6.4, w: 11.5, h: 0.4,
    fontFace: BODY_FONT, fontSize: 13, color: "9AA6B4", margin: 0 });
}

// 2. Study design
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "Study design", "Each co-culture is compared against E. histolytica grown alone.", "Design");
  const rows = [
    [{ text: "Element", options: { bold: true, color: WHITE, fill: { color: INK }, fontSize: 12 } },
     { text: "Detail", options: { bold: true, color: WHITE, fill: { color: INK }, fontSize: 12 } }],
    ["Control", "E. histolytica trophozoites grown alone"],
    ["Co-cultures (5)", "E. coli · L. rhamnosus · Mixed species · E. faecalis · B. subtilis"],
    ["Contrast", "One two-group contrast per co-culture (co-culture vs control)"],
    ["Replicates", "4 per group (4 vs 4), matched samples"],
    ["Omics", "RNA-seq (DESeq2) and proteomics (limma, DIA-NN)"],
    ["Organism", "E. histolytica, non-model (EHI_ loci for RNA, XP_ groups for protein)"],
  ].map((r, i) => Array.isArray(r) && typeof r[0] === "object" ? r
    : r.map((c) => ({ text: c, options: { fontSize: 12, color: BODY,
        fill: { color: i % 2 ? PANEL : WHITE } } })));
  s.addTable(rows, { x: 0.4, y: 1.5, w: 12.5, colW: [3.2, 9.3], border: { type: "solid", color: RULE, pt: 0.5 },
    rowH: 0.42, valign: "middle", margin: [3, 6, 3, 6] });
  sourceLine(s, "config/multiomics_serge_<coculture>.yaml (design, sample_filter)");
}

// 3. How to read these slides + analysis overview
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "How to read these slides", "One idea per slide; every data slide names its source file.", "Reader guide");
  const gloss = [
    ["Co-culture vs control", "E. histolytica with a bacterium, compared to E. histolytica alone."],
    ["Differential expression", "A gene/protein whose level differs between the two groups."],
    ["MOFA / DIABLO", "Methods that find the axes of variation shared between the RNA and protein data."],
    ["Cross-omics enrichment", "Pathways that stand out when RNA and protein evidence are combined (fgsea)."],
    ["NES / direction", "Enrichment score sign: higher (up) or lower (down) in the co-culture."],
    ["Jointly significant", "A pathway significant in RNA and in protein, and after combining the two p-values."],
  ];
  let y = 1.55;
  for (const [term, def] of gloss) {
    s.addText([
      { text: term + " — ", options: { fontFace: BODY_FONT, fontSize: 13, color: STEEL } },
      { text: def, options: { fontFace: BODY_FONT, fontSize: 13, color: BODY } },
    ], { x: 0.5, y, w: 12.3, h: 0.5, valign: "top", margin: 0 });
    y += 0.62;
  }
  s.addText("Every result here is from n=4 per group and a single contrast per co-culture. Treat all of it as hypothesis-generating.",
    { x: 0.5, y: y + 0.1, w: 12.3, h: 0.5, fontFace: BODY_FONT, fontSize: 12, italic: true, color: MUTED, margin: 0 });
}

// =====================================================================
// B. PER CO-CULTURE (2 slides each)
// =====================================================================
for (const cc of COCULTURES) {
  const dir = findMultiomicsDir(cc.token);
  const sm = summaryBy[cc.token];
  const nJoint = sm ? sm.n_joint : "0";
  const pctAgree = sm ? sm.pct_agree : "n/a";

  // Slide 1 — integration view
  {
    const s = pres.addSlide({ masterName: "CONTENT" });
    addTitle(s, cc.label + " — integration", "How the RNA and protein samples relate under MOFA and DIABLO.", cc.label);
    img(s, path.join(dir, "integration/mofa/mofa_variance_heatmap.png"),
      { x: 0.4, y: 1.55, w: 6.1, h: 4.4, cap: "MOFA variance explained per factor (mofa_variance_heatmap.png)" });
    img(s, path.join(dir, "integration/diablo/diablo_sample_plot.png"),
      { x: 6.8, y: 1.55, w: 6.1, h: 4.4, cap: "DIABLO sample projection (diablo_sample_plot.png)" });
    sourceLine(s, cc.token + " run: integration/mofa/, integration/diablo/");
  }

  // Slide 2 — cross-omics pathways
  {
    const s = pres.addSlide({ masterName: "CONTENT" });
    addTitle(s, cc.label + " — cross-omics pathways",
      "Where RNA and protein enrichment line up for this co-culture.", cc.label);
    img(s, path.join(dir, "cross_enrichment/cross_omics_pathway_heatmap.png"),
      { x: 0.4, y: 1.5, w: 7.4, h: 4.7, cap: "Cross-omics pathway NES heatmap (cross_omics_pathway_heatmap.png)" });
    const seen = cc.token === "MixSpp"
      ? "No pathway is significant in both RNA and protein at the 0.05 threshold."
      : nJoint + " pathways are significant in both RNA and protein; " + pctAgree +
        "% of them move in the same direction, the rest disagree on sign.";
    readerBlock(s, 8.1, 1.6, 4.8,
      "Do RNA and protein point to the same pathways here?",
      seen, cc.note);
    sourceLine(s, cc.token + " run: cross_enrichment/cross_omics_pathways_meta_analysis.csv; per_coculture_summary.csv");
  }
}

// =====================================================================
// C. SYNTHESIS — shared vs unique
// =====================================================================
// 14. Shared themes matrix
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "Shared themes across co-cultures",
    "A cell is shown only where a pathway is jointly significant in that co-culture.", "Shared vs unique");
  img(s, path.join(SUMFIG, "pathway_by_coculture_matrix.png"),
    { x: 0.4, y: 1.5, w: 8.0, h: 4.9 });
  s.addText([
    { text: "Recurring across co-cultures: ", options: { fontFace: BODY_FONT, fontSize: 12.5, color: STEEL } },
    { text: "ribosome biogenesis, nucleolar, ribosome and rRNA-processing terms, which score higher (up) in several co-cultures. Nuclear-import terms and nucleotide-binding/GTPase terms also recur, but there RNA and protein disagree on direction (opposite).",
      options: { fontFace: BODY_FONT, fontSize: 12.5, color: BODY } },
  ], { x: 8.6, y: 1.6, w: 4.3, h: 3.6, valign: "top", lineSpacingMultiple: 1.18, margin: 0 });
  sourceLine(s, "cross_coculture_summary/tables/pathway_by_coculture_matrix.csv (fgsea meta-analysis per run)");
}

// 15. Shared pathways bar (>=2 co-cultures)
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "Pathways shared by two or more co-cultures",
    "Most jointly significant pathways are specific to one co-culture; these are the exceptions.", "Shared vs unique");
  img(s, path.join(SUMFIG, "shared_pathways_bar.png"), { x: 0.4, y: 1.5, w: 8.2, h: 4.9 });
  s.addText("The two pathways shared by three co-cultures are nucleolus and ribosome biogenesis. The remainder are shared by two. Sharing here means jointly significant in each; it does not require the same direction.",
    { x: 8.8, y: 1.7, w: 4.1, h: 3.5, fontFace: BODY_FONT, fontSize: 12.5, color: BODY, valign: "top", lineSpacingMultiple: 1.18, margin: 0 });
  sourceLine(s, "cross_coculture_summary/tables/pathway_by_coculture_matrix.csv (n_cocultures column)");
}

// 16a-b. Pathview of shared KEGG pathways regulated across several co-cultures.
const PATHVIEWS = [
  { ehi: "ehi03008", name: "Ribosome biogenesis in eukaryotes",
    tokens: ["L.rhamnosus", "MixSpp"],
    note: "Jointly significant (RNA + protein) in three co-cultures and scoring higher; a recurring theme rather than a single-co-culture effect." },
  { ehi: "ehi03010", name: "Ribosome",
    tokens: ["E.coli", "L.rhamnosus"],
    note: "Jointly significant in two co-cultures; part of the same ribosome / translation-machinery signal." },
];
for (const pv of PATHVIEWS) {
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, pv.name, "KEGG " + pv.ehi + " — E. histolytica; nodes coloured by fold change.", "Shared pathway");
  const maps = pv.tokens.map((t) => ({ t, p: pathviewMap(t, pv.ehi) })).filter((m) => m.p);
  const w = maps.length >= 2 ? 5.9 : 8.0;
  maps.slice(0, 2).forEach((m, i) => {
    img(s, m.p, { x: 0.4 + i * 6.2, y: 1.5, w, h: 4.5,
      cap: m.t + "  (" + path.basename(m.p) + ")" });
  });
  s.addText(pv.note, { x: 0.4, y: 6.15, w: 12.5, h: 0.5, fontFace: BODY_FONT, fontSize: 12,
    color: BODY, italic: true, margin: 0 });
  sourceLine(s, "per-run KEGG pathview (" + pv.ehi + "); shared status from pathway_by_coculture_matrix.csv", SLIDE_H - 0.4);
}

// 17. Per-co-culture response counts + table
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "Cross-omics response per co-culture",
    "How many pathways reach joint significance, and whether the two omics agree on direction.", "Shared vs unique");
  img(s, path.join(SUMFIG, "per_coculture_counts.png"), { x: 0.4, y: 1.5, w: 7.4, h: 3.4 });
  const head = ["Co-culture", "Joint", "Same up", "Same down", "Opposite", "% agree"];
  const keys = ["coculture", "n_joint", "n_same_up", "n_same_down", "n_opposite", "pct_agree"];
  const trows = [head.map((h) => ({ text: h, options: { bold: true, color: WHITE, fill: { color: INK }, fontSize: 11, align: "center" } }))];
  summary.forEach((r, i) => {
    trows.push(keys.map((k, j) => ({ text: String(r[k]),
      options: { fontSize: 11, color: BODY, align: j === 0 ? "left" : "center", fill: { color: i % 2 ? PANEL : WHITE } } })));
  });
  s.addTable(trows, { x: 0.5, y: 5.2, w: 7.2, colW: [1.7, 0.9, 1.2, 1.3, 1.2, 0.9],
    border: { type: "solid", color: RULE, pt: 0.5 }, rowH: 0.34, valign: "middle" });
  s.addText([
    { text: "MixSpp is not shown — it has no jointly significant pathway. ", options: { fontFace: BODY_FONT, fontSize: 12.5, color: BODY } },
    { text: "Across co-cultures a large share of jointly significant pathways move in opposite directions in RNA and protein, so the joint significance reflects combined evidence, not a single coherent direction.",
      options: { fontFace: BODY_FONT, fontSize: 12.5, color: BODY } },
  ], { x: 8.1, y: 1.7, w: 4.8, h: 4.5, valign: "top", lineSpacingMultiple: 1.18, margin: 0 });
  sourceLine(s, "cross_coculture_summary/tables/per_coculture_summary.csv");
}

// =====================================================================
// C2. MARKERS — top robust integration features (shared vs unique)
// =====================================================================
// 18. Top robust features per co-culture
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "Top robust features per co-culture",
    "Genes ranked robust across the MOFA and DIABLO integrations (consensus).", "Markers");
  const byCC = {};
  robustTop.forEach((r) => { (byCC[r.coculture] = byCC[r.coculture] || []).push(r); });
  const head = [{ text: "Co-culture", options: { bold: true, color: WHITE, fill: { color: INK }, fontSize: 11 } },
                { text: "Total robust", options: { bold: true, color: WHITE, fill: { color: INK }, fontSize: 11, align: "center" } },
                { text: "Top robust features (EHI_ locus tags)", options: { bold: true, color: WHITE, fill: { color: INK }, fontSize: 11 } }];
  const rows = [head];
  COCULTURES.forEach((cc, i) => {
    const list = (byCC[cc.token] || []).slice(0, 8).map((r) => r.feature).join(", ");
    const tot = (robustSumm.find((r) => r.coculture === cc.token) || {}).n_robust_total || "0";
    rows.push([
      { text: cc.label, options: { fontSize: 10.5, color: BODY, fill: { color: i % 2 ? PANEL : WHITE } } },
      { text: String(tot), options: { fontSize: 10.5, color: BODY, align: "center", fill: { color: i % 2 ? PANEL : WHITE } } },
      { text: list || "—", options: { fontSize: 9.5, color: BODY, fill: { color: i % 2 ? PANEL : WHITE } } },
    ]);
  });
  s.addTable(rows, { x: 0.4, y: 1.6, w: 12.5, colW: [2.0, 1.4, 9.1],
    border: { type: "solid", color: RULE, pt: 0.5 }, rowH: 0.55, valign: "middle", margin: [3, 6, 3, 6] });
  s.addText("These are transcriptomics-layer features; their XP_ accessions and full ranking are in the table file. Most are uncharacterised E. histolytica loci not confidently detected at the protein level, so they read best at the pathway level (earlier slides).",
    { x: 0.4, y: 5.3, w: 12.5, h: 0.9, fontFace: BODY_FONT, fontSize: 11.5, color: BODY, italic: true, valign: "top", margin: 0 });
  sourceLine(s, "cross_coculture_summary/tables/robust_top_by_coculture.csv (per-run consensus/robust_features.csv)");
}

// 19. Robust markers are co-culture-specific
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "Top robust features are co-culture-specific",
    "Little to no overlap at the gene level, in contrast to the shared pathway themes.", "Markers");
  img(s, path.join(SUMFIG, "robust_features_specificity.png"), { x: 0.4, y: 1.5, w: 7.6, h: 3.6 });
  s.addText("None of the 53 distinct top features is shared by two or more co-cultures. The shared signal in this study is at the pathway level (ribosome biogenesis, nucleolar and rRNA-processing terms), not at the level of individual genes.",
    { x: 8.2, y: 1.7, w: 4.7, h: 3.5, fontFace: BODY_FONT, fontSize: 12.5, color: BODY, valign: "top", lineSpacingMultiple: 1.18, margin: 0 });
  sourceLine(s, "cross_coculture_summary/tables/robust_feature_overlap.csv; robust_feature_summary.csv");
}

// =====================================================================
// D. CLOSE
// =====================================================================
// 17. Limitations
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "Limitations", null, "Close");
  const pts = [
    "n = 4 per group and a single contrast per co-culture; all findings are hypothesis-generating.",
    "Proteomics differential signal is modest; it is near-null for E. faecalis and absent as a joint signal for the mixed-species contrast.",
    "The B. subtilis proteomics run is coverage-driven (heavy imputation, samples separate on missingness); its protein side is not a reliable co-culture signal.",
    "Joint significance combines RNA and protein p-values; many jointly significant pathways still disagree on direction between the two omics.",
    "Pathway sets are GO and KEGG for E. histolytica; enrichment is fgsea with a Fisher combination across omics.",
  ];
  s.addText(pts.map((t) => ({ text: t, options: { bullet: { code: "2013" }, breakLine: true,
    fontFace: BODY_FONT, fontSize: 14, color: BODY, paraSpaceAfter: 10 } })),
    { x: 0.6, y: 1.7, w: 12.1, h: 4.8, valign: "top", lineSpacingMultiple: 1.15, margin: 0 });
  sourceLine(s, "per-run config_used.yaml; de_summary_counts (proteomics); aggregation tables");
}

// 18. Resources
{
  const s = pres.addSlide({ masterName: "CONTENT" });
  addTitle(s, "Where the full detail lives", null, "Resources");
  const items = [
    ["Per-co-culture reports", "outputs/Serge_multiomics/<coculture>/.../multiomics/report_multiomics.html"],
    ["Shared-vs-unique tables", "outputs/Serge_multiomics/cross_coculture_summary/tables/"],
    ["Shared-vs-unique figures", "outputs/Serge_multiomics/cross_coculture_summary/figures/"],
    ["Aggregation script", "scripts/serge_cross_coculture_aggregate.R"],
    ["This deck", "scripts/build_serge_cross_coculture_deck.js"],
  ];
  let y = 1.7;
  for (const [k, v] of items) {
    s.addText([
      { text: k + "\n", options: { fontFace: BODY_FONT, fontSize: 13.5, color: INK, bold: true } },
      { text: v, options: { fontFace: "Consolas", fontSize: 11, color: MUTED } },
    ], { x: 0.6, y, w: 12.1, h: 0.8, valign: "top", margin: 0 });
    y += 0.92;
  }
}

pres.writeFile({ fileName: OUT }).then((f) => {
  console.log("Wrote " + f);
  console.log("Slides: " + pres.slides.length);
});
