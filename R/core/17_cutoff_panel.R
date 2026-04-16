# ============================================================================
# 17_cutoff_panel.R — Interactive cutoff controls for HTML reports
#
# Provides an inline control panel at the top of the DE section that lets
# users adjust log2FC and p-value cutoffs in the browser, AND switch between
# raw p-value and adjusted p-value for significance.  Volcano plots, MA plots,
# and DE summary counts update live via Plotly.react().
# ============================================================================

#' Emit the inline cutoff control panel (HTML + CSS + JS engine)
#'
#' Place this at the beginning of the Differential Expression section.
#'
#' @param default_lfc   Numeric; default log2FC cutoff (e.g. 1)
#' @param default_p     Numeric; default p-value cutoff (e.g. 0.05)
#' @param default_ptype Character; "pval" or "padj" — which p-value type is
#'                      selected by default
#' @return Raw HTML string (use cat(...) in an asis chunk)
cutoff_panel_html <- function(default_lfc = 1, default_p = 0.05,
                               default_ptype = "padj") {
  lfc_r  <- round(default_lfc, 4)
  p_r    <- signif(default_p, 4)
  ptype  <- match.arg(default_ptype, c("pval", "padj"))
  pval_checked <- if (ptype == "pval") "checked" else ""
  padj_checked <- if (ptype == "padj") "checked" else ""

  paste0(
'<!-- ===== Cutoff Control Panel (inline, top of DE section) ===== -->
<style>
#cutoff-panel {
  margin: 16px 0 24px 0; padding: 16px 22px;
  border: 1px solid #d0d7de; border-radius: 8px;
  background: #f6f8fa;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  font-size: 14px;
}
#cutoff-panel .cp-title {
  margin: 0 0 12px 0; font-size: 15px; font-weight: 600; color: #24292f;
}
#cutoff-panel .cp-grid {
  display: flex; flex-wrap: wrap; gap: 20px; align-items: flex-end;
}
#cutoff-panel .cp-field { display: flex; flex-direction: column; gap: 4px; }
#cutoff-panel .cp-field > label {
  font-size: 12px; font-weight: 600; color: #57606a; text-transform: uppercase;
  letter-spacing: 0.3px;
}
#cutoff-panel .cp-radio-group {
  display: flex; gap: 14px; padding-top: 2px;
}
#cutoff-panel .cp-radio-group label {
  font-size: 13px; font-weight: 400; color: #24292f; cursor: pointer;
  display: flex; align-items: center; gap: 5px;
}
#cutoff-panel .cp-radio-group input[type="radio"] {
  accent-color: #0969da; margin: 0;
}
#cutoff-panel input[type="number"] {
  width: 100px; padding: 5px 8px; border: 1px solid #d0d7de; border-radius: 5px;
  background: #fff; font-size: 13px; color: #24292f;
}
#cutoff-panel input[type="number"]:focus {
  outline: none; border-color: #0969da; box-shadow: 0 0 0 3px rgba(9,105,218,0.15);
}
#cutoff-panel .cp-actions {
  display: flex; align-items: center; gap: 16px; margin-top: 14px;
  padding-top: 10px; border-top: 1px solid #d8dee4;
}
#cutoff-panel button.cp-reset {
  padding: 5px 14px; border: 1px solid #d0d7de; border-radius: 5px;
  background: #fff; color: #24292f; cursor: pointer; font-size: 12px;
  font-weight: 500;
}
#cutoff-panel button.cp-reset:hover { background: #f3f4f6; }
#cutoff-panel .cp-summary {
  font-size: 13px; color: #57606a;
}
#cutoff-panel .cp-summary .de-up   { color: #cf222e; font-weight: 600; }
#cutoff-panel .cp-summary .de-down { color: #0969da; font-weight: 600; }
#cutoff-panel .cp-summary .de-total { font-weight: 600; color: #24292f; }
</style>

<div id="cutoff-panel">
<div class="cp-title">DE Cutoff Controls</div>
<div class="cp-grid">
<div class="cp-field">
<label>Significance metric</label>
<div class="cp-radio-group">
<label><input type="radio" name="cp-ptype" value="pval" ', pval_checked, '> P-value</label>
<label><input type="radio" name="cp-ptype" value="padj" ', padj_checked, '> Adjusted P-value</label>
</div>
</div>
<div class="cp-field">
<label>P-value threshold</label>
<input type="number" id="cp-pval-num" min="0.0001" max="1" step="0.005" value="', p_r, '">
</div>
<div class="cp-field">
<label>log2FC threshold</label>
<input type="number" id="cp-lfc-num" min="0" max="10" step="0.1" value="', lfc_r, '">
</div>
</div>
<div class="cp-actions">
<button class="cp-reset" onclick="window.__cutoffReset()">Reset to defaults</button>
<span class="cp-summary" id="cp-global-summary"></span>
</div>
</div>

<script>
(function() {
  "use strict";

  // --- Registry ---
  window.__cutoffPlots = window.__cutoffPlots || {};

  var DEFAULT_LFC   = ', lfc_r, ';
  var DEFAULT_PVAL  = ', p_r, ';
  var DEFAULT_PTYPE = "', ptype, '";
  var currentLfc    = DEFAULT_LFC;
  var currentPval   = DEFAULT_PVAL;
  var currentPtype  = DEFAULT_PTYPE;  // "pval" or "padj"
  var debounceTimer = null;

  // --- DOM refs ---
  var lfcNum, pvalNum;

  function initControls() {
    lfcNum   = document.getElementById("cp-lfc-num");
    pvalNum  = document.getElementById("cp-pval-num");
    if (!lfcNum) return;

    lfcNum.addEventListener("input", function() {
      currentLfc = parseFloat(this.value) || 0;
      debouncedUpdate();
    });
    pvalNum.addEventListener("input", function() {
      currentPval = parseFloat(this.value) || 0.05;
      debouncedUpdate();
    });

    // Radio buttons for pval/padj
    var radios = document.querySelectorAll("input[name=\'cp-ptype\']");
    for (var r = 0; r < radios.length; r++) {
      radios[r].addEventListener("change", function() {
        currentPtype = this.value;
        debouncedUpdate();
      });
    }
  }

  function debouncedUpdate() {
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(updateAllPlots, 100);
  }

  // --- Get the p-value for a point based on current ptype ---
  function getPval(data, idx) {
    if (currentPtype === "padj" && data.adjPval && data.adjPval[idx] !== null) {
      return data.adjPval[idx];
    }
    // Fall back to raw pval
    return data.pval[idx];
  }

  // --- Classify a point ---
  // When passFlag is available (imputation-consensus pipelines), a point must
  // also have passFlag===1 to be called significant.  This keeps the
  // interactive counts consistent with pipeline results at default thresholds.
  function classify(logFC, pval, lfcCut, pCut, passFlag) {
    if (passFlag !== undefined && passFlag !== null && passFlag !== 1) return "NS";
    if (pval !== null && !isNaN(pval) && pval <= pCut &&
        Math.abs(logFC) >= lfcCut) {
      return logFC > 0 ? "Up" : "Down";
    }
    return "NS";
  }

  // --- Rebuild a single plot ---
  function rebuildPlot(plotId) {
    var reg = window.__cutoffPlots[plotId];
    if (!reg) return { up: 0, down: 0 };

    var el = document.getElementById(plotId);
    if (!el) return { up: 0, down: 0 };
    // Wrapper div fallback (metabolomics wraps plotly in a container div)
    if (!el.data && !el._fullLayout) {
      var child = el.querySelector(".js-plotly-plot");
      if (child) el = child;
    }
    if (!el || typeof Plotly === "undefined") return { up: 0, down: 0 };

    var data  = reg.data;
    var type  = reg.plotType;     // "volcano" or "ma"
    var label = reg.entityLabel;  // "Protein", "Gene", "Metabolite"
    var ycap  = reg.ycap || null;

    // Separate points into 3 groups
    var groups = { NS: [], Down: [], Up: [] };
    for (var i = 0; i < data.name.length; i++) {
      var lfc_i = data.logFC[i];
      if (lfc_i === null || isNaN(lfc_i)) continue;
      var pval_i = getPval(data, i);
      var pf_i = (data.passFlag && data.passFlag[i] !== undefined) ? data.passFlag[i] : undefined;
      var dir = classify(lfc_i, pval_i, currentLfc, currentPval, pf_i);
      groups[dir].push(i);
    }

    var colors = { NS: "rgba(179,179,179,0.6)", Down: "rgba(0,0,255,0.6)", Up: "rgba(255,0,0,0.6)" };
    var traces = [];

    // Determine y-axis p-value label
    var pLabel = currentPtype === "padj" ? "adj.P.Val" : "P.Value";

    ["NS", "Down", "Up"].forEach(function(dir) {
      var idx = groups[dir];
      var x = [], y = [], text = [];
      for (var j = 0; j < idx.length; j++) {
        var ii = idx[j];
        var lfc_j  = data.logFC[ii];
        var pval_j = getPval(data, ii);
        var nm     = data.name[ii] || "";

        if (type === "volcano") {
          x.push(lfc_j);
          var nlogp = -Math.log10(Math.max(pval_j || 1, 1e-300));
          if (ycap !== null) nlogp = Math.min(nlogp, ycap);
          y.push(nlogp);
        } else {
          // MA plot: x = avgExpr, y = logFC
          x.push(data.avgExpr ? data.avgExpr[ii] : 0);
          y.push(lfc_j);
        }

        // Build hover text
        var hoverParts = [label + ": " + nm,
                          "log2FC: " + lfc_j.toFixed(3)];
        // Show both p-values in hover when available
        if (data.pval[ii] !== null) {
          hoverParts.push("P.Value: " + data.pval[ii].toExponential(3));
        }
        if (data.adjPval && data.adjPval[ii] !== null) {
          hoverParts.push("adj.P.Val: " + data.adjPval[ii].toExponential(3));
        }
        if (type === "ma" && data.avgExpr) {
          hoverParts.push((data.avgExprLabel || "AvgExpr") + ": " + (data.avgExpr[ii] || 0).toFixed(3));
        }
        text.push(hoverParts.join("<br>"));
      }

      var legendLabel = dir + " (" + idx.length + ")";
      traces.push({
        x: x, y: y, text: text,
        type: "scatter", mode: "markers",
        marker: { color: colors[dir], size: 5 },
        name: legendLabel,
        hoverinfo: "text",
        showlegend: true
      });
    });

    // Build threshold shapes
    var shapes = [];
    if (type === "volcano") {
      shapes.push(
        { type: "line", x0: -currentLfc, x1: -currentLfc, y0: 0, y1: 1, yref: "paper",
          line: { dash: "dash", color: "rgba(0,0,0,0.4)", width: 1 } },
        { type: "line", x0: currentLfc, x1: currentLfc, y0: 0, y1: 1, yref: "paper",
          line: { dash: "dash", color: "rgba(0,0,0,0.4)", width: 1 } },
        { type: "line", x0: 0, x1: 1, xref: "paper",
          y0: -Math.log10(Math.max(currentPval, 1e-300)),
          y1: -Math.log10(Math.max(currentPval, 1e-300)),
          line: { dash: "dash", color: "rgba(0,0,0,0.4)", width: 1 } }
      );
    } else {
      // MA: horizontal lines at +/- lfc, plus y=0
      shapes.push(
        { type: "line", x0: 0, x1: 1, xref: "paper", y0: -currentLfc, y1: -currentLfc,
          line: { dash: "dash", color: "rgba(0,0,0,0.3)", width: 1 } },
        { type: "line", x0: 0, x1: 1, xref: "paper", y0: currentLfc, y1: currentLfc,
          line: { dash: "dash", color: "rgba(0,0,0,0.3)", width: 1 } },
        { type: "line", x0: 0, x1: 1, xref: "paper", y0: 0, y1: 0,
          line: { color: "rgba(0,0,0,0.2)", width: 1 } }
      );
    }

    // Preserve existing annotations (from gene search/highlight)
    var existingAnnotations = [];
    try {
      if (el._fullLayout && el._fullLayout.annotations) {
        existingAnnotations = el._fullLayout.annotations.map(function(a) {
          return Object.assign({}, a);
        });
      }
    } catch(e) {}

    var layout = Object.assign({}, el.layout || {});
    layout.shapes = shapes;
    layout.annotations = existingAnnotations;
    // Update y-axis label for volcano
    if (type === "volcano") {
      layout.yaxis = Object.assign({}, layout.yaxis || {});
      layout.yaxis.title = "-log10(" + pLabel + ")";
    }

    Plotly.react(el, traces, layout);

    return { up: groups.Up.length, down: groups.Down.length };
  }

  // --- Update all registered plots + summary spans ---
  function updateAllPlots() {
    var totalUp = 0, totalDown = 0;
    var contrastCounts = {};

    for (var plotId in window.__cutoffPlots) {
      var reg = window.__cutoffPlots[plotId];
      var counts = rebuildPlot(plotId);

      var cname = reg.contrast || "";
      if (cname) {
        if (!contrastCounts[cname]) contrastCounts[cname] = { up: 0, down: 0 };
        // Only count volcano plots to avoid double-counting (MA has same data)
        if (reg.plotType === "volcano") {
          contrastCounts[cname].up   += counts.up;
          contrastCounts[cname].down += counts.down;
          totalUp   += counts.up;
          totalDown += counts.down;
        }
      }
    }

    // Update global summary line
    var summEl = document.getElementById("cp-global-summary");
    if (summEl) {
      var ptypeLabel = currentPtype === "padj" ? "adj. p" : "p";
      summEl.innerHTML =
        "<span class=\\"de-up\\">" + totalUp + " Up</span> &middot; " +
        "<span class=\\"de-down\\">" + totalDown + " Down</span> &middot; " +
        "<span class=\\"de-total\\">" + (totalUp + totalDown) + " Total</span>" +
        " <span style=\\"color:#8b949e;\\">(|log2FC| &ge; " + currentLfc.toFixed(2) +
        ", " + ptypeLabel + " &le; " + (currentPval < 0.01 ? currentPval.toExponential(2) : currentPval.toFixed(4)) + ")</span>";
    }

    // Update per-contrast DE summary spans
    var spans = document.querySelectorAll(".cutoff-de-count");
    for (var s = 0; s < spans.length; s++) {
      var span = spans[s];
      var contrast = span.getAttribute("data-contrast");
      if (contrast && contrastCounts[contrast]) {
        var cc = contrastCounts[contrast];
        var upEl = span.querySelector(".de-up-count");
        var dnEl = span.querySelector(".de-down-count");
        var ttEl = span.querySelector(".de-total-count");
        if (upEl) upEl.textContent = cc.up;
        if (dnEl) dnEl.textContent = cc.down;
        if (ttEl) ttEl.textContent = cc.up + cc.down;
      }
    }
  }

  // --- Reset ---
  window.__cutoffReset = function() {
    currentLfc   = DEFAULT_LFC;
    currentPval  = DEFAULT_PVAL;
    currentPtype = DEFAULT_PTYPE;
    if (lfcNum) {
      lfcNum.value  = currentLfc;
      pvalNum.value = currentPval < 0.01 ? currentPval : currentPval;
    }
    // Reset radio buttons
    var radios = document.querySelectorAll("input[name=\'cp-ptype\']");
    for (var r = 0; r < radios.length; r++) {
      radios[r].checked = (radios[r].value === DEFAULT_PTYPE);
    }
    updateAllPlots();
  };

  // --- Init on DOM ready ---
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initControls);
  } else {
    initControls();
  }
})();
</script>
')
}


#' Register a plot for interactive cutoff updates
#'
#' Emits a <script> tag that stores the plot's raw point data as JSON into
#' window.__cutoffPlots[plot_id], enabling the cutoff panel JS to rebuild
#' traces when thresholds change.
#'
#' @param plot_id       Character; the HTML element ID of the plotly div
#' @param point_data_df Data frame with columns: name, logFC, pval;
#'                      optionally: adjPval, avgExpr, avgExprLabel
#' @param plot_type     "volcano" or "ma"
#' @param entity_label  "Protein", "Gene", or "Metabolite"
#' @param contrast      Character; contrast name (for linking to summary spans)
#' @param ycap          Numeric or NULL; y-axis cap for volcano plots
#' @return Raw HTML string (use cat(...) in asis chunk)
cutoff_register_plot <- function(plot_id, point_data_df, plot_type = "volcano",
                                  entity_label = "Protein", contrast = "",
                                  ycap = NULL) {
  df <- point_data_df
  js_name    <- jsonlite::toJSON(as.character(df$name), auto_unbox = FALSE)
  js_logFC   <- jsonlite::toJSON(as.numeric(df$logFC), auto_unbox = FALSE, na = "null")
  js_pval    <- jsonlite::toJSON(as.numeric(df$pval), auto_unbox = FALSE, na = "null")

  js_avgExpr <- "null"
  js_avgExprLabel <- "null"
  if ("avgExpr" %in% names(df)) {
    js_avgExpr <- jsonlite::toJSON(as.numeric(df$avgExpr), auto_unbox = FALSE, na = "null")
    js_avgExprLabel <- jsonlite::toJSON(
      if ("avgExprLabel" %in% names(df)) df$avgExprLabel[1] else "AvgExpr",
      auto_unbox = TRUE
    )
  }

  js_adjPval <- "null"
  if ("adjPval" %in% names(df)) {
    js_adjPval <- jsonlite::toJSON(as.numeric(df$adjPval), auto_unbox = FALSE, na = "null")
  }

  # Optional imputation-consensus pass flag (proteomics: pass.imputs column).
  # When present, the cutoff panel ANDs this flag with the slider thresholds
  # so that interactive counts stay consistent with pipeline results.
  js_passFlag <- "null"
  if ("passFlag" %in% names(df)) {
    # Convert to 1/0/null for JS
    pf <- df$passFlag
    pf_num <- ifelse(is.na(pf), NA_real_, ifelse(pf == 1 | pf == TRUE, 1, 0))
    js_passFlag <- jsonlite::toJSON(as.numeric(pf_num), auto_unbox = FALSE, na = "null")
  }

  js_ycap <- if (is.null(ycap)) "null" else as.character(ycap)

  sprintf('
<script>
(function() {
  window.__cutoffPlots = window.__cutoffPlots || {};
  window.__cutoffPlots["%s"] = {
    plotType: "%s",
    entityLabel: "%s",
    contrast: "%s",
    ycap: %s,
    data: {
      name: %s,
      logFC: %s,
      pval: %s,
      avgExpr: %s,
      avgExprLabel: %s,
      adjPval: %s,
      passFlag: %s
    }
  };
})();
</script>',
    plot_id, plot_type, entity_label, contrast, js_ycap,
    js_name, js_logFC, js_pval, js_avgExpr, js_avgExprLabel, js_adjPval, js_passFlag)
}


#' Emit a DE summary count span that the cutoff panel JS can update
#'
#' @param contrast_name  Character; the contrast identifier
#' @param n_up           Integer; initial up-regulated count
#' @param n_down         Integer; initial down-regulated count
#' @param entity         Character; "proteins", "genes", or "metabolites"
#' @return Raw HTML string
cutoff_summary_span <- function(contrast_name, n_up, n_down, entity = "features") {
  total <- n_up + n_down
  sprintf(
    paste0('<span class="cutoff-de-count" data-contrast="%s">',
           '<span class="de-up-count" style="color:#c0392b;font-weight:600;">%d</span> up / ',
           '<span class="de-down-count" style="color:#2980b9;font-weight:600;">%d</span> down / ',
           '<span class="de-total-count" style="font-weight:600;">%d</span> total %s',
           '</span>'),
    htmltools::htmlEscape(contrast_name), n_up, n_down, total, entity
  )
}
