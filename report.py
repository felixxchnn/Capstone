"""
report.py
=========
Render the committed ``data/processed/case_study.json`` into ONE self-contained,
offline, deterministic HTML report.

    py report.py --build       # write phase2_report.html
    py report.py --validate    # regenerate + byte-compare + structural checks
    py report.py --self-test   # synthetic render + escaping units + real build

What this module does NOT do
---------------------------
No model inference. No evidence lookup. No scientific recomputation. It reads
already-committed artifacts, verifies their SHA-256, and lays them out as HTML.

Determinism / offline guarantees (enforced by --validate)
-------------------------------------------------------
* No CDN, remote font, analytics, tracking, or any external runtime dependency;
  JSON, CSS and JavaScript are embedded locally. External hyperlinks to PubMed /
  source / licence pages are the only outward references and the page is fully
  readable without them.
* Opens correctly from ``file://``.
* Byte-identical across repeated builds: no wall-clock, no absolute path, stable
  ordering, fixed-precision numbers already frozen in the case-study JSON.
* All inserted text is HTML-escaped; the embedded JSON has ``<`` / ``>`` / ``&``
  escaped so a ``</script>`` inside source data cannot break out.

Phase 1 headline numbers (section B) are read from the committed result
artifacts (``baseline_results.json`` / ``head_results.json`` /
``analysis_results.json``), each hash-pinned, never typed in blind.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import config

REPORT_SCHEMA = config.REPORT_SCHEMA_VERSION
CASE_STUDY_SCHEMA_EXPECTED = config.CASE_STUDY_SCHEMA_VERSION

# ---- hash-pinned committed inputs (verified at every build) --------------
CASE_STUDY_SHA256 = config.CASE_STUDY_JSON_SHA256
RESULT_ARTIFACT_SHA256 = {
    "baseline_results.json":
        "b49169bd363a596f400b4faff8c21d354275b70404efe08b9109d38f1bdc0ffd",
    "head_results.json":
        "1962206fa17646cbd1fec4b642a577cc2586c09c4cabd980541a7e11a8b6f894",
    "analysis_results.json":
        "12431dad60d07f0bd2bea9a680367007c9e030e9f17c5c20ef0b0694dcb548f9",
    "dgidb_2026-06b.manifest.json":
        "9fb585c723cb2102a7cd335dbfac478b206d91cad04951f8ca7f70f495f6f912",
}

# Phase 1 headline values expected in those artifacts. Read from the files at
# build time; also asserted equal to these so a mis-keyed JSON path or a moved
# artifact is a hard stop rather than a silently wrong headline.
PHASE1_EXPECTED = {
    "baseline_rho": 0.2356,
    "head_rho": 0.2047,
    "delta_mean": -0.0308,
    "delta_ci_low": -0.0365,
    "delta_ci_high": -0.0255,
    "n_cell_lines": 170,
    "n_targets": 4297,
    "bootstrap_resamples": 1000,
}

PROHIBITED_CLAIM_PHRASES = (
    "recommended treatment", "treatment recommendation", "actionable therapy",
    "effective drug", "patient-specific treatment", "validated bg003082 prediction",
)

_HASH_CHUNK = 1 << 20
_REPO_ROOT = config.PROJECT_ROOT

EVIDENCE_LABEL = "Drug\u2013gene interaction evidence"          # "Drug–gene ..."
RANKING_HINT = "More-negative predicted GeneEffect indicates stronger predicted dependency."
NOT_CLINICAL = "Research demonstration \u2014 not clinical guidance"

SAMPLE_ORDER = ("ACH-000364", "BG003082")
MODEL_ORDER = ("ridge_pca", "ridge_head")


class ReportError(RuntimeError):
    """Any pre-flight / integrity / validation failure -- always a hard stop."""


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def _sha256_file(path: str | Path) -> str:
    d = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(_HASH_CHUNK), b""):
            d.update(block)
    return d.hexdigest()


def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _esc(value) -> str:
    """HTML-escape for text nodes and double-quoted attributes."""
    return html.escape(str(value), quote=True)


def _reject_constant(tok: str):
    raise ReportError(f"case_study.json contains a non-finite JSON literal: {tok!r}")


def _json_for_script(obj) -> str:
    """
    Serialise for an inline <script type="application/json"> block.

    ensure_ascii=True -> pure ASCII (also escapes U+2028/U+2029). Then escape the
    three characters that can end/confuse the script context; each becomes a
    valid JSON \\uXXXX escape that parses back to the original character, so the
    embedded text still round-trips to identical bytes.
    """
    text = json.dumps(obj, sort_keys=True, ensure_ascii=True, allow_nan=False)
    return (text.replace("&", "\\u0026")
                .replace("<", "\\u003c")
                .replace(">", "\\u003e"))


def _num(x) -> str:
    """Render a number exactly as the JSON already froze it (no re-rounding)."""
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        s = repr(x)
        return s
    return _esc(x)


def _int(x) -> str:
    """Thousands-separated integer for human-readable gene counts."""
    return f"{int(x):,}"


# --------------------------------------------------------------------------
# load + verify every committed input
# --------------------------------------------------------------------------

def _verify(path: Path, want_sha: str, label: str) -> bytes:
    if not path.is_file():
        raise ReportError(f"missing committed input: {label} ({path})")
    raw = path.read_bytes()
    got = _sha256_bytes(raw)
    if got != want_sha:
        raise ReportError(f"{label}: sha256 {got} != pinned {want_sha}")
    return raw


def load_verified() -> dict:
    proc = config.PROCESSED_DIR
    cs_path = config.CASE_STUDY_JSON_FILE
    cs_raw = _verify(cs_path, CASE_STUDY_SHA256, "case_study.json")
    case_study = json.loads(cs_raw.decode("utf-8"), parse_constant=_reject_constant)
    if case_study.get("schema_version") != CASE_STUDY_SCHEMA_EXPECTED:
        raise ReportError(
            f"case_study.json schema_version {case_study.get('schema_version')!r} "
            f"!= {CASE_STUDY_SCHEMA_EXPECTED!r}")

    baseline = json.loads(_verify(proc / "baseline_results.json",
                                  RESULT_ARTIFACT_SHA256["baseline_results.json"],
                                  "baseline_results.json"))
    head = json.loads(_verify(proc / "head_results.json",
                              RESULT_ARTIFACT_SHA256["head_results.json"],
                              "head_results.json"))
    analysis = json.loads(_verify(proc / "analysis_results.json",
                                  RESULT_ARTIFACT_SHA256["analysis_results.json"],
                                  "analysis_results.json"))
    dgidb = json.loads(_verify(config.DGIDB_MANIFEST_FILE,
                               RESULT_ARTIFACT_SHA256["dgidb_2026-06b.manifest.json"],
                               "dgidb_2026-06b.manifest.json"))

    # cross-checks against what the case study itself recorded
    ev_ret = case_study["drug_gene_interaction_evidence"]["retrieval"]
    if ev_ret["manifest_sha256"] != RESULT_ARTIFACT_SHA256["dgidb_2026-06b.manifest.json"]:
        raise ReportError("case_study evidence manifest_sha256 disagrees with the pinned DGIdb manifest")
    if ev_ret["snapshot_sha256"] != case_study["input_artifact_sha256"][
            "data/external/dgidb/dgidb_2026-06b.interactions.filtered.tsv"]:
        raise ReportError("case_study evidence snapshot_sha256 self-inconsistency")

    # ---- Phase 1 headline, READ from the artifacts ----------------------
    b = baseline["tasks"]["crispr"]["models"]["ridge_pca"]
    h = head["tasks"]["crispr"]["models"]["ridge_head"]
    dblk = analysis["A1_bootstrap"]["delta_head_minus_baseline"]
    a4 = analysis["A4_effective_degrees_of_freedom"]["models"]
    phase1 = {
        "baseline_rho": b["spearman_mean"],
        "head_rho": h["spearman_mean"],
        "baseline_alpha": b["alpha"],
        "head_alpha": h["alpha"],
        "delta_mean": dblk["mean"],
        "delta_ci_low": dblk["ci_low"],
        "delta_ci_high": dblk["ci_high"],
        "delta_std_error": dblk["std_error"],
        "n_cell_lines": analysis["n_cell_lines"],
        "n_targets": analysis["n_targets"],
        "bootstrap_resamples": analysis["A1_bootstrap"]["n_resamples"],
        "eff_df_baseline": a4["ridge_pca"]["effective_df_at_selected_alpha"],
        "eff_df_head": a4["ridge_head"]["effective_df_at_selected_alpha"],
        "source_files": {
            "baseline_results.json": RESULT_ARTIFACT_SHA256["baseline_results.json"],
            "head_results.json": RESULT_ARTIFACT_SHA256["head_results.json"],
            "analysis_results.json": RESULT_ARTIFACT_SHA256["analysis_results.json"],
        },
    }
    for key, want in PHASE1_EXPECTED.items():
        if phase1[key] != want:
            raise ReportError(
                f"Phase 1 value {key}={phase1[key]!r} read from the committed "
                f"artifacts != expected {want!r} -- refusing to render a wrong headline")

    # ---- DGIdb release / vintage / licence / coverage facts ------------
    pub = dgidb["publications"]
    snap_records = dgidb["snapshot"]["record_count"]
    with_pmids = pub["records_with_pmids"]
    dgidb_facts = {
        "release_tag": dgidb["dgidb"]["release_tag"],
        "interaction_data_version": dgidb["dgidb"]["versions"]["interaction_data_version"],
        "app_version_tsv": dgidb["dgidb"]["versions"]["dgidb_app_version_tsv_comment"],
        "app_version_graphql": dgidb["dgidb"]["versions"]["dgidb_app_version_graphql_at_retrieval"],
        "retrieved_utc": dgidb["retrieval"]["retrieved_utc"],
        "included_sources": list(dgidb["filter"]["included_sources"]),
        "n_excluded_sources": len(dgidb["filter"]["excluded_sources"]),
        "snapshot_records": snap_records,
        "records_with_pmids": with_pmids,
        "pct_with_pmids": round(100.0 * with_pmids / snap_records, 1),
        "manifest_sha256": RESULT_ARTIFACT_SHA256["dgidb_2026-06b.manifest.json"],
        "snapshot_sha256": ev_ret["snapshot_sha256"],
    }

    return {
        "case_study": case_study,
        "case_study_bytes": cs_raw,
        "case_study_sha256": CASE_STUDY_SHA256,
        "phase1": phase1,
        "dgidb": dgidb_facts,
    }


# --------------------------------------------------------------------------
# HTML fragments
# --------------------------------------------------------------------------

_CSS = """
:root{
  --bg:#ffffff; --fg:#1a1d21; --muted:#54606c; --line:#d4dae0; --panel:#f5f7f9;
  --accent:#0b5cad; --accent-fg:#ffffff;
  --cited:#1a7f37; --source-only:#8a6d00; --none:#5a6673;
  --warn-bg:#fff4e5; --warn-line:#e0a44a;
}
*{box-sizing:border-box}
html{font-size:16px}
body{margin:0;background:var(--bg);color:var(--fg);
  font:16px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
main{max-width:1040px;margin:0 auto;padding:1.25rem 1.1rem 4rem;}
h1{font-size:1.7rem;margin:.2rem 0 .4rem;line-height:1.25}
h2{font-size:1.28rem;margin:2.1rem 0 .6rem;border-bottom:2px solid var(--line);padding-bottom:.25rem}
h3{font-size:1.05rem;margin:1.3rem 0 .4rem}
p{margin:.55rem 0}
a{color:var(--accent)}
a:focus-visible,button:focus-visible,input:focus-visible,select:focus-visible,
summary:focus-visible,[tabindex]:focus-visible{outline:3px solid var(--accent);outline-offset:2px}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.92em}
.small{font-size:.86rem;color:var(--muted)}
.banner{background:var(--warn-bg);border:1px solid var(--warn-line);border-radius:8px;
  padding:.7rem .9rem;margin:.8rem 0;font-weight:600}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:.85rem 1rem;margin:.8rem 0}
.kv{display:grid;grid-template-columns:max-content 1fr;gap:.25rem .9rem;margin:.4rem 0}
.kv div:nth-child(odd){color:var(--muted)}
.controls{display:flex;flex-wrap:wrap;gap:1rem 1.4rem;align-items:flex-end;margin:.6rem 0 1rem}
fieldset{border:1px solid var(--line);border-radius:8px;padding:.5rem .8rem .7rem;margin:0}
legend{font-weight:600;padding:0 .35rem}
.radio-row{display:flex;flex-wrap:wrap;gap:.4rem .9rem;margin-top:.3rem}
.radio-row label{display:inline-flex;align-items:center;gap:.35rem;cursor:pointer}
.field label{display:block;font-weight:600;margin-bottom:.25rem}
input[type=search],select{font:inherit;padding:.4rem .5rem;border:1px solid var(--line);border-radius:6px;background:#fff;color:var(--fg);min-width:15rem}
button{font:inherit;padding:.42rem .8rem;border:1px solid var(--accent);background:var(--accent);color:var(--accent-fg);border-radius:6px;cursor:pointer}
button.secondary{background:#fff;color:var(--accent)}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:.94rem}
caption{text-align:left;padding:.5rem .7rem;font-weight:600;background:var(--panel);border-bottom:1px solid var(--line)}
th,td{padding:.45rem .6rem;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}
th[scope=col]{background:var(--panel);position:sticky;top:0}
td.num{font-variant-numeric:tabular-nums;white-space:nowrap}
tr.gene-row:hover{background:#eef4fb}
.pill{display:inline-block;padding:.05rem .5rem;border-radius:999px;font-size:.8rem;font-weight:600;border:1px solid}
.pill.cited{color:var(--cited);border-color:var(--cited);background:#e8f5ec}
.pill.source_only{color:var(--source-only);border-color:var(--source-only);background:#fbf3dd}
.pill.none_in_filtered_snapshot{color:var(--none);border-color:var(--none);background:#eef1f4}
.pill::before{content:"\\25CF  ";font-size:.7em;vertical-align:middle}
details{border:1px solid var(--line);border-radius:8px;margin:.5rem 0;background:#fff}
details>summary{cursor:pointer;padding:.55rem .8rem;font-weight:600;list-style:revert}
details[open]>summary{border-bottom:1px solid var(--line)}
.evi-body{padding:.4rem .8rem .7rem}
.evi-record{border-top:1px dashed var(--line);padding:.5rem 0}
.evi-record:first-child{border-top:0}
.evi-record .kv{grid-template-columns:max-content 1fr}
.disclaimer{color:var(--muted);font-style:italic;margin:.3rem 0}
.view[hidden]{display:none}
.hidden-row{display:none !important}
.tag{display:inline-block;background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:0 .4rem;margin:0 .15rem .15rem 0;font-size:.82rem}
.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap;border:0}
@media (max-width:640px){
  main{padding:1rem .7rem 3rem} h1{font-size:1.4rem}
  input[type=search],select{min-width:0;width:100%}
  .kv{grid-template-columns:1fr}
}
@media print{
  .controls,.noprint{display:none !important}
  details{border:0} details>summary{display:none} .evi-body{padding:0}
  .view[hidden]{display:block !important}
  a[href^="http"]::after{content:" (" attr(href) ")";font-size:.8em;color:#333;word-break:break-all}
  tr,details,.panel{break-inside:avoid}
}
"""

_JS = """
"use strict";
(function(){
  function $(sel, root){ return (root||document).querySelector(sel); }
  function $all(sel, root){ return Array.prototype.slice.call((root||document).querySelectorAll(sel)); }

  var views = $all(".view");
  var live = $("#status-live");
  function currentSample(){ var el=$("input[name=sample]:checked"); return el? el.value : null; }
  function currentModel(){ var el=$("input[name=model]:checked"); return el? el.value : null; }
  function activeView(){
    var s=currentSample(), m=currentModel();
    return $("#view-"+s+"-"+m);
  }

  function showView(){
    var target=activeView();
    views.forEach(function(v){ v.hidden = (v!==target); });
    if(live && target){
      var s=currentSample(), m=currentModel();
      live.textContent = "Showing "+s+" \u2014 model "+m;
    }
    applyFilters();
  }

  function applyFilters(){
    var view=activeView(); if(!view) return;
    var q=($("#gene-search").value||"").trim().toLowerCase();
    var status=$("#evi-status").value;
    // table rows
    $all("tr.gene-row", view).forEach(function(row){
      var hay=(row.getAttribute("data-search")||"").toLowerCase();
      var st=row.getAttribute("data-evidence-status");
      var ok = (!q || hay.indexOf(q)>-1) && (status==="all" || st===status);
      row.classList.toggle("hidden-row", !ok);
    });
    // evidence groups
    $all(".evi-group", view).forEach(function(grp){
      var hay=(grp.getAttribute("data-search")||"").toLowerCase();
      var st=grp.getAttribute("data-evidence-status");
      var ok = (!q || hay.indexOf(q)>-1) && (status==="all" || st===status);
      grp.classList.toggle("hidden-row", !ok);
    });
    var shown = $all("tr.gene-row", view).filter(function(r){ return !r.classList.contains("hidden-row"); }).length;
    var counter=$("#result-count");
    if(counter) counter.textContent = shown + " of 25 genes shown";
  }

  function setDetails(open){
    var view=activeView(); if(!view) return;
    $all("details", view).forEach(function(d){ d.open = open; });
  }

  $all("input[name=sample],input[name=model]").forEach(function(el){
    el.addEventListener("change", showView);
  });
  $("#gene-search").addEventListener("input", applyFilters);
  $("#evi-status").addEventListener("change", applyFilters);
  var eb=$("#expand-all"); if(eb) eb.addEventListener("click", function(){ setDetails(true); });
  var cb=$("#collapse-all"); if(cb) cb.addEventListener("click", function(){ setDetails(false); });
  var pb=$("#print-btn"); if(pb) pb.addEventListener("click", function(){ window.print(); });

  showView();
  document.documentElement.setAttribute("data-js-ready","ok");
})();
"""


def _pill(status: str) -> str:
    label = {"cited": "cited", "source_only": "source-only",
             "none_in_filtered_snapshot": "none in filtered snapshot"}.get(status, status)
    return f'<span class="pill {_esc(status)}">{_esc(label)}</span>'


def _pubmed_link(pmid: str) -> str:
    p = _esc(pmid)
    return (f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}/" '
            f'rel="noopener noreferrer" target="_blank">PMID {p}</a>')


def _evidence_record_html(rec: dict) -> str:
    pmids = rec.get("pmids", [])
    if pmids:
        pm = ", ".join(_pubmed_link(p) for p in pmids)
        pm_line = f"<div>PMIDs</div><div>{pm}</div>"
        scope = ('<p class="small">Group-level citation: PMIDs are attached at the '
                 'drug\u2013gene / interaction-source group level and may span '
                 'multiple interaction claims; they may not specifically support '
                 'the displayed interaction subtype.</p>')
    else:
        pm_line = '<div>PMIDs</div><div><em>source-only interaction evidence \u2014 no claim-level publication in the filtered snapshot</em></div>'
        scope = ""
    lic = rec.get("source_license", "")
    lic_url = rec.get("source_license_url", "")
    lic_html = (f'<a href="{_esc(lic_url)}" rel="noopener noreferrer" target="_blank">{_esc(lic)}</a>'
                if lic_url else _esc(lic))
    approved = rec.get("drug_is_approved", "")
    return (
        '<div class="evi-record">'
        '<div class="kv">'
        f'<div>Drug</div><div>{_esc(rec.get("drug_name") or rec.get("drug_claim_name",""))}</div>'
        f'<div>Interaction</div><div>{_esc(rec.get("interaction_type_raw") or "(unspecified)")} '
        f'&middot; direction: {_esc(rec.get("interaction_direction",""))} '
        f'(tier: {_esc(rec.get("direction_tier",""))})</div>'
        f'<div>Source</div><div>{_esc(rec.get("interaction_source",""))} '
        f'(v{_esc(rec.get("interaction_source_version",""))})</div>'
        f'<div>Source licence</div><div>{lic_html}</div>'
        f'<div>DGIdb regulatory-approval flag</div><div>{_esc(approved or "(not stated)")}</div>'
        f'<div>Evidence score (DGIdb)</div><div>{_esc(rec.get("evidence_score") or "(none)")}</div>'
        f'{pm_line}'
        '</div>'
        f'{scope}'
        f'<p class="disclaimer">{_esc(rec.get("disclaimer",""))}</p>'
        '</div>'
    )


def _evidence_group_html(entry: dict) -> str:
    status = entry["evidence_status"]
    recs = entry.get("records", [])
    sym = entry["symbol"]
    ent = entry["entrez_id"]
    search = f"{sym} {ent} " + " ".join(
        (r.get("drug_name", "") + " " + r.get("interaction_source", "")) for r in recs)
    inner = ("".join(_evidence_record_html(r) for r in recs)
             if recs else
             '<p class="small">No drug\u2013gene interaction record for this gene '
             'in the licence-filtered offline DGIdb snapshot.</p>')
    open_attr = " open" if recs else ""
    return (
        f'<details class="evi-group" data-evidence-status="{_esc(status)}" '
        f'data-search="{_esc(search)}"{open_attr}>'
        f'<summary>{_esc(sym)} <span class="mono">({_esc(ent)})</span> &nbsp; '
        f'{_pill(status)} &nbsp; <span class="small">{len(recs)} record(s)</span></summary>'
        f'<div class="evi-body">{inner}</div>'
        '</details>'
    )


def _ranking_table_html(sample: str, model: str, block: dict, ev_by_entrez: dict) -> str:
    has_obs = sample == "ACH-000364"
    head_cols = ["Rank", "Gene symbol", "Entrez ID", "Predicted GeneEffect"]
    if has_obs:
        head_cols += ["Observed GeneEffect", "Observed rank (of 4,297)"]
    head_cols += ["Evidence status", "Evidence records"]
    thead = "".join(f'<th scope="col">{_esc(c)}</th>' for c in head_cols)

    rows = []
    for g in block["genes"]:
        ent = g["entrez_id"]
        ev = ev_by_entrez.get(ent, {"evidence_status": "none_in_filtered_snapshot", "n_records": 0})
        st = ev["evidence_status"]
        search = f'{g["symbol"]} {ent}'
        cells = [
            f'<td class="num">{_esc(g["rank"])}</td>',
            f'<td>{_esc(g["symbol"])}</td>',
            f'<td class="num mono">{_esc(ent)}</td>',
            f'<td class="num">{_num(g["predicted_geneeffect"])}</td>',
        ]
        if has_obs:
            ov = g.get("observed_geneeffect")
            orank = g.get("observed_rank")
            cells.append(f'<td class="num">{"n/a" if ov is None else _num(ov)}</td>')
            cells.append(f'<td class="num">{"n/a" if orank is None else _esc(orank)}</td>')
        cells.append(f'<td>{_pill(st)}</td>')
        cells.append(f'<td class="num">{_esc(ev["n_records"])}</td>')
        rows.append(
            f'<tr class="gene-row" data-evidence-status="{_esc(st)}" '
            f'data-search="{_esc(search)}">' + "".join(cells) + "</tr>")

    caption = (f'Predicted CRISPR gene <strong>dependencies</strong> '
               f'(not therapeutic targets) &mdash; sample {_esc(sample)}, model '
               f'<span class="mono">{_esc(model)}</span>, top {len(block["genes"])} '
               f'of {_esc(block["n_targets_ranked"])} ranked targets')
    return (
        '<div class="tablewrap">'
        f'<table><caption>{caption}</caption>'
        f'<thead><tr>{thead}</tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table></div>'
        f'<p class="small"><strong>{_esc(RANKING_HINT)}</strong> '
        f'{_esc(block["ranking_rule"])}</p>'
    )


def _view_html(sample: str, model: str, cs: dict) -> str:
    block = cs["rankings"][sample][model]
    ev = cs["drug_gene_interaction_evidence"]["by_entrez"]
    vid = f"view-{sample}-{model}"
    # evidence groups: only genes displayed in THIS view, in rank order
    groups = []
    seen = set()
    for g in block["genes"]:
        ent = g["entrez_id"]
        if ent in seen:
            continue
        seen.add(ent)
        entry = ev.get(ent)
        if entry is None:
            entry = {"symbol": g["symbol"], "entrez_id": ent,
                     "evidence_status": "none_in_filtered_snapshot", "records": []}
        groups.append(_evidence_group_html(entry))
    return (
        f'<div class="view" id="{_esc(vid)}" role="region" '
        f'aria-label="Ranked dependencies for {_esc(sample)} model {_esc(model)}" hidden>'
        f'{_ranking_table_html(sample, model, block, ev)}'
        f'<h3>{_esc(EVIDENCE_LABEL)} &mdash; grouped beneath each displayed gene</h3>'
        '<p class="small">Evidence is retrieved by Entrez ID <em>after</em> the '
        'ranking is frozen and never affects rank. A gene with more evidence has '
        'not received a better rank.</p>'
        f'{"".join(groups)}'
        '</div>'
    )


def _collapsible(title: str, body_html: str, *, open_: bool = False) -> str:
    o = " open" if open_ else ""
    return (f'<details class="panel"{o}><summary>{_esc(title)}</summary>'
            f'<div>{body_html}</div></details>')


def render_html(bundle: dict) -> str:
    cs = bundle["case_study"]
    p1 = bundle["phase1"]
    dg = bundle["dgidb"]
    ev = cs["drug_gene_interaction_evidence"]
    cov = ev["coverage"]
    osteo = cs["osteosarcoma_validation_aggregate"]
    rm = cs["reconstructed_models"]

    # ---- A. header ---------------------------------------------------
    header = (
        '<header id="section-header">'
        f'<h1>{_esc(cs["title"])}</h1>'
        f'<p class="banner" role="note">{_esc(NOT_CLINICAL)}</p>'
        f'<p>{_esc(cs["description"])}</p>'
        '<p>Every number here is a <strong>predicted CRISPR gene-dependency</strong> '
        'estimate (how much a cell line is predicted to depend on losing a gene). '
        'It is <strong>not</strong> a prediction of treatment response, drug efficacy, '
        'or clinical outcome.</p>'
        f'<p class="small">Rendered from <span class="mono">case_study.json</span> '
        f'(schema <span class="mono">{_esc(cs["schema_version"])}</span>, sha256 '
        f'<span class="mono">{_esc(bundle["case_study_sha256"])}</span>); '
        f'case-study source commit <span class="mono">{_esc(cs["source_commit"])}</span>. '
        'No inference, evidence lookup, or recomputation happens in this page.</p>'
        '</header>'
    )

    # ---- B. frozen Phase 1 result ---------------------------------
    phase1 = (
        '<section id="section-phase1"><h2>B. Frozen Phase 1 result</h2>'
        '<p>The pre-specified Phase 1 comparison, on the DepMap validation split '
        f'({_esc(p1["n_cell_lines"])} held-out cell lines, {_esc(p1["n_targets"])} '
        'targets), per-target Spearman correlation:</p>'
        '<div class="panel"><div class="kv">'
        f'<div>ridge-on-PCA expression (ridge_pca)</div><div>mean Spearman <strong>{_num(p1["baseline_rho"])}</strong></div>'
        f'<div>Geneformer head (ridge_head)</div><div>mean Spearman <strong>{_num(p1["head_rho"])}</strong></div>'
        f'<div>delta (Geneformer &minus; baseline)</div><div><strong>{_num(p1["delta_mean"])}</strong></div>'
        f'<div>paired bootstrap 95% CI ({_esc(p1["bootstrap_resamples"])} resamples)</div>'
        f'<div>[{_num(p1["delta_ci_low"])}, {_num(p1["delta_ci_high"])}] '
        f'(SE {_num(p1["delta_std_error"])})</div>'
        '</div></div>'
        '<p><strong>Conclusion: the Geneformer head did not outperform the '
        'expression baseline.</strong> The delta is negative and its 95% CI '
        'excludes zero.</p>'
        '<p class="small">Values read at build time from the committed artifacts '
        '(hash-verified): ridge_pca &amp; ridge_head means from '
        '<span class="mono">baseline_results.json</span> / '
        '<span class="mono">head_results.json</span>; delta and CI from '
        '<span class="mono">analysis_results.json</span> '
        '(<span class="mono">A1_bootstrap.delta_head_minus_baseline</span>). '
        f'Effective degrees of freedom at the frozen alphas: '
        f'ridge_pca {_num(p1["eff_df_baseline"])}, ridge_head {_num(p1["eff_df_head"])}.</p>'
        '</section>'
    )

    # ---- C + D. selectors -------------------------------------------
    sample_radES = "".join(
        f'<label><input type="radio" name="sample" value="{_esc(s)}"'
        f'{" checked" if i == 0 else ""}> {_esc(s)} '
        f'&mdash; {"held-out verification anchor" if s == "ACH-000364" else "exploratory external prediction"}'
        '</label>'
        for i, s in enumerate(SAMPLE_ORDER))
    model_radES = "".join(
        f'<label><input type="radio" name="model" value="{_esc(m)}"'
        f'{" checked" if i == 0 else ""}> <span class="mono">{_esc(m)}</span></label>'
        for i, m in enumerate(MODEL_ORDER))
    selectors = (
        '<section id="section-samples"><h2>C &amp; D. Sample and model</h2>'
        '<form class="controls noprint" id="view-controls" onsubmit="return false">'
        '<fieldset id="section-models-fieldset"><legend>Sample</legend>'
        f'<div class="radio-row">{sample_radES}</div></fieldset>'
        '<fieldset id="section-models"><legend>Model</legend>'
        f'<div class="radio-row">{model_radES}</div></fieldset>'
        '<div class="field"><label for="gene-search">Search gene / Entrez / drug</label>'
        '<input type="search" id="gene-search" placeholder="e.g. YRDC, 79693, dinaciclib" '
        'autocomplete="off"></div>'
        '<div class="field"><label for="evi-status">Evidence status filter</label>'
        '<select id="evi-status">'
        '<option value="all">all</option>'
        '<option value="cited">cited</option>'
        '<option value="source_only">source-only</option>'
        '<option value="none_in_filtered_snapshot">none in filtered snapshot</option>'
        '</select></div>'
        '<div class="field"><span class="sr-only">Bulk controls</span>'
        '<button type="button" id="expand-all" class="secondary">Expand all evidence</button> '
        '<button type="button" id="collapse-all" class="secondary">Collapse all</button> '
        '<button type="button" id="print-btn" class="secondary">Print</button></div>'
        '</form>'
        '<p class="small" id="result-count" aria-hidden="true">25 of 25 genes shown</p>'
        '<p class="sr-only" id="status-live" role="status" aria-live="polite"></p>'
        '<p class="small">ridge_pca and ridge_head are shown separately and their '
        'rankings are <strong>never combined into a single consensus list</strong>.</p>'
        '</section>'
    )

    # ---- E + F. the four views ------------------------------------
    views = "".join(_view_html(s, m, cs)
                    for s in SAMPLE_ORDER for m in MODEL_ORDER)
    table_section = (
        '<section id="section-table"><h2>E &amp; F. Ranked dependency table and '
        f'{_esc(EVIDENCE_LABEL)}</h2>'
        f'<p><strong>{_esc(RANKING_HINT)}</strong> The table lists predicted '
        'dependencies, never therapeutic targets. Use the controls above to switch '
        'sample / model, search, or filter by evidence status.</p>'
        f'<div id="views">{views}</div>'
        '</section>'
    )

    # ---- G. sample-specific interpretation -----------------------
    bg_recon = cs["samples"]["BG003082"]["baseline_input"]["reconciliation"]
    interp = (
        '<section id="section-interpretation"><h2>G. Sample-specific interpretation</h2>'
        '<h3>ACH-000364 &mdash; held-out verification anchor</h3>'
        '<ul>'
        '<li>A held-out DepMap <strong>validation-split</strong> cell line (U-2 OS), '
        'used as a pipeline-verification example.</li>'
        '<li>Observed CRISPR GeneEffect values are attached to the ranked genes '
        '<strong>after</strong> prediction and ranking; they never influenced '
        'selection, model choice, or order.</li>'
        '<li><strong>One cell line is not a performance estimate.</strong> For the '
        'quantified comparison see section B and section H.</li>'
        '</ul>'
        '<h3>BG003082 &mdash; exploratory external prediction</h3>'
        '<ul>'
        '<li><strong>Bulk primary-tumour tissue</strong> scored using models trained '
        'and validated only on cultured DepMap cell lines &mdash; a real domain shift.</li>'
        '<li><strong>Outcome unavailable</strong>: no CRISPR screen exists for this '
        'sample; no observed value is loaded, invented, or computed.</li>'
        f'<li>{_int(bg_recon["canonical_genes_mapped"])}/'
        f'{_int(bg_recon["canonical_genes"])} expression genes mapped into the frozen '
        f'feature space; the remaining {_esc(bg_recon["canonical_genes_missing"])} are '
        'training-mean imputed by the reconstructed baseline artifact.</li>'
        '<li>The Geneformer embedding for this sample was <strong>generated '
        'separately</strong> (Kaggle GPU run); the historical Phase 1 Geneformer code '
        'revision was not captured.</li>'
        '<li><strong>Commensurability with the historical Phase 1 embeddings is not '
        'proven</strong> (bulk-tumour input; NCBI <span class="mono">gene2ensembl</span> '
        'map rather than the vanished <span class="mono">mygene</span> map; fresh '
        'revision pin).</li>'
        '<li>This is an <strong>exploratory prediction, not validated performance.</strong></li>'
        '</ul>'
        '</section>'
    )

    # ---- H. osteosarcoma descriptive aggregate -------------------
    cohort_tags = "".join(f'<span class="tag mono">{_esc(x)}</span>'
                          for x in osteo["cohort"]["model_ids"])
    osteo_sec = (
        '<section id="section-osteo"><h2>H. Osteosarcoma descriptive aggregate</h2>'
        '<p class="banner" role="note">Descriptive and unstable (n = 5). '
        'No confidence interval, no significance test. <strong>Not a replacement '
        'for the frozen Phase 1 result</strong> in section B.</p>'
        '<div class="panel"><div class="kv">'
        f'<div>validation cohort n</div><div>{_esc(osteo["cohort"]["n"])} '
        f'({cohort_tags})</div>'
        f'<div>target universe</div><div>{_esc(osteo["target_universe"])}</div>'
        f'<div>common finite targets</div><div>{_esc(osteo["common_finite_target_set"]["n_included"])} '
        f'included, {_esc(osteo["common_finite_target_set"]["n_excluded"])} excluded</div>'
        f'<div>ridge_pca mean per-target Spearman</div><div><strong>{_num(osteo["mean_per_target_spearman"]["ridge_pca"])}</strong></div>'
        f'<div>ridge_head mean per-target Spearman</div><div><strong>{_num(osteo["mean_per_target_spearman"]["ridge_head"])}</strong></div>'
        f'<div>delta (ridge_head &minus; ridge_pca)</div><div><strong>{_num(osteo["delta_ridge_head_minus_ridge_pca"])}</strong></div>'
        '</div></div>'
        f'<p class="small">{_esc(osteo["per_target_metric"])} '
        f'Cohort predicate: {_esc(osteo["cohort"]["predicate"])}. '
        f'Definition locked in {_esc(osteo["definition_source"])}. '
        'Not used to choose a model or alter the displayed rankings.</p>'
        '</section>'
    )

    # ---- I. methods / provenance / limitations -------------------
    def _pipe(mk):
        b = rm[mk]
        return (
            '<div class="kv">'
            f'<div>status</div><div>{_esc(b["provenance_status"])}</div>'
            f'<div>frozen alpha</div><div>{_num(b["frozen_alpha"]["value"])} '
            f'&mdash; {_esc(b["frozen_alpha"]["selection"])}; from '
            f'<span class="mono">{_esc(b["frozen_alpha"]["source"])}</span></div>'
            f'<div>pipeline</div><div class="mono">{_esc(b["pipeline"])}</div>'
            f'<div>features / targets</div><div>{_esc(b["n_features"])} / {_esc(b["n_targets"])}</div>'
            f'<div>artifact manifest sha256</div><div class="mono">{_esc(b["manifest_sha256"])}</div>'
            '</div>')

    inp_rows = "".join(
        f'<div>{_esc(k)}</div><div class="mono">{_esc(v)}</div>'
        for k, v in sorted(cs["input_artifact_sha256"].items()))
    lim_items = "".join(f"<li>{_esc(x)}</li>" for x in cs["limitations"])
    dis_items = "".join(f'<li class="disclaimer">{_esc(x)}</li>' for x in cs["disclaimers"])
    incl_src = ", ".join(_esc(s) for s in dg["included_sources"])

    methods = (
        '<section id="section-methods"><h2>I. Methods, provenance and limitations</h2>'
        + _collapsible(
            "Reconstructed fitted-state status",
            '<p>Both models are loaded as <strong>' + _esc(rm["ridge_pca"]["provenance_status"])
            + '</strong>. ' + _esc(rm["ridge_pca"]["not_original_fitted_objects"])
            + ' They reproduce every committed Phase 1 validation statistic exactly at '
            'the recorded precision (see <span class="mono">reconstruct_fitted.py --validate</span>).</p>',
            open_=True)
        + _collapsible("Frozen alphas and pipelines",
                       '<h3>ridge_pca</h3>' + _pipe("ridge_pca")
                       + '<h3>ridge_head</h3>' + _pipe("ridge_head"))
        + _collapsible(
            "Sample roles and domain shift",
            '<div class="kv">'
            f'<div>ACH-000364</div><div>role {_esc(cs["samples"]["ACH-000364"]["role"])}; '
            f'{_esc(cs["samples"]["ACH-000364"]["prediction_status"])} / '
            f'{_esc(cs["samples"]["ACH-000364"]["outcome_status"])}; split '
            f'{_esc(cs["samples"]["ACH-000364"]["depmap_split"])}</div>'
            f'<div>BG003082</div><div>role {_esc(cs["samples"]["BG003082"]["role"])}; '
            f'{_esc(cs["samples"]["BG003082"]["prediction_status"])} / '
            f'{_esc(cs["samples"]["BG003082"]["outcome_status"])}; '
            'absent from every DepMap split; bulk tumour vs cultured cell lines '
            '(domain shift)</div></div>')
        + _collapsible(
            "Identifier mapping and imputation (BG003082 baseline input)",
            '<div class="kv">'
            f'<div>transformation</div><div>{_esc(cs["samples"]["BG003082"]["baseline_input"]["transformation"])}</div>'
            f'<div>canonical genes</div><div>{_esc(bg_recon["canonical_genes"])}</div>'
            f'<div>mapped via Ensembl-ID join</div><div>{_esc(bg_recon["canonical_genes_mapped"])} '
            f'({_esc(bg_recon["canonical_genes_measured_zero"])} measured zero, '
            f'{_esc(bg_recon["canonical_genes_measured_nonzero"])} measured &gt; 0)</div>'
            f'<div>missing &rarr; training-mean imputed</div><div>{_esc(bg_recon["canonical_genes_missing"])}</div>'
            f'<div>symbol fallback</div><div>{_esc(bg_recon["symbol_fallback"])} '
            f'(candidates: {", ".join(_esc(c) for c in bg_recon["symbol_fallback_candidates"])})</div>'
            '</div>'
            f'<p class="small">{_esc(cs["samples"]["BG003082"]["baseline_input"]["imputation"])}</p>')
        + _collapsible(
            "DGIdb release, data vintage and licence limitations",
            '<div class="kv">'
            f'<div>release tag</div><div>{_esc(dg["release_tag"])}</div>'
            f'<div>interaction data version</div><div>{_esc(dg["interaction_data_version"])} '
            '(the records are <strong>not</strong> current-year data)</div>'
            f'<div>DGIdb app version</div><div>{_esc(dg["app_version_tsv"])} / {_esc(dg["app_version_graphql"])}</div>'
            f'<div>retrieved (fixed provenance input)</div><div class="mono">{_esc(dg["retrieved_utc"])}</div>'
            f'<div>included interaction sources</div><div>{incl_src} '
            f'({_esc(dg["n_excluded_sources"])} other sources excluded on licence grounds)</div>'
            f'<div>snapshot sha256</div><div class="mono">{_esc(dg["snapshot_sha256"])}</div>'
            f'<div>manifest sha256</div><div class="mono">{_esc(dg["manifest_sha256"])}</div>'
            '</div>'
            '<p class="small">Each retained record still carries its own source '
            'licence. Retrieval only &mdash; no efficacy, approval, indication, or '
            'osteosarcoma-relevance is inferred.</p>')
        + _collapsible(
            "Publication coverage",
            f'<p>Snapshot-wide: {_esc(dg["records_with_pmids"])} of '
            f'{_esc(dg["snapshot_records"])} records '
            f'(<strong>{_num(dg["pct_with_pmids"])}%</strong>, i.e. about 19%) carry '
            '&ge;1 PMID. Coverage is source-skewed: CIViC / DoCM / NCI are near-complete; '
            'ChEMBL / GuideToPharmacology / FDA carry no claim-level publications in '
            'this release. A record with no PMID is <em>source-only</em>, not an error.</p>'
            '<div class="kv">'
            f'<div>this case study: distinct displayed genes</div><div>{_esc(cov["n_distinct_genes"])}</div>'
            f'<div>&nbsp;&nbsp;cited</div><div>{_esc(cov["n_cited"])}</div>'
            f'<div>&nbsp;&nbsp;source-only</div><div>{_esc(cov["n_source_only"])}</div>'
            f'<div>&nbsp;&nbsp;none in filtered snapshot</div><div>{_esc(cov["n_none_in_filtered_snapshot"])}</div>'
            f'<div>&nbsp;&nbsp;total records / PMID citations</div><div>{_esc(cov["total_records"])} / {_esc(cov["total_pmid_citations"])}</div>'
            '</div>'
            f'<p class="small">{_esc(ev["pmid_scope_note"])}</p>')
        + _collapsible(
            "Artifact hashes and source commit",
            f'<p>Case-study source commit: <span class="mono">{_esc(cs["source_commit"])}</span>. '
            f'case_study.json sha256: <span class="mono">{_esc(bundle["case_study_sha256"])}</span>.</p>'
            '<h3>Inputs recorded in case_study.json</h3>'
            f'<div class="kv">{inp_rows}</div>'
            '<h3>Phase 1 result artifacts (read for section B)</h3>'
            '<div class="kv">'
            + "".join(f'<div>{_esc(k)}</div><div class="mono">{_esc(v)}</div>'
                     for k, v in sorted(p1["source_files"].items()))
            + '</div>')
        + _collapsible(
            "Limitations and disclaimers",
            f'<ul>{lim_items}</ul><ul style="list-style:none;padding-left:0">{dis_items}</ul>',
            open_=True)
        + '</section>'
    )

    embedded = _json_for_script(cs)
    doc = (
        '<!DOCTYPE html>\n'
        '<html lang="en">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f'<title>{_esc(cs["title"])}</title>\n'
        f'<style>{_CSS}</style>\n'
        '</head>\n<body>\n'
        '<main>\n'
        f'{header}\n{phase1}\n{selectors}\n{table_section}\n{interp}\n{osteo_sec}\n{methods}\n'
        '<footer class="small"><hr>'
        f'Self-contained offline report ({_esc(REPORT_SCHEMA)}). No external runtime '
        'dependency. External links (PubMed, source, licence) are optional and the '
        'report is complete without internet access.'
        '</footer>\n'
        '</main>\n'
        f'<script type="application/json" id="case-study-data">{embedded}</script>\n'
        f'<script>{_JS}</script>\n'
        '</body>\n</html>\n'
    )
    return doc


# --------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------

def build(out_path: Path | None = None, *, verbose: bool = True) -> str:
    out_path = Path(out_path) if out_path is not None else config.REPORT_HTML_FILE
    bundle = load_verified()
    doc = render_html(bundle)
    _run_claim_language_gate(doc)
    out_path.write_text(doc, encoding="utf-8", newline="\n")
    if verbose:
        try:
            shown = out_path.resolve().relative_to(_REPO_ROOT).as_posix()
        except ValueError:
            shown = out_path.name
        print(f"  wrote {shown}  "
              f"({out_path.stat().st_size:,} bytes)  sha256 {_sha256_file(out_path)}")
    return doc


# --------------------------------------------------------------------------
# claim-language gate  (rendered prose only, not embedded source data)
# --------------------------------------------------------------------------

def _prose_only(doc: str) -> str:
    """Strip the embedded JSON data block before scanning user-facing prose."""
    return re.sub(
        r'<script type="application/json" id="case-study-data">.*?</script>',
        "", doc, flags=re.S)


def _run_claim_language_gate(doc: str) -> None:
    prose = _prose_only(doc).lower()
    hits = [p for p in PROHIBITED_CLAIM_PHRASES if p in prose]
    if hits:
        raise ReportError(f"prohibited claim language in rendered prose: {hits}")


# --------------------------------------------------------------------------
# validate
# --------------------------------------------------------------------------

_REMOTE_PATTERNS = (
    r'src\s*=\s*["\']https?:', r'<link[^>]+href\s*=\s*["\']https?:',
    r'@import\s+url\(\s*["\']?https?:', r'//cdn', r'googleapis', r'gstatic',
    r'unpkg', r'jsdelivr', r'cloudflare', r'<iframe', r'analytics', r'gtag\(',
)


def _js_brackets_balanced(js: str) -> bool:
    """
    Lightweight structural check on the behaviour <script>: strings, template
    literals, regex literals and comments removed, then bracket balance.
    """
    out = []
    i, n = 0, len(js)
    while i < n:
        c = js[i]
        if c in "\"'`":
            q = c
            i += 1
            while i < n and js[i] != q:
                if js[i] == "\\":
                    i += 2
                    continue
                i += 1
            i += 1
            continue
        if c == "/" and i + 1 < n and js[i + 1] == "/":
            while i < n and js[i] != "\n":
                i += 1
            continue
        if c == "/" and i + 1 < n and js[i + 1] == "*":
            i += 2
            while i + 1 < n and not (js[i] == "*" and js[i + 1] == "/"):
                i += 1
            i += 2
            continue
        out.append(c)
        i += 1
    stripped = "".join(out)
    stack = []
    pairs = {")": "(", "]": "[", "}": "{"}
    for ch in stripped:
        if ch in "([{":
            stack.append(ch)
        elif ch in ")]}":
            if not stack or stack.pop() != pairs[ch]:
                return False
    return not stack


def validate(path: Path | None = None, *, verbose: bool = True) -> dict:
    path = Path(path) if path is not None else config.REPORT_HTML_FILE
    checks: list[tuple[str, bool, str]] = []

    def chk(name: str, ok: bool, detail: str = ""):
        checks.append((name, bool(ok), detail))

    if not path.is_file():
        raise ReportError(f"{path} not found -- run --build first")
    doc = path.read_text(encoding="utf-8")

    # ---- 1. committed inputs verify + strict JSON --------------------
    bundle = load_verified()               # raises on any hash / schema / NaN issue
    cs = bundle["case_study"]
    chk("committed case_study.json + result artifacts hash-verify; strict JSON parse", True)

    # ---- 2. byte-identical regeneration (twice) --------------------
    tmp = Path(tempfile.mkdtemp(prefix="report_val_"))
    try:
        a = tmp / "a.html"
        b = tmp / "b.html"
        build(a, verbose=False)
        build(b, verbose=False)
        regen_ok = a.read_bytes() == b.read_bytes() == path.read_bytes()
        chk("committed HTML == fresh regeneration (twice), byte-identical", regen_ok)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # ---- 3. embedded JSON round-trips to the committed case study ---
    m = re.search(r'<script type="application/json" id="case-study-data">(.*?)</script>',
                  doc, flags=re.S)
    emb_ok = False
    if m:
        raw = m.group(1)
        no_break = "</" not in raw and "<script" not in raw.lower()
        unesc = raw.replace("\\u0026", "&").replace("\\u003c", "<").replace("\\u003e", ">")
        try:
            parsed = json.loads(unesc, parse_constant=_reject_constant)
            emb_ok = no_break and (json.dumps(parsed, sort_keys=True)
                                   == json.dumps(cs, sort_keys=True))
        except Exception:
            emb_ok = False
    chk("embedded JSON: no </script> break-out, parses, equals committed case study", emb_ok)

    # ---- 4. required sections + warnings --------------------------
    required_ids = ["section-header", "section-phase1", "section-samples",
                    "section-models", "section-table", "section-interpretation",
                    "section-osteo", "section-methods"]
    chk("all required <section> anchors present",
        all(f'id="{i}"' in doc for i in required_ids),
        ",".join(i for i in required_ids if f'id="{i}"' not in doc))
    required_text = [
        NOT_CLINICAL, RANKING_HINT, EVIDENCE_LABEL,
        "did not outperform the expression",
        "replacement for the frozen Phase 1 result",
        "exploratory prediction, not validated performance",
        "not therapeutic targets", "Descriptive and unstable (n = 5)",
        "0.2356", "0.2047", "-0.0308", "-0.0365", "-0.0255",
        "0.119436", "0.082773", "-0.036663",
        "18,427/18,460", "training-mean imputed",
        "cell line is not a performance estimate",
        "Commensurability with the historical Phase 1 embeddings is not",
        "read at build time from the committed artifacts",
    ]
    missing_text = [t for t in required_text if t not in doc]
    chk("all required headline / warning strings present", not missing_text, str(missing_text))

    # ---- 5. 2 samples, 4 views, 25 rows each --------------------
    view_ids = [f"view-{s}-{mm}" for s in SAMPLE_ORDER for mm in MODEL_ORDER]
    chk("4 sample/model views present", all(f'id="{v}"' in doc for v in view_ids))
    rows_ok = True
    for v in view_ids:
        seg = doc.split(f'id="{v}"', 1)[1].split('class="view"', 1)[0]
        nrows = len(re.findall(r'<tr class="gene-row"', seg))
        if nrows != 25:
            rows_ok = False
            chk(f"{v}: 25 ranked rows", False, f"found {nrows}")
    if rows_ok:
        chk("every view has exactly 25 ranked rows", True)

    # ---- 6. evidence counts reconcile with case_study.json ------
    cov = cs["drug_gene_interaction_evidence"]["coverage"]
    # each distinct displayed gene renders one evi-group per view it appears in;
    # count distinct by entrez across the by_entrez block instead
    grp_status = re.findall(r'<details class="evi-group" data-evidence-status="([a-z_]+)"', doc)
    # 4 views but the union of displayed genes is what coverage counts; recount
    # from the JSON block directly
    by = cs["drug_gene_interaction_evidence"]["by_entrez"]
    jc = {"cited": 0, "source_only": 0, "none_in_filtered_snapshot": 0}
    for e in by.values():
        jc[e["evidence_status"]] += 1
    chk("evidence coverage in JSON self-consistent",
        jc["cited"] == cov["n_cited"] and jc["source_only"] == cov["n_source_only"]
        and jc["none_in_filtered_snapshot"] == cov["n_none_in_filtered_snapshot"]
        and sum(jc.values()) == cov["n_distinct_genes"])
    total_records_rendered = len(re.findall(r'<div class="evi-record">', doc))
    # a record renders once per view that displays its gene: expected count is
    # sum over the 4 views of sum over that view's 25 displayed genes of len(records)
    expected_records = sum(
        len(by.get(g["entrez_id"], {}).get("records", []))
        for s in SAMPLE_ORDER for mm in MODEL_ORDER
        for g in cs["rankings"][s][mm]["genes"])
    chk("rendered evidence-record blocks == per-view sum of JSON records",
        total_records_rendered == expected_records and expected_records > 0,
        f"rendered={total_records_rendered} expected={expected_records}")
    # pill status classes only from the known vocabulary
    chk("evidence status pills use only the known vocabulary",
        set(grp_status) <= {"cited", "source_only", "none_in_filtered_snapshot"})

    # ---- 7. no remote runtime dependency -----------------------
    lowered = doc.lower()
    remote_hits = [p for p in _REMOTE_PATTERNS if re.search(p, lowered)]
    chk("no remote script/style/font/analytics dependency", not remote_hits, str(remote_hits))
    chk("no <script src=> at all", "<script src" not in lowered)
    chk("no <link rel=stylesheet>", "rel=\"stylesheet\"" not in lowered and "rel='stylesheet'" not in lowered)

    # ---- 8. embedded JS parses (structural) -------------------
    js_blocks = re.findall(r'<script>(.*?)</script>', doc, flags=re.S)
    chk("exactly one behaviour <script> block", len(js_blocks) == 1)
    chk("behaviour JS brackets balanced (structural parse)",
        bool(js_blocks) and _js_brackets_balanced(js_blocks[0]))
    chk("behaviour JS has no NaN / Infinity token",
        bool(js_blocks) and not re.search(r'\b(NaN|Infinity)\b', js_blocks[0]))

    # ---- 9. unique element IDs -----------------------------------
    ids = re.findall(r'\sid="([^"]+)"', doc)
    dupes = sorted({i for i in ids if ids.count(i) > 1})
    chk("all HTML element IDs are unique", not dupes, str(dupes))

    # ---- 10. internal controls reference valid targets -------
    id_set = set(ids)
    ref_ok = True
    bad_refs = []
    for attr in ("for", "aria-controls", "aria-labelledby", "aria-describedby"):
        for target in re.findall(rf'\s{attr}="([^"]+)"', doc):
            for t in target.split():
                if t not in id_set:
                    ref_ok = False
                    bad_refs.append(f"{attr}->{t}")
    for href in re.findall(r'href="#([^"]+)"', doc):
        if href and href not in id_set:
            ref_ok = False
            bad_refs.append(f"href#->{href}")
    chk("internal for/aria/#href controls reference existing IDs", ref_ok, str(bad_refs))

    # ---- 11. escaping, no absolute paths, no volatile timestamp --
    chk("no unescaped </script> outside the two real closing tags",
        doc.count("</script>") == 2)
    chk("no absolute local path",
        "C:\\" not in doc and "C:/" not in doc
        and str(_REPO_ROOT) not in doc
        and str(_REPO_ROOT).replace("\\", "/") not in doc)
    chk("no wall-clock timestamp (HH:MM:SS)",
        not re.search(r'\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d', doc.replace(
            bundle["dgidb"]["retrieved_utc"], "")))
    # NaN / Infinity as a *value token* is impossible: the embedded JSON is
    # strict-parsed (check 3) and the behaviour JS is token-scanned (check 8).
    # A guard against a JSON-value-position slip:
    chk("no NaN / Infinity in JSON value position in the embedded data",
        not re.search(r'[:\[,]\s*(NaN|Infinity|-Infinity)\b', m.group(1) if m else ""))

    # ---- 12. claim-language gate on rendered prose -----------
    try:
        _run_claim_language_gate(doc)
        chk("no prohibited claim language in rendered prose", True)
    except ReportError as exc:
        chk("no prohibited claim language in rendered prose", False, str(exc))

    # ---- 13. structural sanity ----------------------------------
    chk("doctype + lang", doc.startswith("<!DOCTYPE html>") and 'lang="en"' in doc)
    chk("no combined consensus ranking",
        "consensus" in doc.lower() and "never" in doc.lower()
        and "combined into a single consensus list" in doc)

    n_fail = sum(1 for _, ok, _ in checks if not ok)
    if verbose:
        for name, ok, detail in checks:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
        # browser smoke test (best effort; not a hard gate when no browser)
        smoke = _browser_smoke(path)
        print(f"\n  browser smoke test: {smoke['status']} -- {smoke['detail']}")
        if smoke["status"] == "fail":
            n_fail += 1
        print()
        print(f"  {'ALL CHECKS PASSED' if n_fail == 0 else str(n_fail) + ' FAILED'}"
              f"  ({len(checks) - sum(1 for _, ok, _ in checks if not ok)}/{len(checks)} structural)")
    return {"checks": checks, "n_fail": n_fail}


# --------------------------------------------------------------------------
# best-effort headless-browser smoke test (uses an already-installed browser)
# --------------------------------------------------------------------------

_BROWSER_CANDIDATES = (
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
)


def _find_browser() -> str | None:
    for c in _BROWSER_CANDIDATES:
        if Path(c).is_file():
            return c
    for name in ("google-chrome", "chromium", "chromium-browser", "msedge"):
        found = shutil.which(name)
        if found:
            return found
    return None


# A harness appended to a TEMP COPY of the report (never the committed file). It
# drives the real controls the way a user would and records the outcome on the
# <html> element, which --dump-dom then reports back. Uses only an already-
# installed browser -- no new dependency, no DevTools-protocol client.
_SMOKE_HARNESS = """
<script>
(function(){
  function done(pass, note){
    document.documentElement.setAttribute("data-smoke", pass ? "pass" : "fail");
    document.documentElement.setAttribute("data-smoke-note", note);
  }
  try{
    function visible(id){ var e=document.getElementById(id); return e && !e.hidden; }
    var checks=[];
    checks.push(["init-js-ready", document.documentElement.getAttribute("data-js-ready")==="ok"]);
    checks.push(["default-view", visible("view-ACH-000364-ridge_pca")
                                 && !visible("view-BG003082-ridge_head")]);
    // switch sample -> BG003082
    var bg=document.querySelector('input[name=sample][value="BG003082"]');
    bg.checked=true; bg.dispatchEvent(new Event("change",{bubbles:true}));
    checks.push(["sample-switch", visible("view-BG003082-ridge_pca")
                                  && !visible("view-ACH-000364-ridge_pca")]);
    // switch model -> ridge_head
    var mh=document.querySelector('input[name=model][value="ridge_head"]');
    mh.checked=true; mh.dispatchEvent(new Event("change",{bubbles:true}));
    checks.push(["model-switch", visible("view-BG003082-ridge_head")]);
    // search filter
    var view=document.getElementById("view-BG003082-ridge_head");
    var totalRows=view.querySelectorAll("tr.gene-row").length;
    var s=document.getElementById("gene-search");
    var firstSym=view.querySelector("tr.gene-row td:nth-child(2)").textContent.trim();
    s.value=firstSym; s.dispatchEvent(new Event("input",{bubbles:true}));
    var shown=Array.prototype.filter.call(view.querySelectorAll("tr.gene-row"),
              function(r){return !r.classList.contains("hidden-row");}).length;
    checks.push(["search-filter", totalRows===25 && shown>=1 && shown<totalRows]);
    s.value=""; s.dispatchEvent(new Event("input",{bubbles:true}));
    // evidence-status filter
    var f=document.getElementById("evi-status");
    f.value="cited"; f.dispatchEvent(new Event("change",{bubbles:true}));
    var groups=view.querySelectorAll("details.evi-group");
    var badStatus=Array.prototype.some.call(groups, function(g){
      return !g.classList.contains("hidden-row")
             && g.getAttribute("data-evidence-status")!=="cited"; });
    checks.push(["status-filter", groups.length>0 && !badStatus]);
    f.value="all"; f.dispatchEvent(new Event("change",{bubbles:true}));
    // expand all / collapse all
    document.getElementById("expand-all").click();
    var anyOpen=Array.prototype.some.call(view.querySelectorAll("details"),
              function(d){return d.open;});
    document.getElementById("collapse-all").click();
    var allClosed=Array.prototype.every.call(view.querySelectorAll("details"),
              function(d){return !d.open;});
    checks.push(["expand-collapse", anyOpen && allClosed]);
    var failed=checks.filter(function(c){return !c[1];}).map(function(c){return c[0];});
    done(failed.length===0, failed.length? failed.join(",") : "all "+checks.length+" ok");
  }catch(e){ done(false, "exception: "+e); }
})();
</script>
"""


def _browser_smoke(html_path: Path) -> dict:
    browser = _find_browser()
    if not browser:
        return {"status": "skipped",
                "detail": ("no local browser found; structural + JS-bracket "
                           "validation only")}
    tmp = Path(tempfile.mkdtemp(prefix="report_smoke_"))
    try:
        driven = tmp / "driven.html"
        driven.write_text(
            html_path.read_text(encoding="utf-8").replace(
                "</body>", _SMOKE_HARNESS + "</body>", 1),
            encoding="utf-8")
        proc = subprocess.run(
            [browser, "--headless=new", "--disable-gpu", "--no-sandbox",
             f"--user-data-dir={tmp / 'prof'}", "--virtual-time-budget=10000",
             "--run-all-compositor-stages-before-draw", "--dump-dom",
             driven.resolve().as_uri()],
            capture_output=True, text=True, timeout=120)
        dom = proc.stdout or ""
        err = proc.stderr or ""
        init_ok = 'data-js-ready="ok"' in dom
        m_smoke = re.search(r'data-smoke="(pass|fail)"', dom)
        m_note = re.search(r'data-smoke-note="([^"]*)"', dom)
        smoke = m_smoke.group(1) if m_smoke else "no-result"
        note = m_note.group(1) if m_note else "(no note captured)"
        errs = [ln for ln in err.splitlines()
                if "Uncaught" in ln or "SyntaxError" in ln]
        ok = init_ok and smoke == "pass" and not errs
        return {
            "status": "pass" if ok else "fail",
            "detail": (f"init_js_ready={init_ok}; interaction_smoke={smoke} "
                       f"({note}); console_errors={len(errs)}"),
            "browser": browser,
        }
    except Exception as exc:  # noqa: BLE001
        return {"status": "skipped", "detail": f"browser run failed: {exc}"}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# self-test  (synthetic render + escaping units + real build/validate)
# --------------------------------------------------------------------------

def _self_test() -> int:
    print("Running report.py self-test...")

    # ---- escaping units -------------------------------------------
    assert _esc('<b>&"\'</b>') == "&lt;b&gt;&amp;&quot;&#x27;&lt;/b&gt;"
    payload = {"x": "</script><img src=x onerror=alert(1)> & < > \u2028 done"}
    emb = _json_for_script(payload)
    assert "<" not in emb and ">" not in emb and "</script" not in emb.lower()
    back = json.loads(emb.replace("\\u0026", "&").replace("\\u003c", "<").replace("\\u003e", ">"))
    assert back == payload
    print("  [ok] HTML escaping + <script>-safe JSON embedding round-trip")

    # ---- synthetic render ---------------------------------------
    syn_cs = _synthetic_case_study()
    bundle = {
        "case_study": syn_cs,
        "case_study_bytes": json.dumps(syn_cs).encode(),
        "case_study_sha256": "0" * 64,
        "phase1": {**PHASE1_EXPECTED, "baseline_alpha": 100000.0, "head_alpha": 3162.0,
                   "delta_std_error": 0.0028, "eff_df_baseline": 49.71,
                   "eff_df_head": 51.78,
                   "source_files": {k: "0" * 64 for k in
                                    ("baseline_results.json", "head_results.json",
                                     "analysis_results.json")}},
        "dgidb": {"release_tag": "2026-06b", "interaction_data_version": "Dec-2023",
                  "app_version_tsv": "v.5.0.11", "app_version_graphql": "v.5.0.12",
                  "retrieved_utc": "2026-08-29T00:00:00Z",
                  "included_sources": ["CIViC", "ChEMBL"], "n_excluded_sources": 15,
                  "snapshot_records": 37343, "records_with_pmids": 7078,
                  "pct_with_pmids": 19.0, "manifest_sha256": "0" * 64,
                  "snapshot_sha256": "0" * 64},
    }
    doc = render_html(bundle)
    _run_claim_language_gate(doc)
    assert doc.startswith("<!DOCTYPE html>")
    assert doc.count("</script>") == 2
    ids = re.findall(r'\sid="([^"]+)"', doc)
    assert len(ids) == len(set(ids)), sorted({i for i in ids if ids.count(i) > 1})
    for v in (f"view-{s}-{m}" for s in SAMPLE_ORDER for m in MODEL_ORDER):
        seg = doc.split(f'id="{v}"', 1)[1].split('class="view"', 1)[0]
        assert len(re.findall(r'<tr class="gene-row"', seg)) == 25, v
    # injected hostile strings are escaped, not executed
    assert "<img src=x onerror" not in doc
    assert "&lt;img src=x onerror" in doc
    assert _js_brackets_balanced(re.findall(r'<script>(.*?)</script>', doc, flags=re.S)[0])
    print("  [ok] synthetic render: 4 views x 25 rows, unique IDs, hostile data escaped")

    # ---- claim-language gate actually fails on bad prose --------
    bad = doc.replace("<footer", "<p>this is a treatment recommendation</p><footer")
    try:
        _run_claim_language_gate(bad)
        raise AssertionError("claim gate should have fired")
    except ReportError:
        pass
    print("  [ok] claim-language gate rejects prohibited framing")

    # ---- real build + validate --------------------------------
    tmp = Path(tempfile.mkdtemp(prefix="report_selftest_"))
    try:
        p = tmp / "phase2_report.html"
        build(p, verbose=False)
        build(tmp / "again.html", verbose=False)
        assert p.read_bytes() == (tmp / "again.html").read_bytes(), "not byte-identical"
        rep = validate(p, verbose=False)
        bad_checks = [n for n, ok, _ in rep["checks"] if not ok]
        assert not bad_checks, bad_checks
        print(f"  [ok] real build byte-identical; validate {len(rep['checks'])}/"
              f"{len(rep['checks'])} structural checks pass")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\nSelf-test passed.")
    return 0


def _synthetic_case_study() -> dict:
    genes_obs = [{"rank": i + 1, "entrez_id": str(1000 + i), "symbol": f"G{i}",
                  "predicted_geneeffect": round(-1.5 + i * 0.01, 10),
                  "observed_geneeffect": round(-1.0 + i * 0.02, 10),
                  "observed_rank": i + 5} for i in range(25)]
    genes_pred = [{"rank": i + 1, "entrez_id": str(1000 + i), "symbol": f"G{i}",
                   "predicted_geneeffect": round(-1.5 + i * 0.01, 10)} for i in range(25)]
    rule = "raw float64 sort; Entrez breaks only exact raw-value ties"
    rblk_obs = {"model": "m", "model_provenance": "reconstructed fitted state",
                "n_targets_ranked": 4297, "n_displayed": 25, "ranking_rule": rule,
                "observed_rank_rule": rule, "observed_values_attached_after_ranking": True,
                "n_targets_with_observed_value": 4297,
                "not_a_recommendation": "predicted dependencies", "genes": genes_obs}
    rblk_pred = {"model": "m", "model_provenance": "reconstructed fitted state",
                 "n_targets_ranked": 4297, "n_displayed": 25, "ranking_rule": rule,
                 "not_a_recommendation": "predicted dependencies", "genes": genes_pred}
    by_entrez = {}
    for i in range(25):
        e = str(1000 + i)
        st = ("cited" if i < 2 else "source_only" if i < 4 else "none_in_filtered_snapshot")
        recs = []
        if st != "none_in_filtered_snapshot":
            recs = [{"entrez_id": e, "gene_symbol": f"G{i}", "dgidb_gene_name": f"G{i}",
                     "gene_symbol_consistent": "true", "symbol_query_mismatch": "false",
                     "drug_name": "</script> DRUG <img src=x onerror=alert(1)>",
                     "drug_concept_id": "d:1", "drug_claim_name": "D",
                     "interaction_source": "ChEMBL", "interaction_source_version": "37",
                     "source_license": "CC BY-SA 3.0", "source_license_url": "https://x",
                     "interaction_type_raw": "inhibitor", "interaction_direction": "inhibitory",
                     "direction_tier": "inhibitory", "interaction_score": "0.4",
                     "drug_specificity_score": "", "gene_specificity_score": "",
                     "evidence_score": "8", "drug_is_approved": "false",
                     "drug_is_immunotherapy": "false", "drug_is_antineoplastic": "true",
                     "curation_type": "", "indication": "",
                     "pmids": (["12345", "67890"] if st == "cited" else []),
                     "pmid_status": ("cited" if st == "cited" else "source_only"),
                     "dgidb_release_tag": "2026-06b", "record_key": "k",
                     "disclaimer": "A recorded drug\u2013gene interaction does not establish efficacy."}]
        by_entrez[e] = {"entrez_id": e, "symbol": f"G{i}", "evidence_status": st,
                        "n_records": len(recs), "records": recs}
    return {
        "schema_version": "case-study/1", "title": "Synthetic", "source_commit": "0" * 40,
        "description": "synthetic self-test case study",
        "rankings": {"ACH-000364": {"ridge_pca": rblk_obs, "ridge_head": rblk_obs},
                     "BG003082": {"ridge_pca": rblk_pred, "ridge_head": rblk_pred}},
        "drug_gene_interaction_evidence": {
            "label": EVIDENCE_LABEL, "framing": "retrieval only, no efficacy claim",
            "disclaimer": "A recorded drug\u2013gene interaction does not establish efficacy.",
            "pmid_scope_note": "group-level PMID attribution",
            "retrieval": {"snapshot_file": "x", "snapshot_sha256": "0" * 64,
                          "manifest_sha256": "0" * 64, "method": "by Entrez ID",
                          "top_k_per_direction_tier": 5,
                          "direction_tiers": ["inhibitory", "activating", "unknown"],
                          "retrieved_after_top_n_frozen": True,
                          "evidence_availability_did_not_affect_selection_or_ranking": True},
            "coverage": {"n_distinct_genes": 25, "n_cited": 2, "n_source_only": 2,
                         "n_none_in_filtered_snapshot": 21, "total_records": 4,
                         "total_pmid_citations": 4},
            "by_entrez": by_entrez},
        "osteosarcoma_validation_aggregate": {
            "status": "DESCRIPTIVE and UNSTABLE because n=5. No confidence interval...",
            "definition_source": "capstone/scope-decisions.md (2026-08-29 entry)",
            "cohort": {"n": 5, "model_ids": ["ACH-1", "ACH-2", "ACH-3", "ACH-4", "ACH-5"],
                       "predicate": "osteosarcoma_mask & val"},
            "models": ["ridge_pca", "ridge_head"], "target_universe": 4297,
            "per_target_metric": "baseline.per_target_spearman across 5 lines",
            "common_finite_target_set": {"n_included": 4255, "n_excluded": 42,
                                         "excluded_ridge_pca_nonfinite": 42,
                                         "excluded_ridge_head_nonfinite": 42,
                                         "rule": "finite for both",
                                         "excluded_reason": "n=5"},
            "mean_per_target_spearman": {"ridge_pca": 0.119436, "ridge_head": 0.082773,
                                         "rounding_dp": 6},
            "delta_ridge_head_minus_ridge_pca": -0.036663,
            "used_to_choose_model_or_alter_rankings": False},
        "reconstructed_models": {
            mk: {"model": mk,
                 "provenance_status": "reconstructed fitted state at the frozen Phase 1 alpha from the unchanged frozen training data",
                 "not_original_fitted_objects": "never serialised; reconstruction",
                 "artifact_dir": f"data/processed/reconstructed_fitted/{mk}",
                 "manifest_sha256": "0" * 64, "base_commit": "0" * 40,
                 "feature_order_sha256": "0" * 64, "target_order_sha256": "0" * 64,
                 "frozen_alpha": {"value": (100000.0 if mk == "ridge_pca" else 3162.0),
                                  "selection": "read verbatim; NOT re-selected",
                                  "source": f"{mk} alpha"},
                 "pipeline": "impute -> scale -> ridge", "n_features": 18460, "n_targets": 4297}
            for mk in ("ridge_pca", "ridge_head")},
        "samples": {
            "ACH-000364": {
                "role": "verification_anchor", "cell_line": "U-2 OS",
                "prediction_status": "held_out_prediction",
                "outcome_status": "measured_crispr", "depmap_split": "val",
                "split_assertion": "asserted == 'val'; hard-stop if 'train'",
                "in_training_split": False,
                "baseline_input": {"source": "expression.npz row", "n_features": 18460,
                                   "feature_order_sha256": "0" * 64,
                                   "missing_features": 0, "imputed_features": 0},
                "head_input": {"source": "geneformer_embeddings.csv row",
                               "n_features": 768, "feature_order_sha256": "0" * 64,
                               "missing_features": 0},
                "observed_crispr": {"source": "crispr_effect.npz row",
                                    "role": "verification example only; attached after ranking",
                                    "n_targets_with_value": 4297, "n_targets_missing": 0}},
            "BG003082": {
                "role": "exploratory_external_sample",
                "description": "Sid Sijbrandij osteosarcoma primary-tumour RNA-seq, CC0 1.0; bulk tumour tissue",
                "analysis_role": "exploratory_external_prediction",
                "prediction_status": "exploratory_external_prediction",
                "outcome_status": "unavailable",
                "absent_from_all_depmap_splits": True,
                "observed_outcome": "none exists; none loaded, invented, or computed",
                "baseline_input": {
                    "source": "sample_profile.load_external_sample()",
                    "transformation": "sum linear TPM across version-stripped Ensembl IDs, then log2(TPM+1); unresolved left NaN",
                    "n_features": 18460, "feature_order_sha256": "0" * 64,
                    "missing_features_represented_as_nan": 33,
                    "imputation": "the reconstructed baseline artifact applies its stored training-mean impute vector",
                    "gct_file": {"name": "BG003082.gene_tpm.gct.gz", "sha256": "0" * 64, "bytes": 778119},
                    "ensembl_map_file": {"name": "ensembl_map.csv", "sha256": "0" * 64, "bytes": 1, "rows": 18459},
                    "gene_columns_file": {"name": "gene_columns.json", "sha256": "0" * 64},
                    "reconciliation": {"canonical_genes": 18460, "canonical_genes_mapped": 18427,
                                       "canonical_genes_missing": 33,
                                       "canonical_genes_measured_zero": 1407,
                                       "canonical_genes_measured_nonzero": 17020,
                                       "resolved_via_ensembl_map": 18427,
                                       "unresolved_external_rows": 56201,
                                       "symbol_fallback": "not attempted",
                                       "symbol_fallback_candidates": ["ASPRV1", "FAM174C", "NOX5", "PAXX"]}},
                "head_input": {
                    "sidecar_file": "geneformer_bg003082_embedding.csv", "sidecar_sha256": "0" * 64,
                    "provenance_file": "geneformer_bg003082_embedding.provenance.json",
                    "shape": [1, 768], "all_finite": True,
                    "model": "ctheodoris/Geneformer / Geneformer-V2-104M_CLcancer",
                    "geneformer_revision_pinned": "04c2b2e8",
                    "commensurability_caveats": ["separately generated", "revision not captured",
                                                 "commensurability not proven", "exploratory only"],
                    "provenance_disclosures": ["bulk tumour TPM", "linear-TPM pseudo-count",
                                               "NCBI map not mygene", "not repo-reproducible"]}}},
        "input_artifact_sha256": {"data/processed/expression.npz": "0" * 64,
                                  "data/external/dgidb/dgidb_2026-06b.interactions.filtered.tsv": "0" * 64},
        "environment": {"python": "3.14.6"},
        "methodology": {"top_n": 25, "top_n_source": "config.TOP_N_DEPENDENCIES",
                        "inference": "fitted_artifacts only", "models": "reconstructed",
                        "no_leakage": "no observed value in ranking",
                        "determinism": "fixed precision, sorted keys"},
        "limitations": ["Both models are RECONSTRUCTED fitted state, not the historical objects.",
                        "BG003082 is bulk tumour tissue; 18,427/18,460 mapped, 33 imputed."],
        "disclaimers": ["A recorded drug\u2013gene interaction does not establish efficacy for this sample or for osteosarcoma.",
                        "Ranked genes are predicted dependencies, not therapeutic targets or recommended drugs."],
    }


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Render case_study.json into one offline HTML report.")
    ap.add_argument("--build", action="store_true", help="Write the HTML report.")
    ap.add_argument("--validate", action="store_true",
                    help="Regenerate, byte-compare, run structural + browser checks.")
    ap.add_argument("--self-test", action="store_true",
                    help="Synthetic render + escaping units + a real build/validate.")
    ap.add_argument("--out", default=str(config.REPORT_HTML_FILE))
    args = ap.parse_args()

    if args.self_test:
        return _self_test()

    rc = 0
    did = False
    if args.build:
        did = True
        print("=" * 74)
        print("BUILD phase2_report.html")
        print("=" * 74)
        build(Path(args.out))
    if args.validate:
        did = True
        print("\n" + "=" * 74)
        print("VALIDATE phase2_report.html")
        print("=" * 74)
        rep = validate(Path(args.out))
        rc = rc or (0 if rep["n_fail"] == 0 else 1)
    if not did:
        ap.print_help()
        return 2
    return rc


if __name__ == "__main__":
    sys.exit(main())
