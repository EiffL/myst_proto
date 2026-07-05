# DESI DR1 BAO — astra-theme live example

This project renders the real DESI DR1 BAO analysis through **astra-theme**,
authored with the [MySTRA](https://github.com/LightconeResearch/MySTRA) plugin
(`@astra-spec/mystra`).

## Run it

```bash
myst start                             # → http://localhost:3000  (or pass --port 3111)
```

How it's wired (see `myst.yml`):
- `project.plugins:` — the bundled MySTRA release
  (`https://github.com/LightconeResearch/MySTRA/releases/download/v0.0.7/mystra.mjs`).
  It emits the neutral `astra-*` markers + the per-page resolved store.
- `site.template: https://github.com/Eiffl/astra-theme.git` — the theme (a
  standalone fork of @myst-theme/book with the ASTRA renderers compiled in).
  Point it at a local checkout instead while developing the theme.

The full rich experience renders here: inline hover preview cards, the interactive decision
panel, output figures with provenance drawers, finding/insight cards, registry tables, and live
value tokens.

---

The rest of this README is a quick tour of the **ASTRA authoring grammar** used
in this project. The authoritative guide is the
[MySTRA documentation](https://github.com/LightconeResearch/MySTRA/tree/main/docs).

## Authoring grammar

Every reference is a **dot-separated path** that mirrors `astra.yaml`
(collections = `inputs outputs decisions findings prior_insights analyses
universes`; a sub-analysis id may be written directly, `analyses.` implied).
Paths always resolve from the **root analysis**, on every page. One name —
`astra` — is both a role (inline) and a directive (block).

Block embeds (`:::{astra} <path>`):

```markdown
:::{astra} decisions.covariance_source
:::
:::{astra} outputs.bao_fit_plot
:::
:::{astra} findings.bao_detected_post_recon
:::
:::{astra} prior_insights.spline_broadband_fiducial
:::
:::{astra} inputs
:::                                   # full inputs registry table (root scope)
:::{astra} clustering.outputs
:::                                   # outputs registry, clustering sub-analysis
:::{astra} reconstruction
:::                                   # nav card to the sub-analysis page
```

Inline references (the `{astra}` role) — each renders inline as a **label**
(neutral text; a rich theme like `astra-theme` adds the kind glyph and a hover
**preview card** built from the resolved store).

```markdown
{astra}`decisions.covariance_source`              ◇ Decision
{astra}`outputs.hubble_diagram_plot`              ◆ Output
{astra}`findings.subpercent_alpha_iso_precision`  ● Finding
{astra}`prior_insights.recon_reduces_bao_damping` ◈ Prior insight
{astra}`reconstruction`                           ◐ Sub-analysis (bare path)
```

**Neutral by design.** The plugin emits *only* a semantic span + a join key — no
glyphs, colours, or card markup are baked into the AST. The theme builds the card
from the resolved store, keyed by `data.astra`:

```
span.astra-ref.astra-ref--<kind>[ .astra-ref--<subtype>]   children: [ text(label) ]
  data.astra: { kind, id, path }
```

**Display text** — MyST's `text <path>` override (the card still shows the
element's own label/claim):

```markdown
{astra}`few-fold precision loss <prior_insights.precision_loss_factor_three>`
```

**Numbered cross-references** — `{astra:ref}` (like `{ref}`; supports `%s`;
`{astra:numref}` is accepted as an alias):

```markdown
{astra:ref}`outputs.bao_fit_plot`            → Figure 1
{astra:ref}`see Fig. %s <outputs.bao_fit_plot>`
```

Inline **value interpolation** — never hard-type a measured number; pull it
live from a result product. The body is the path; the selection is expressed as
**role options** (MyST inline options):

```markdown
{astra:value col=DV_over_rd where="tracer=lrg3_elg1" pm=true}`outputs.bao_distance_table`  → 19.88 ± 0.17
{astra:value col=alpha1_mean where="tracer=lrg3_elg1 recon=Post" err=alpha1_std}`outputs.bao_alpha_values`
{astra:value col=alpha1_std where="tracer=elg1 recon=Pre"}`outputs.bao_alpha_values`       → 0.0696
{astra:value}`decisions.smoothing_radius`                                                  → the selected option
```

Options: `col=<column>`, `where="<key>=<val> …"` (row filters), `pm=true`
(append `± <col>_std`), `err=<column>` (explicit uncertainty column),
`sig=<N>` (significant figures). Metric outputs interpolate directly with no
options. `index.md` uses this for every number in the prose, so nothing is
hand-typed.

Extras:

- `:::{astra} outputs.<id>` **figures** carry the `output-<id>` anchor; **tables**
  render as a clean numbered `container[table]`. Findings accept `:compact:`
  (claim + notes + scope, no evidence), and outputs accept `:caption:`.
- A path that stops at a collection (`outputs`, `clustering.inputs`) renders the
  whole **registry**; children are addressable too, with the child collection
  implied (`decisions.<id>.<option>`, `findings.<id>.<evidence>`).
- **Citations**: `{astra:cite}` / `{astra:cite:t}` turn DOI-backed evidence into
  real author–year citations through MyST's bibliography pipeline.

**Scoping.** A sub-analysis is a path segment
(`reconstruction.decisions.algorithm`, `clustering.outputs.xi_multipoles_plot`)
and paths mean the same thing on every page. A page's *store scope* follows the
dotted-filename convention (`reconstruction.md` → the `reconstruction`
sub-analysis); cross-page links use plain MyST anchors (`[](#output-<id>)`).
