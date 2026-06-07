# DESI DR1 BAO — astra-theme live example

This project renders the real DESI DR1 BAO analysis through **astra-theme** on the
actual `myst start` runtime (Phase 1: the Vellum design system as a CSS layer over
stock book-theme — see `../../STATUS.md`).

## Run it

```bash
# 1. build the MySTRA plugin once (sibling repo): cd ../../../MySTRA && npm i && npm run build
# 2. build the astra-theme Remix app once:        cd ../.. && npm install && npm run build
# 3. start MyST from this project:
cd examples/desi-dr1
myst start                             # → http://localhost:3000  (or pass --port 3111)
```

How it's wired (just `myst.yml` + the plugin shim):
- `project.plugins: [local/astra-plugin.mjs]` — re-exports the compiled MySTRA plugin from
  the sibling `../../../MySTRA/dist/index.js` (it emits the neutral `astra-*` markers + store).
- `site.template: ../..` — this repo IS the theme (a standalone fork of @myst-theme/book with
  the ASTRA renderers compiled in). MyST resolves a local template by path; once astra-theme is
  published to the registry the one-liner `template: astra-theme` works for everyone. The theme
  carries the ASTRA design system + renderers itself — no per-project `style` option.

The full rich experience renders here: inline hover preview cards, the interactive decision
panel, output figures with provenance drawers, finding/insight cards, registry tables, and live
value tokens. Edit the design in `../../styles/astra.css` or renderers in `../../app/astra/`,
re-run `npm run build` at the repo root, restart `myst start`.

---

The rest of this README is the **ASTRA authoring grammar** — how the Markdown in
this project imports and cites ASTRA components via the `@astra-spec/mystra`
plugin (the same neutral markers `astra-theme` renders richly).

## Authoring grammar

Block "import" directives (one component, by id):

```markdown
:::{astra:decision} covariance_source
:::
:::{astra:output} bao_fit_plot
:::
:::{astra:finding} bao_detected_post_recon
:::
:::{astra:prior-insight} spline_broadband_fiducial
:::
:::{astra:inputs}
:::                                   # full inputs registry table (root scope)
:::{astra:outputs} clustering
:::                                   # outputs table, clustering sub-analysis
:::{astra:subanalysis} reconstruction
:::                                   # nav card to the sub-analysis page
```

Inline "cite" roles — each renders inline as a **label** (neutral text; a rich
theme like `astra-theme` adds the kind glyph and a hover **preview card** built
from the resolved store: claim / rationale / selected option / exact quote +
citation — focused, not a whole-table dump).

```markdown
{astra:decision}`covariance_source`           ◇ Decision  — selected option + rationale
{astra:output}`hubble_diagram_plot`           ▦ Output    — type + description
{astra:finding}`subpercent_alpha_iso_precision`  ● Finding — claim + scope
{astra:prior-insight}`recon_reduces_bao_damping`  ◈ Prior insight — claim + quote + citation
{astra:analysis}`reconstruction`              ◐ Sub-analysis — summary + counts
```

**Neutral by design.** The plugin emits *only* a semantic span + a join key — no
glyphs, colours, card markup, or inline styles are baked into the AST. The theme
builds the card from the resolved store, keyed by `data.astra`:

```
span.astra-ref.astra-ref--<kind>[ .astra-ref--<subtype>]   children: [ text(label) ]
  data.astra: { kind, id, path }
```

**Optional display text** — append `|text` to set the inline label (the card
still shows the element's own label/claim); the id fallback is humanised:

```markdown
{astra:prior-insight}`precision_loss_factor_three|few-fold precision loss`
```

Inline **value interpolation** — never hard-type a measured number; pull it
live from a result product:

```markdown
{astra:value}`bao_distance_table tracer=lrg3_elg1 col=DV_over_rd pm`   → 19.88 ± 0.17
{astra:value}`bao_alpha_values tracer=lrg3_elg1 recon=Post col=alpha1_mean err=alpha1_std`
{astra:value}`bao_alpha_values tracer=elg1 recon=Pre col=alpha1_std`   → 0.0696
```

Grammar: `<output-path> col=<column> [<key>=<val> …] [pm] [err=<col>] [sig=N]`.
It reads the materialised CSV/JSON, filters rows by the `key=val` pairs, and
renders the cell (with `± std` via `pm`/`err=`). `index.md` uses this for every
number in the prose, so nothing is hand-typed.

Each interpolated number renders as a glyph token (▦/▤ by source type) with a
focused hover card naming its **source product, column, and row filters** — so
the reader sees at a glance that a number is sourced data and exactly where it
comes from, without a whole-table overlay.

Extras:

- `:::{astra:output}` **figures** carry the `output-<id>` anchor; **tables**
  render as a clean numbered `container[table]` (not a collapsible). Both get a
  collapsed **ASTRA provenance** disclosure (type, upstream products, decisions,
  recipe command) — a traceable, first-class representation of the product.
- `:::{astra:finding} <id>` accepts `:compact:` to render claim + notes + scope
  only (no evidence figures) — used for the register's hover targets.
- Reference a figure/table by number with a plain link: `[](#output-bao_fit_plot)`.

**Scoping.** A path is `<id>` (root analysis) or `<sub>.<id>` (sub-analysis),
e.g. `reconstruction.algorithm`, `clustering.xi_multipoles_plot`. Sub-analysis
pages (`reconstruction.md`, `clustering.md`) use the scoped prefix.

The plugin reads `astra.yaml` once (cached) and renders each component via the
per-component helpers in `../src/transform/`. Source: `../src/index.ts` (the
plugin *is* the package entry; compiled to `../dist/index.js`, re-exported by
`local/astra-plugin.mjs`).

## Known limitations (prototype scope)

- **`astra.yaml` live reload.** `myst start` watches Markdown, not `astra.yaml`,
  and the plugin caches the parse — edit `astra.yaml` ⇒ restart the server (or
  touch a `.md`). Follow-up: cache invalidation / add `astra.yaml` to a watch.
- **Citations are plain DOI links.** Evidence DOIs render as `doi.org` links —
  no author–year label and no reference list yet. Wiring a MyST-native
  bibliography is a follow-up (see `../followup.md`).
- **Cross-scope insight references.** A sub-analysis decision whose options list
  `insights:` in the *parent* scope emits cross-references
  (`prior_insight-…`) that only resolve if that insight is also imported on the
  page. On `reconstruction.md` these show as warnings (the text still renders).
- **`prior-insight` rendering.** The plugin renders a prior insight as a
  `seealso` admonition (a custom `container` kind with no caption is rejected by
  the stock theme). A bundled theme could render a richer treatment.
