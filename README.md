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

Every reference is a **path** that mirrors `astra.yaml` (`<collection>/<id>[/…]`,
collections = `inputs outputs decisions findings prior_insights analyses
universes`). One name — `astra` — is both a role (inline) and a directive (block).

Block embeds (`:::{astra} <path>`):

```markdown
:::{astra} decisions/covariance_source
:::
:::{astra} outputs/bao_fit_plot
:::
:::{astra} findings/bao_detected_post_recon
:::
:::{astra} prior_insights/spline_broadband_fiducial
:::
:::{astra} inputs
:::                                   # full inputs registry table (root scope)
:::{astra} clustering/outputs
:::                                   # outputs registry, clustering sub-analysis
:::{astra} reconstruction
:::                                   # nav card to the sub-analysis page
```

Inline references (the `{astra}` role) — each renders inline as a **label**
(neutral text; a rich theme like `astra-theme` adds the kind glyph and a hover
**preview card** built from the resolved store).

```markdown
{astra}`decisions/covariance_source`              ◇ Decision
{astra}`outputs/hubble_diagram_plot`              ▦ Output
{astra}`findings/subpercent_alpha_iso_precision`  ● Finding
{astra}`prior_insights/recon_reduces_bao_damping` ◈ Prior insight
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
{astra}`few-fold precision loss <prior_insights/precision_loss_factor_three>`
```

**Numbered cross-references** — `{astra:numref}` (like `{numref}`; supports `%s`):

```markdown
{astra:numref}`outputs/bao_fit_plot`            → Figure 1
{astra:numref}`see Fig. %s <outputs/bao_fit_plot>`
```

Inline **value interpolation** — never hard-type a measured number; pull it
live from a result product:

```markdown
{astra:value}`outputs/bao_distance_table col=DV_over_rd tracer=lrg3_elg1 ±`   → 19.88 ± 0.17
{astra:value}`outputs/bao_alpha_values col=alpha1_mean tracer=lrg3_elg1 recon=Post err=alpha1_std`
{astra:value}`outputs/bao_alpha_values col=alpha1_std tracer=elg1 recon=Pre`  → 0.0696
{astra:value}`decisions/smoothing_radius`                                     → the selected option
```

Grammar: `<path> [col=<column>] [<key>=<val> …] [±|pm] [err=<col>] [sig=N]`. It
reads the materialised CSV/JSON, filters rows by the `key=val` pairs, and renders
the cell (with `± std` via `±`/`err=`). `index.md` uses this for every number in
the prose, so nothing is hand-typed.

Extras:

- `:::{astra} outputs/<id>` **figures** carry the `output-<id>` anchor; **tables**
  render as a clean numbered `container[table]`. Findings accept `:compact:`
  (claim + notes + scope, no evidence), and outputs accept `:caption:`.
- A path that stops at a collection (`outputs`, `clustering/inputs`) renders the
  whole **registry**; children are addressable too (`decisions/<id>/options/<o>`,
  `findings/<id>/evidence/<e>`).

**Scoping.** Role/directive paths resolve from the **root analysis**; a
sub-analysis is a path segment (`reconstruction/decisions/algorithm`,
`clustering/outputs/xi_multipoles_plot`), so sub-analysis pages
(`reconstruction.md`, `clustering.md`) write the scoped path. The `#astra:<path>`
link scheme (e.g. `[](#astra:outputs/bao_fit_plot)`) resolves relative to the
current page.

> **Plugin version.** These conventions need a MySTRA build that postdates the
> unified-path-grammar redesign; the `v0.0.2` pin in `myst.yml` predates them.

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
