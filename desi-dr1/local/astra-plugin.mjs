// In-project shim: re-export the compiled ASTRA plugin from the sibling MySTRA
// build so myst.yml can reference a plugin path inside the project.
// Assumes the repo layout  repo/{astra-theme,MySTRA}  and that MySTRA has been
// built (`npm run build` at the MySTRA repo root → dist/index.js).
export { default } from '../../../../MySTRA/dist/index.js';
