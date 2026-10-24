MolDockPipe offline molecular viewer

3Dmol.js version: 2.5.5
Source: https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js
SHA256: f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016
Upstream license: LICENSE-3Dmol.txt
Bundled dependency notices: 3Dmol-min.js.LICENSE.txt

viewer.html and viewer.js are MolDockPipe's local frontend. No CDN, network
structure lookup, or local HTTP server is used at runtime. Python supplies PDB
and SDF strings over QWebChannel rather than embedding them in the HTML URL.
