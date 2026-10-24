/* Local-only 3Dmol frontend. No structure text is evaluated as JavaScript. */
"use strict";
new QWebChannel(qt.webChannelTransport, function (channel) {
    const bridge = channel.objects.bridge;
    let viewer, receptor, ligand, resetView, token = 0;
    let contacts = [];
    const contactStyles = {
        hbond: ["H-bond", "#0072B2"], pi_stack: ["π–π", "#7B3294"],
        cation_pi: ["Cation–π", "#D55E00"], salt_bridge: ["Salt", "#009E73"],
        hydrophobic: ["Hydrophobic", "#777777"]
    };
    const host = document.getElementById("molecule");
    function style(message) {
        receptor.setStyle({}, {});
        if (message.show_receptor) {
            if (message.representation === "Sticks") {
                receptor.setStyle({}, {stick: {colorscheme: "Jmol", radius: 0.12}});
            } else {
                receptor.setStyle({}, {cartoon: {color: "spectrum"}});
                receptor.setStyle({hetflag: true}, {stick: {colorscheme: "Jmol"}, sphere: {scale: 0.3}});
            }
        }
        ligand.setStyle({}, {stick: {colorscheme: "Jmol", radius: 0.22}});
        viewer.removeAllShapes();
        viewer.removeAllLabels();
        const active = new Set(message.interaction_types || []);
        const serials = new Set();
        const labelOffsets = new Map();
        if (message.show_receptor) for (const contact of contacts) {
            if (!active.has(contact.kind)) continue;
            const [name, color] = contactStyles[contact.kind];
            for (const serial of contact.protein_serials) serials.add(serial);
            // Cylinder dashes stay legible on Windows drivers that force WebGL
            // lines to one pixel, including in the higher-resolution export.
            const vector = {x: contact.end.x - contact.start.x, y: contact.end.y - contact.start.y,
                            z: contact.end.z - contact.start.z};
            const length = Math.hypot(vector.x, vector.y, vector.z);
            const point = t => ({x: contact.start.x + vector.x * t / length,
                                 y: contact.start.y + vector.y * t / length,
                                 z: contact.start.z + vector.z * t / length});
            for (let d = 0; d < length; d += 0.30) {
                viewer.addCylinder({start: point(d), end: point(Math.min(d + 0.18, length)),
                                    color, radius: 0.035, fromCap: 1, toCap: 1});
            }
            if (message.interaction_labels) {
                const candidate = contact.detail.includes("candidate") ? " candidate" : "";
                const row = labelOffsets.get(contact.residue) || 0;
                labelOffsets.set(contact.residue, row + 1);
                viewer.addLabel(`${contact.residue} · ${name}${candidate} ${contact.distance.toFixed(2)} Å`, {
                    position: {x: (contact.start.x + contact.end.x) / 2,
                               y: (contact.start.y + contact.end.y) / 2,
                               z: (contact.start.z + contact.end.z) / 2},
                    screenOffset: {x: 8, y: row * 22},
                    fontSize: 12, fontColor: color, backgroundColor: "white", backgroundOpacity: 0.85,
                    borderThickness: 0, inFront: true
                });
            }
        }
        // Make contacted side chains visible even when the protein is a cartoon.
        const atoms = receptor.selectedAtoms({});
        const residueKey = atom => `${atom.chain}|${atom.resi}|${atom.icode || ""}`;
        const residues = new Set(atoms.filter(atom => serials.has(atom.serial)).map(residueKey));
        const highlight = atoms.filter(atom => residues.has(residueKey(atom))).map(atom => atom.serial);
        if (highlight.length) viewer.addStyle({model: receptor.getID(), serial: highlight},
                                             {stick: {colorscheme: "Jmol", radius: 0.13}});
    }
    try {
        if (typeof $3Dmol === "undefined") throw new Error("The bundled 3Dmol library is missing");
        const probe = document.createElement("canvas");
        if (!(probe.getContext("webgl2") || probe.getContext("webgl"))) throw new Error("WebGL is unavailable");
        viewer = $3Dmol.createViewer(host, {backgroundColor: "white", antialias: true});
        host.addEventListener("webglcontextlost", function (event) {
            event.preventDefault();
            bridge.error(-1, "The WebGL graphics context was lost");
        }, true);
        new ResizeObserver(() => {viewer.resize(); viewer.render();}).observe(host);
    } catch (error) {
        bridge.error(-1, String(error.message || error));
        return;
    }
    bridge.command.connect(function (json) {
        let message;
        try {
            message = JSON.parse(json);
            if (message.kind === "clear" || message.kind === "load") token = message.token;
            if (message.token !== token) return;
            if (message.kind === "clear") {
                viewer.clear(); receptor = ligand = resetView = null; contacts = []; viewer.render();
            } else if (message.kind === "load") {
                viewer.clear();
                contacts = message.interactions || [];
                receptor = viewer.addModel(message.receptor_pdb, "pdb");
                ligand = viewer.addModel(message.ligand_sdf, "sdf");
                if (!receptor.selectedAtoms({}).length || !ligand.selectedAtoms({}).length) throw new Error("The structures contain no displayable atoms");
                style(message);
                viewer.zoomTo({model: ligand.getID()});
                viewer.zoom(0.65);
                resetView = viewer.getView().slice();
                viewer.render();
                bridge.rendered(token);
            } else if (message.kind === "style" && receptor && ligand) {
                style(message); viewer.render();
            } else if (message.kind === "focus" && ligand) {
                viewer.zoomTo({model: ligand.getID()}); viewer.render();
            } else if (message.kind === "reset" && resetView) {
                viewer.setView(resetView.slice()); viewer.render();
            } else if (message.kind === "export" && ligand) {
                const view = viewer.getView().slice();
                const width = host.clientWidth, height = host.clientHeight;
                // 3Dmol already supersamples on ordinary and high-DPI displays.
                // Compensate so the file is exactly 2x the CSS viewport, not 4x.
                const pixelRatio = host.querySelector("canvas").width / width;
                let png;
                try {
                    viewer.setWidth(width * 2 / pixelRatio); viewer.setHeight(height * 2 / pixelRatio);
                    viewer.setView(view); viewer.render();
                    png = viewer.pngURI();
                } finally {
                    viewer.setWidth(width); viewer.setHeight(height);
                    viewer.setView(view); viewer.render();
                }
                bridge.png(token, png);
            }
        } catch (error) {
            bridge.error(message ? message.token : token, String(error.message || error));
        }
    });
    bridge.ready();
});
