// ---- CSV preview loader ----
function loadCSVPreview(filePath, tableId, maxRows = 10) {
  fetch(filePath)
    .then(response => {
      if (!response.ok) throw new Error("Not found");
      return response.text();
    })
    .then(text => {
      const rows = text.trim().split("\n").map(r => {
        // simple CSV split: this assumes no commas inside fields
        return r.split(",");
      });
      const tbody = document.querySelector(`#${tableId} tbody`);
      tbody.innerHTML = "";
      rows.slice(1, maxRows + 1).forEach(row => {
        const tr = document.createElement("tr");
        row.forEach(cell => {
          const td = document.createElement("td");
          td.textContent = cell;
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      });
    })
    .catch(err => {
      // file missing or parse error — leave table empty
      console.warn(`CSV ${filePath} not loaded:`, err);
    });
}

async function checkAndShowVisuals() {
  const pngPath = "edna_clusters_2d_Unknown_Sequences.png";
  const htmlPath = "interactive_edna_clusters(Unknown_Sequences).html";

  try {
    // Check PNG
    let resp = await fetch(pngPath, { method: "HEAD" });
    if (resp.ok) {
      document.getElementById("clusterImg").src = pngPath;
      document.getElementById("clusterImg").style.display = "block";

      const dl = document.getElementById("downloadImg");
      dl.href = pngPath;
      dl.style.display = "inline-block";
    }

    // Check interactive HTML
    let resp2 = await fetch(htmlPath, { method: "HEAD" });
    if (resp2.ok) {
      const link = document.getElementById("interactiveLink");
      link.href = htmlPath;
      link.style.display = "inline-block";
    }

    document.getElementById("visualBlock").style.display = "block";
  } catch (e) {
    console.warn("Visual check error:", e);
  }
}

// Run loaders on page load
loadCSVPreview("known_sequences.csv", "knownTable", 10);
loadCSVPreview("unknown_sequences.csv", "unknownTable", 10);
checkAndShowVisuals();
