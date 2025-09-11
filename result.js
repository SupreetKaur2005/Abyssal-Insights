function loadCSVPreview(filePath, tableId, maxRows = 10) {
  fetch(filePath)
    .then(response => response.text())
    .then(text => {
      const rows = text.trim().split("\n").map(r => r.split(","));
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
      console.error("Error loading CSV:", err);
    });
}

// Load previews for both CSVs
loadCSVPreview("known_sequences.csv", "knownTable");
loadCSVPreview("unknown_sequences.csv", "unknownTable");
