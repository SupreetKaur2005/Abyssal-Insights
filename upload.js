const fileInput = document.getElementById("fileInput");
const fileName = document.getElementById("fileName");
const uploadForm = document.getElementById("uploadForm");
const uploadResult = document.getElementById("uploadResult");

// Show selected file name
fileInput.addEventListener("change", () => {
  fileName.textContent = fileInput.files.length > 0 
    ? fileInput.files[0].name 
    : "Click here to select a file";
});

// Handle form submit
uploadForm.addEventListener("submit", (e) => {
  e.preventDefault();

  if (fileInput.files.length === 0) {
    uploadResult.textContent = "⚠️ Please select a file before uploading.";
    return;
  }

  const file = fileInput.files[0];
  
  // Just simulate upload result for now
  uploadResult.innerHTML = `
    ✅ File Uploaded Successfully! <br>
    <b>File Name:</b> ${file.name} <br>
    <b>Size:</b> ${(file.size / 1024).toFixed(2)} KB <br>
    <b>Type:</b> ${file.type || "Unknown"}
  `;
});
