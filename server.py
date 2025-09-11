from flask import Flask, request, redirect, send_from_directory, flash, url_for
import os
import subprocess
import pathlib

app = Flask(__name__)
app.secret_key = "replace-with-your-secret"  
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route("/<path:filename>")
def static_files(filename):
    return send_from_directory(".", filename)

@app.route("/", methods=["GET"])
def index():
    return redirect("/upload.html")

@app.route("/process", methods=["POST"])
def process():
    if "file" not in request.files:
        flash("No file part")
        return redirect("/upload.html")
    file = request.files["file"]
    if file.filename == "":
        flash("No selected file")
        return redirect("/upload.html")

    safe_name = pathlib.Path(file.filename).name
    upload_path = os.path.join(UPLOAD_FOLDER, safe_name)
    file.save(upload_path)

    try:
    
        subprocess.run(["python", "loading.py", upload_path], check=True)
    except subprocess.CalledProcessError as e:

        print("Pipeline failed:", e)
        flash("Server error while processing your file.")
        return redirect("/upload.html")

    return redirect("/result.html")

if __name__ == "__main__":

    app.run(debug=True, host="0.0.0.0", port=5000)
