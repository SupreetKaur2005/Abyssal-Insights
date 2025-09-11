# from fastapi import FastAPI, Depends, Request
# from fastapi.middleware.cors import CORSMiddleware
# from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse
# from fastapi.staticfiles import StaticFiles
# from database.models import Base, engine
# from routes.auth import router as auth_router
# from routes.data import router as data_router
# from routes.data import get_current_user, User
# from typing import Optional
# import os

# # Create tables
# Base.metadata.create_all(bind=engine)

# app = FastAPI(title="Abyssal Insights API")

# app.add_middleware(
#     CORSMiddleware,
#     allow_origins=["*"],
#     allow_credentials=True,
#     allow_methods=["*"],
#     allow_headers=["*"],
# )

# # Mount static files (HTML, CSS, JS, etc.)
# app.mount("/static", StaticFiles(directory="static"), name="static")

# # Include routers
# app.include_router(auth_router)
# app.include_router(data_router)

# async def get_optional_user(request: Request) -> Optional[User]:
#     """Try to get current user, return None if not authenticated"""
#     try:
#         return await get_current_user(request)
#     except:
#         return None

# @app.get("/", response_class=HTMLResponse)
# async def read_root(request: Request, current_user: Optional[User] = Depends(get_optional_user)):
#     """
#     Serve home page if authenticated, otherwise login page.
#     """
#     if current_user:
#         return FileResponse("static/home-page.html")
#     else:
#         return FileResponse("static/login_signup.html")

# @app.get("/me")
# def get_me(current_user: User = Depends(get_current_user)):
#     """
#     Get info about the current logged-in user (protected route).
#     """
#     return {"id": current_user.id, "email": current_user.email, "username": current_user.username}

# # Add a catch-all route to serve the login page for any other routes
# @app.get("/{full_path:path}")
# async def catch_all(request: Request, full_path: str, current_user: Optional[User] = Depends(get_optional_user)):
#     """
#     Catch all routes and redirect to appropriate page based on authentication.
#     """
#     if current_user:
#         return FileResponse("static/home-page.html")
#     else:
#         return FileResponse("static/login_signup.html")

# if __name__ == "__main__":
#     import uvicorn
#     uvicorn.run(app, host="0.0.0.0", port=8000)


# from fastapi import FastAPI, Depends, Request
# from fastapi.middleware.cors import CORSMiddleware
# from fastapi.responses import FileResponse, HTMLResponse
# from fastapi.staticfiles import StaticFiles
# from database.models import Base, engine
# from routes.auth import router as auth_router
# from routes.data import router as data_router
# from typing import Optional
# import os

# # Create tables
# Base.metadata.create_all(bind=engine)

# app = FastAPI(title="Abyssal Insights API")

# app.add_middleware(
#     CORSMiddleware,
#     allow_origins=["*"],
#     allow_credentials=True,
#     allow_methods=["*"],
#     allow_headers=["*"],
# )

# # Mount static files (HTML, CSS, JS, etc.)
# app.mount("/static", StaticFiles(directory="static"), name="static")

# # Include routers
# app.include_router(auth_router)
# app.include_router(data_router)

# @app.get("/", response_class=HTMLResponse)
# async def read_root():
#     """
#     Serve home page
#     """
#     return FileResponse("static/home-page.html")

# @app.get("/login", response_class=HTMLResponse)
# async def login_page():
#     """
#     Serve login page
#     """
#     return FileResponse("static/login.html")

# @app.get("/upload", response_class=HTMLResponse)
# async def upload_page():
#     """
#     Serve upload page
#     """
#     return FileResponse("static/upload.html")

# if __name__ == "__main__":
#     import uvicorn
#     uvicorn.run(app, host="127.0.0.1", port=8000)

from fastapi import FastAPI, Depends, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from database.models import Base, engine
from routes.auth import router as auth_router
from routes.data import router as data_router
from typing import Optional
import os

# Create tables
Base.metadata.create_all(bind=engine)

app = FastAPI(title="Abyssal Insights API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount current directory for static files
app.mount("/static", StaticFiles(directory="."), name="static")

# Include routers
app.include_router(auth_router)
app.include_router(data_router)

# Main pages
@app.get("/", response_class=HTMLResponse)
async def read_root():
    return FileResponse("home-page.html")

@app.get("/login", response_class=HTMLResponse)
async def login_page():
    return FileResponse("login.html")

@app.get("/upload", response_class=HTMLResponse)
async def upload_page():
    return FileResponse("upload.html")

@app.get("/interactive_edna_clusters", response_class=HTMLResponse)
async def interactive_clusters():
    return FileResponse("interactive_edna_clusters(Unknown_Sequences).html")

# Static file endpoints for all your files
@app.get("/{filename}.css")
async def get_css(filename: str):
    return FileResponse(f"{filename}.css")

@app.get("/{filename}.js")
async def get_js(filename: str):
    return FileResponse(f"{filename}.js")

@app.get("/{filename}.png")
async def get_png(filename: str):
    return FileResponse(f"{filename}.png")

@app.get("/{filename}.csv")
async def get_csv(filename: str):
    return FileResponse(f"{filename}.csv")

@app.get("/{filename}.joblib")
async def get_joblib(filename: str):
    return FileResponse(f"{filename}.joblib")

@app.get("/{filename}.html")
async def get_html(filename: str):
    return FileResponse(f"{filename}.html")

# Catch-all for any other file requests
@app.get("/{path:path}")
async def catch_all(path: str):
    if os.path.exists(path):
        return FileResponse(path)
    else:
        return {"error": "File not found"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)