from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from database.models import Base, engine
from routes.auth import router as auth_router
from routes.data import router as data_router  # Import data router

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

app.include_router(auth_router)
app.include_router(data_router)  # Register data routes

@app.get("/")
def read_root():
    return {"message": "Welcome to Abyssal Insights"}

