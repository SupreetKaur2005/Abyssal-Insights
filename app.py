from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from database.models import Base, engine
from routes.auth import router as auth_router

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

@app.get("/")
def read_root():
    return {"message": "Welcome to Abyssal Insights"}

