# from fastapi import APIRouter, Depends, HTTPException, status
# from pydantic import BaseModel, EmailStr, Field
# from sqlalchemy.orm import Session
# from passlib.context import CryptContext
# from fastapi.security import OAuth2PasswordRequestForm
# from jose import JWTError, jwt
# from datetime import datetime, timedelta
# import os

# router = APIRouter(prefix="/auth", tags=["Authentication"])

# pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

# SECRET_KEY = os.getenv("JWT_SECRET_KEY", "supersecretkey")
# ALGORITHM = "HS256"
# ACCESS_TOKEN_EXPIRE_MINUTES = 30

# class UserCreate(BaseModel):
#     email: EmailStr
#     password: str = Field(..., min_length=8, description="Password must be at least 8 characters long")
#     full_name: str

# from database.models import User, get_db

# def create_access_token(data: dict):
#     to_encode = data.copy()
#     expire = datetime.utcnow() + timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
#     to_encode.update({"exp": expire})
#     # jwt.encode returns a string token
#     encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
#     return encoded_jwt

# @router.post("/signup", status_code=status.HTTP_201_CREATED)
# def signup(user_in: UserCreate, db: Session = Depends(get_db)):
#     existing = db.query(User).filter(User.email == user_in.email).first()
#     if existing:
#         raise HTTPException(status_code=400, detail="Email already registered")
#     hashed_password = pwd_context.hash(user_in.password)
#     user = User(email=user_in.email, hashed_password=hashed_password, full_name=user_in.full_name)
#     db.add(user)
#     db.commit()
#     db.refresh(user)
#     return {"msg": "User created successfully", "user_id": user.id}

# @router.post("/login")
# def login(form_data: OAuth2PasswordRequestForm = Depends(), db: Session = Depends(get_db)):
#     user = db.query(User).filter(User.email == form_data.username).first()
#     if not user or not pwd_context.verify(form_data.password, user.hashed_password):
#         raise HTTPException(
#             status_code=status.HTTP_401_UNAUTHORIZED,
#             detail="Incorrect email or password",
#             headers={"WWW-Authenticate": "Bearer"},
#         )
#     access_token = create_access_token({"user_id": user.id})
#     return {"access_token": access_token, "token_type": "bearer"}

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from passlib.context import CryptContext
from fastapi.security import OAuth2PasswordRequestForm
from jose import JWTError, jwt
from datetime import datetime, timedelta
import os
from schemas.schemas import UserCreate, UserResponse, Token, TokenData

router = APIRouter(prefix="/auth", tags=["Authentication"])

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

SECRET_KEY = os.getenv("JWT_SECRET_KEY", "supersecretkey")
ALGORITHM = "HS256"
ACCESS_TOKEN_EXPIRE_MINUTES = 30

# Remove the duplicate UserCreate definition here
# class UserCreate(BaseModel):
#     email: EmailStr
#     password: str = Field(..., min_length=8, description="Password must be at least 8 characters long")
#     full_name: str

from database.models import User, get_db

def create_access_token(data: dict):
    to_encode = data.copy()
    expire = datetime.utcnow() + timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    to_encode.update({"exp": expire})
    # jwt.encode returns a string token
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt

@router.post("/signup", status_code=status.HTTP_201_CREATED)
def signup(user_in: UserCreate, db: Session = Depends(get_db)):
    existing = db.query(User).filter(User.email == user_in.email).first()
    if existing:
        raise HTTPException(status_code=400, detail="Email already registered")
    hashed_password = pwd_context.hash(user_in.password)
    user = User(email=user_in.email, hashed_password=hashed_password, full_name=user_in.full_name)
    db.add(user)
    db.commit()
    db.refresh(user)
    return {"msg": "User created successfully", "user_id": user.id}

@router.post("/login")
def login(form_data: OAuth2PasswordRequestForm = Depends(), db: Session = Depends(get_db)):
    user = db.query(User).filter(User.email == form_data.username).first()
    if not user or not pwd_context.verify(form_data.password, user.hashed_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect email or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token = create_access_token({"user_id": user.id})
    return {"access_token": access_token, "token_type": "bearer"}