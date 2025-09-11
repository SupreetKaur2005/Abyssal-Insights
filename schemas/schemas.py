from pydantic import BaseModel
from typing import Optional
from datetime import datetime

# User schemas
class UserBase(BaseModel):
    email: str
    full_name: str  # Changed from username to full_name

class UserCreate(UserBase):
    password: str

class UserLogin(BaseModel):
    email: str
    password: str

class UserResponse(UserBase):
    id: int
    created_at: datetime
    
    class Config:
        from_attributes = True

# Token schemas
class Token(BaseModel):
    access_token: str
    token_type: str

class TokenData(BaseModel):
    id: Optional[str] = None

# Data schemas (for your main application data)
class DataBase(BaseModel):
    title: str
    content: str
    category: str

class DataCreate(DataBase):
    pass

class DataResponse(DataBase):
    id: int
    created_at: datetime
    user_id: int
    
    class Config:
        from_attributes = True

# Update schemas
class UserUpdate(BaseModel):
    email: Optional[str] = None  # Changed from EmailStr to str
    username: Optional[str] = None

class DataUpdate(BaseModel):
    title: Optional[str] = None
    content: Optional[str] = None
    category: Optional[str] = None