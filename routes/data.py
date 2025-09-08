from fastapi import APIRouter, Depends, File, UploadFile
from sqlalchemy.orm import Session

from dependencies.auth import get_current_user
from database.models import get_db, User

router = APIRouter(prefix="/data", tags=["Data Operations"])

@router.post("/upload")
async def upload_data(
    file: UploadFile = File(...),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Only authenticated users can upload files.
    """
    # Save the file to disk
    file_location = f"uploads/{current_user.id}_{file.filename}"
    with open(file_location, "wb") as buffer:
        buffer.write(await file.read())
    return {
        "msg": "Upload successful",
        "filename": file.filename,
        "path": file_location
    }

@router.post("/test")
async def test_data(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Only authenticated users can test their data.
    """
    # Placeholder for actual testing logic
    return {"status": "success", "detail": "Your data was tested."}

