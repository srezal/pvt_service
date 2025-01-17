import uvicorn
from fastapi import FastAPI
from routes.pvt import router as pvt_router

app = FastAPI()

app.include_router(pvt_router)

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8001)