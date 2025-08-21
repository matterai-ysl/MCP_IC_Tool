from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def hello():
    return {"msg": "Hello from ln2!"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=9528)