"""Local-only ASGI service; scientific work is isolated from request handling."""

import asyncio
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from contextlib import asynccontextmanager
from functools import partial
from multiprocessing import get_context
from pathlib import Path
from datetime import datetime
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import Response, FileResponse
from fastapi.staticfiles import StaticFiles
from .models import Scene, Orbit, ExportRequest
from .service import compute, load_weather, export_report


@asynccontextmanager
async def lifespan(app):
    app.state.pool = ProcessPoolExecutor(max_workers=1, mp_context=get_context("spawn"))
    app.state.heavy_pool = ProcessPoolExecutor(
        max_workers=1, mp_context=get_context("spawn")
    )
    app.state.heavy_lock = asyncio.Lock()
    app.state.io_pool = ThreadPoolExecutor(max_workers=1)
    try:
        yield
    finally:
        app.state.pool.shutdown(wait=False, cancel_futures=True)
        app.state.heavy_pool.shutdown(wait=False, cancel_futures=True)
        app.state.io_pool.shutdown(wait=False, cancel_futures=True)


app = FastAPI(title="digitalQQ", version="0.1.0", lifespan=lifespan)


async def run(fn, *args, io=False, heavy=False, timeout=120, request=None):
    try:
        pool = (
            app.state.io_pool
            if io
            else app.state.heavy_pool if heavy else app.state.pool
        )

        async def submit():
            if request is not None and await request.is_disconnected():
                raise HTTPException(499, "请求已被更新或关闭")
            return await asyncio.get_running_loop().run_in_executor(
                pool, partial(fn, *args)
            )

        async def serialized():
            async with app.state.heavy_lock:
                return await submit()

        return await asyncio.wait_for(serialized() if heavy else submit(), timeout)
    except asyncio.TimeoutError:
        raise HTTPException(504, "计算或数据读取超时，请减小精度或检查数据源")
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        raise HTTPException(422, str(exc)) from exc


@app.get("/api/health")
def health():
    return {"status": "ok", "version": "0.1.0"}


@app.post("/api/attitude")
async def attitude(scene: Scene):
    return await run(compute, "attitude", scene.model_dump(mode="json"))


@app.post("/api/solve-pointing")
async def solve(scene: Scene):
    return await run(compute, "solve", scene.model_dump(mode="json"))


@app.post("/api/environment")
async def environment(
    scene: Scene,
    request: Request,
    part: str = Query("all", pattern="^(all|boundaries)$"),
):
    return await run(
        compute,
        "boundaries" if part == "boundaries" else "environment",
        scene.model_dump(mode="json"),
        heavy=part == "all",
        request=request,
        timeout=300 if scene.quality == "fine" else 120,
    )


@app.post("/api/scene/validate")
def validate_scene(scene: Scene):
    return scene.model_dump(mode="json")


@app.post("/api/orbit/validate")
def orbit(orbit: Orbit):
    return orbit.model_dump(mode="json")


@app.get("/api/space-weather")
async def weather(time: datetime, download: bool = False):
    if time.tzinfo is None:
        raise HTTPException(422, "时间必须包含时区")
    return await run(load_weather, time.isoformat(), download, io=True)


@app.post("/api/export")
async def export(req: ExportRequest, request: Request):
    data = await run(
        export_report,
        req.model_dump(mode="json"),
        timeout=300,
        heavy=True,
        request=request,
    )
    return Response(
        data,
        media_type="application/zip",
        headers={
            "Content-Disposition": 'attachment; filename="digitalQQ_observation.zip"'
        },
    )


DIST = Path(__file__).resolve().parents[1] / "frontend" / "dist"
if DIST.exists():
    app.mount("/assets", StaticFiles(directory=DIST / "assets"), name="assets")

    @app.get("/")
    def index():
        return FileResponse(DIST / "index.html")
