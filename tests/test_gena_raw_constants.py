import os
import subprocess
import sys


def _read_save_pdf(env):
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from gena.raw import constants; print(constants.save_pdf)",
        ],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    return int(result.stdout.strip())


def test_save_pdf_defaults_to_enabled():
    env = os.environ.copy()
    env.pop("GENA_SAVE_PDF", None)
    assert _read_save_pdf(env) == 1


def test_save_pdf_reads_environment_override():
    env = os.environ.copy()
    env["GENA_SAVE_PDF"] = "0"
    assert _read_save_pdf(env) == 0
