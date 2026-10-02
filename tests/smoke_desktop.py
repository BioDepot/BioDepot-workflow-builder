#!/usr/bin/env python3
"""Start the actual desktop stack in an isolated container without publishing ports."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import time
import uuid

PROBE = r'''
import glob
import pickle
import urllib.request
paths = glob.glob("/tmp/bwb-desktop-cache/**/widget-registry.pck", recursive=True)
assert len(paths) == 1, paths
with open(paths[0], "rb") as handle:
    registry = pickle.load(handle)
assert len(registry.widgets()) == 34, len(registry.widgets())
assert urllib.request.urlopen("http://127.0.0.1:6080/", timeout=3).status == 200
print("PASS: actual desktop registry has 34 widgets; web endpoint responds")
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image")
    args = parser.parse_args()
    name = "bwb-desktop-smoke-" + uuid.uuid4().hex
    with tempfile.TemporaryDirectory(prefix="bwb-desktop-smoke-") as folder:
        data = Path(folder) / "data"
        data.mkdir()
        for path in (data, data / ".bwb", data / ".bwbshare"):
            path.mkdir(exist_ok=True)
            path.chmod(0o777)
        log_sentinel = data / ".bwb" / "other-job"
        share_sentinel = data / ".bwbshare" / "other-job"
        log_sentinel.write_text("keep")
        share_sentinel.write_text("keep")
        cid = subprocess.check_output([
            "docker", "run", "--rm", "-d", "--name", name,
            "-v", str(data) + ":/data",
            "-v", "/var/run/docker.sock:/var/run/docker.sock",
            "-e", "XDG_CACHE_HOME=/tmp/bwb-desktop-cache", args.image],
            universal_newlines=True).strip()
        try:
            deadline = time.monotonic() + 90
            while time.monotonic() < deadline:
                result = subprocess.run(
                    ["docker", "exec", "-e", "QT_QPA_PLATFORM=offscreen",
                     cid, "python3", "-c", PROBE], stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, universal_newlines=True, timeout=10)
                if result.returncode == 0:
                    processes = subprocess.check_output(["docker", "top", cid],
                                                        universal_newlines=True)
                    assert "orange-canvas" in processes, processes
                    assert log_sentinel.read_text() == "keep"
                    assert share_sentinel.read_text() == "keep"
                    print(result.stdout, end="", flush=True)
                    return
                time.sleep(2)
            logs = subprocess.check_output(["docker", "logs", "--tail", "120", cid],
                                           stderr=subprocess.STDOUT,
                                           universal_newlines=True)
            raise AssertionError((result.stdout, logs))
        finally:
            subprocess.run(["docker", "stop", "-t", "5", cid], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, timeout=15)


if __name__ == "__main__":
    main()
