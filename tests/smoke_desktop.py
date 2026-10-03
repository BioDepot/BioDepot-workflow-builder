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
import os
from pathlib import Path
import pickle
import urllib.request
import xml.etree.ElementTree as ET
paths = glob.glob("/tmp/bwb-desktop-cache/**/widget-registry.pck", recursive=True)
assert len(paths) == 1, paths
with open(paths[0], "rb") as handle:
    registry = pickle.load(handle)
workflow = os.environ.get("STARTING_WORKFLOW")
if workflow:
    expected = {node.attrib["qualified_name"] for node in ET.parse(workflow).findall("./nodes/node")}
    discovered = {widget.qualified_name for widget in registry.widgets()}
    assert expected <= discovered, expected - discovered
else:
    assert len(registry.widgets()) == 34, len(registry.widgets())
canvases = []
for path in Path("/proc").glob("[0-9]*/cmdline"):
    try:
        command = path.read_bytes().decode().strip("\0").split("\0")
    except (OSError, UnicodeError):
        continue
    if not any(arg.endswith("/orange-canvas") for arg in command):
        continue
    if workflow and (workflow not in command or "__init" in command):
        continue
    canvases.append(path.parent.name)
assert len(canvases) == 1, canvases
assert urllib.request.urlopen("http://127.0.0.1:6080/", timeout=3).status == 200
print("PASS: desktop canvas pid {} with {}; registry and web endpoint ready".format(
    canvases[0], workflow or "34 standard widgets"))
'''


def run_desktop(args, workflow=None):
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
        command = [
            "docker", "run", "--rm", "-d", "--name", name,
            "-v", str(data) + ":/data",
            "-v", "/var/run/docker.sock:/var/run/docker.sock",
            "-e", "XDG_CACHE_HOME=/tmp/bwb-desktop-cache"]
        if workflow:
            command += ["-e", "STARTING_WORKFLOW=" + workflow]
        cid = subprocess.check_output(command + [args.image],
            universal_newlines=True).strip()
        try:
            deadline = time.monotonic() + args.startup_timeout
            previous_probe = None
            while time.monotonic() < deadline:
                result = subprocess.run(
                    ["docker", "exec", "-e", "QT_QPA_PLATFORM=offscreen",
                     cid, "python3", "-c", PROBE], stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, universal_newlines=True,
                    timeout=args.probe_timeout)
                if result.returncode == 0 and result.stdout == previous_probe:
                    # Two probes must see the same relaunched process alive.
                    assert log_sentinel.read_text() == "keep"
                    assert share_sentinel.read_text() == "keep"
                    print(result.stdout, end="", flush=True)
                    return
                previous_probe = result.stdout if result.returncode == 0 else None
                time.sleep(2)
            logs = subprocess.check_output(["docker", "logs", "--tail", "120", cid],
                                           stderr=subprocess.STDOUT,
                                           universal_newlines=True)
            desktop_logs = subprocess.run(
                ["docker", "exec", cid, "sh", "-c",
                 "tail -n 100 /var/log/supervisor/fluxbox* /var/log/web.log 2>/dev/null"],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                universal_newlines=True, timeout=10)
            processes = subprocess.check_output(["docker", "top", cid],
                                                universal_newlines=True)
            raise AssertionError((result.stdout, logs, desktop_logs.stdout, processes))
        finally:
            subprocess.run(["docker", "stop", "-t", "5", cid], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, timeout=15)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image")
    parser.add_argument("--startup-timeout", type=float, default=90,
                        help="seconds to wait for each desktop (allow more under emulation)")
    parser.add_argument("--probe-timeout", type=float, default=10,
                        help="seconds allowed for one registry/web probe")
    args = parser.parse_args()
    run_desktop(args)
    run_desktop(args, "/workflows/Demo_kallisto/Demo_kallisto.ows")


if __name__ == "__main__":
    main()
