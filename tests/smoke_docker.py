#!/usr/bin/env python3
"""Opt-in real Docker integration test of an already-built Bwb image.

Usage: python3 tests/smoke_docker.py IMAGE
Only temporary test directories and --rm test containers are used. The Bwb
container runs without root and starts real sibling containers via the socket.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile


IN_CONTAINER = r'''
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import uuid

sys.path.insert(0, "/coreutils")
from DockerClient import DockerClient

client = DockerClient("unix:///var/run/docker.sock", "local")
assert os.environ["BWBSHARE"] == "/data/.bwbshare", os.environ["BWBSHARE"]
assert os.environ["BWBHOSTSHARE"] == os.environ["EXPECTED_HOST_SHARE"]
sentinel = Path("/data/.bwbshare/another-active-job")
assert sentinel.read_text() == "keep"
# The published desktop runs as root and installs jsonpickle in root's user
# site. Test its full GUI imports in that normal configuration separately.
if os.environ.get("BWB_SMOKE_GUI"):
    from PyQt5.QtWidgets import QApplication
    app = QApplication([])
    import BwBase
    assert sentinel.read_text() == "keep"
    print("PASS: complete production GUI module imports")
    sys.exit(0)
assert not os.access("/data", os.W_OK), "test parent must really be unwritable"
os.environ["BWBHOSTSHARE"] = "/stale-desktop-hash"
client.findShareMountPoint(overwrite=True)
assert os.environ["BWBHOSTSHARE"] == os.environ["EXPECTED_HOST_SHARE"]
assert sentinel.read_text() == "keep"

image = shlex.quote(os.environ["BWB_SMOKE_IMAGE"])
def run_job(command, expected_code):
    proc = "proc.smoke-" + uuid.uuid4().hex
    output = "/tmp/" + proc + ".json"
    result = subprocess.run(["/usr/local/bin/runDockerJob.sh", output, proc,
                             "/data/.bwb", "--entrypoint /bin/sh " + image + " -c " + shlex.quote(command)],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, timeout=120)
    assert result.returncode == expected_code, (result.returncode, result.stdout)
    if expected_code == 0:
        assert json.loads(Path(output).read_text()) == [{"result": "ok"}]
    else:
        assert Path("/tmp/" + proc + "/errors/job0.42").is_dir()
        assert not Path(output).exists()
    assert sentinel.read_text() == "keep"
    assert not Path("/data/.bwbshare/" + proc).exists()

run_job("printf ok > /tmp/output/result", 0)
run_job("exit 42", 1)
print("PASS: actual module imports, workspace mapping, restart, Docker output, failure propagation")
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image")
    args = parser.parse_args()
    if os.geteuid() == 0:
        parser.error("run as a non-root Docker user so parent permissions are meaningful")
    socket_gid = os.stat("/var/run/docker.sock").st_gid
    with tempfile.TemporaryDirectory(prefix="bwb-docker-smoke-") as folder:
        base = Path(folder)
        drive = base / "drive with spaces"
        drive.mkdir()
        workspace = drive / ".bwbshare"
        workspace.mkdir()
        (workspace / "another-active-job").write_text("keep")
        (drive / ".bwb").mkdir()
        drive.chmod(0o555)
        try:
            command = ["docker", "run", "--rm", "-i", "--cap-drop", "ALL",
                       "--user", "{}:{}".format(os.getuid(), os.getgid()),
                       "--group-add", str(socket_gid),
                       "-v", str(drive) + ":/data",
                       "-v", "/var/run/docker.sock:/var/run/docker.sock",
                       "-e", "QT_QPA_PLATFORM=offscreen", "-e", "PYTHONDONTWRITEBYTECODE=1",
                       "-e", "EXPECTED_HOST_SHARE=" + str(workspace),
                       "-e", "BWB_SMOKE_IMAGE=" + args.image,
                       "--entrypoint", "python3", args.image, "-"]
            # Two separate Bwb processes, like closing and reopening the GUI.
            for attempt in (1, 2):
                print("Bwb container start {}".format(attempt), flush=True)
                subprocess.run(command, input=IN_CONTAINER, universal_newlines=True,
                               check=True, timeout=300)
            gui_command = list(command)
            gui_command[gui_command.index("--user") + 1] = "0:0"
            cap_index = gui_command.index("--cap-drop")
            del gui_command[cap_index:cap_index + 2]
            gui_command[gui_command.index("--entrypoint"):gui_command.index("--entrypoint")] = [
                "-e", "BWB_SMOKE_GUI=1"]
            subprocess.run(gui_command, input=IN_CONTAINER, universal_newlines=True,
                           check=True, timeout=120)
            assert (workspace / "another-active-job").read_text() == "keep"
        finally:
            drive.chmod(0o755)


if __name__ == "__main__":
    main()
