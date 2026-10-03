#!/usr/bin/env python3
"""Test startup failures and the complete GUI's widget registry in a built image.

Usage: python3 tests/smoke_startup.py IMAGE
Uses only disposable containers and temporary data; does not publish ports.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

GUI = r'''
import os
import subprocess
import sys
sys.path.insert(0, "/coreutils")
if os.environ.get("BWB_TEST_HIDE_MOUNTINFO"):
    import builtins
    import io
    import DockerClient
    DockerClient.open = lambda path, *args, **kwargs: (
        io.StringIO("unrecognized mount layout") if path == "/proc/self/mountinfo"
        else builtins.open(path, *args, **kwargs))
from Orange.canvas import __main__ as canvas_main
from Orange.canvas.registry import global_registry

def checked_event_loop(app):
    widgets = global_registry().widgets()
    assert len(widgets) == 34, [w.name for w in widgets]
    assert os.environ["BWBSHARE"] == "/data/.bwbshare"
    host_share = os.environ["BWBHOSTSHARE"]
    assert host_share != "/stale-daemon-path", host_share
    # Docker Desktop may report a daemon path different from the client path.
    visible = subprocess.check_output([
        "docker", "run", "--rm", "--network", "none",
        "-v", host_share + ":/shared:ro", "--entrypoint", "cat",
        os.environ["BWB_SMOKE_IMAGE"], "/shared/another-active-job"],
        universal_newlines=True, timeout=60)
    assert visible == "keep", visible
    # Discovery alone is insufficient: also construct a standard Bwb widget.
    from importlib import import_module
    from unittest import mock
    from BwbStartup import get_docker_client
    description = next(w for w in widgets if w.name == "bash_utils")
    module, name = description.qualified_name.rsplit(".", 1)
    widget_class = getattr(import_module(module), name)
    with mock.patch.object(widget_class, "startJob", side_effect=AssertionError("Unexpected job launch")):
        widget = widget_class()
        assert widget.dockerClient is get_docker_client()
        assert not widget.jobRunning
        widget.deleteLater()
    print("PASS: full GUI startup discovered all 34 standard widgets", flush=True)
    return 0

canvas_main.CanvasApplication.exec_ = checked_event_loop
canvas_main.show_survey = lambda: None
canvas_main.check_for_updates = lambda: None
status = canvas_main.main(["orange-canvas", "--force-discovery", "--no-welcome", "--no-splash"])
print("GUI_EXIT=" + str(status), flush=True)
sys.exit(status)
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="bwb-startup-smoke-") as folder:
        data = Path(folder) / "data with spaces"
        data.mkdir()
        logs = data / ".bwb"
        share = data / ".bwbshare"
        logs.mkdir()
        share.mkdir()
        sentinel = share / "another-active-job"
        sentinel.write_text("keep")
        log_sentinel = logs / "another-active-job"
        log_sentinel.write_text("keep logs")
        x11 = Path(folder) / "x11"
        x11.mkdir()
        cases = ("writable", "permission_denied", "read_only", "logs_denied",
                 "share_denied", "missing_data", "missing_socket", "readonly_x11")
        for case in cases:
            for path in (data, logs, share):
                path.chmod(0o777)
            if case == "permission_denied":
                data.chmod(0o555)
            elif case == "logs_denied":
                logs.chmod(0o555)
            elif case == "share_denied":
                share.chmod(0o555)
            try:
                command = ["docker", "run", "--rm", "-i",
                           "--cap-drop", "DAC_OVERRIDE", "--cap-drop", "DAC_READ_SEARCH",
                           "-e", "QT_QPA_PLATFORM=offscreen",
                           "-e", "PYTHONDONTWRITEBYTECODE=1",
                           "-e", "BWBHOSTSHARE=/stale-daemon-path",
                           "-e", "BWB_SMOKE_IMAGE=" + args.image]
                if case != "missing_data":
                    command += ["-v", str(data) + ":/data" + (":ro" if case == "read_only" else "")]
                if case != "missing_socket":
                    command += ["-v", "/var/run/docker.sock:/var/run/docker.sock"]
                if case == "readonly_x11":
                    command += ["-v", str(x11) + ":/tmp/.X11-unix:ro"]
                # Successful GUI starts are repeated to exercise cached settings
                # and ensure that preflight never destroys another job's files.
                modes = ("preflight", "gui", "gui", "hostname_fallback") if case == "writable" else ("entrypoint", "gui")
                if case == "readonly_x11":
                    modes = ("preflight", "gui")
                for mode in modes:
                    mode_command = command + (["-e", "BWB_TEST_HIDE_MOUNTINFO=1"]
                                              if mode == "hostname_fallback" else [])
                    if mode == "entrypoint":
                        invocation = mode_command + [args.image]
                    elif mode == "preflight":
                        invocation = mode_command + ["--entrypoint", "python3", args.image,
                                                "/coreutils/BwbStartup.py"]
                    else:
                        invocation = mode_command + ["--entrypoint", "python3", args.image, "-"]
                    result = subprocess.run(
                        invocation, input=GUI if mode in ("gui", "hostname_fallback") else "",
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        universal_newlines=True, timeout=90)
                    expected = 0 if case in ("writable", "readonly_x11") else 1
                    assert result.returncode == expected, (case, mode, result.returncode, result.stdout)
                    if expected:
                        assert "Bwb cannot start:" in result.stdout, result.stdout
                        assert "-v /path/to/writable-directory:/data" in result.stdout, result.stdout
                        assert "PASS: full GUI" not in result.stdout, result.stdout
                        assert "The widget will not be shown" not in result.stdout, result.stdout
                    elif mode in ("gui", "hostname_fallback"):
                        assert "PASS: full GUI startup discovered all 34" in result.stdout, result.stdout
                    assert sentinel.read_text() == "keep"
                    assert log_sentinel.read_text() == "keep logs"
                    assert not list(data.rglob(".bwb-probe-*"))
                    print("PASS: {} / {} / exit {}".format(case, mode, result.returncode), flush=True)
            finally:
                for path in (data, logs, share):
                    if path.exists():
                        path.chmod(0o777)


if __name__ == "__main__":
    main()
