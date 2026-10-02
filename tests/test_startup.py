"""Tests of the startup contract without loading the legacy GUI dependencies."""
import ast
import contextlib
import importlib.util
import io
from pathlib import Path
import subprocess
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


class StartupContractTests(unittest.TestCase):
    def test_builds_share_a_digest_pinned_modern_docker_client(self):
        pins = []
        for filename in ("Dockerfile", "Dockerfile.workspace-fix"):
            source = (ROOT / filename).read_text()
            pins.append(next(line for line in source.splitlines()
                             if line.startswith("ARG DOCKER_CLI_IMAGE=")))
            self.assertIn("docker:29.8.2-cli@sha256:", pins[-1])
            self.assertIn("COPY --from=docker_cli /usr/local/bin/docker /usr/bin/docker", source)
            self.assertIn("COPY --from=docker_cli /usr/local/libexec/docker/cli-plugins/docker-buildx", source)
        self.assertEqual(pins[0], pins[1])

    def load_startup(self, directory):
        spec = importlib.util.spec_from_file_location(
            "tested_startup", str(ROOT / directory / "BwbStartup.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_preflight_success_returns_the_checked_shared_client(self):
        for directory in ("coreutils", "VM/coreutils"):
            module = self.load_startup(directory)
            client = mock.Mock()
            module._client = client
            self.assertIs(module.startup_preflight(), client)
            client.checkStartupWorkspace.assert_called_once_with()

    def test_permission_or_discovery_failure_is_actionable_and_exits_nonzero(self):
        for directory in ("coreutils", "VM/coreutils"):
            for reason in ("Permission denied", "shared mount discovery failed"):
                module = self.load_startup(directory)
                client = mock.Mock()
                client.checkStartupWorkspace.side_effect = RuntimeError(reason)
                module._client = client
                stderr = io.StringIO()
                with contextlib.redirect_stderr(stderr):
                    self.assertEqual(module.main(), 1)
                message = stderr.getvalue()
                self.assertIn(reason, message)
                self.assertIn("Bwb cannot start", message)
                self.assertIn("-v /path/to/writable-directory:/data", message)
                self.assertIn("No workflow has been started", message)

    def test_widget_definition_does_not_construct_a_docker_client(self):
        for directory in ("coreutils", "VM/coreutils"):
            tree = ast.parse((ROOT / directory / "BwBase.py").read_text())
            widget = next(n for n in tree.body
                          if isinstance(n, ast.ClassDef) and n.name == "OWBwBWidget")
            assignments = [n for n in widget.body if isinstance(n, ast.Assign)]
            for assignment in assignments:
                self.assertFalse(any(
                    isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                    and n.func.id in ("DockerClient", "get_docker_client")
                    for n in ast.walk(assignment)))

    def test_preflight_is_before_canvas_and_widget_discovery(self):
        for filename in ("orangePatches/__main__.py",
                         "VM/orange3/Orange/canvas/__main__.py"):
            source = (ROOT / filename).read_text()
            preflight = source.index("        startup_preflight()")
            self.assertLess(preflight, source.index("    canvas_window = CanvasMainWindow()"))
            self.assertLess(preflight, source.index("        widget_discovery.run("))
            self.assertIn('QMessageBox.critical(None, "Bwb cannot start", message)', source)

    def test_launcher_preserves_failed_program_exit_code(self):
        for filename in ("scripts/startBwb.sh", "scripts/startSingleBwb.sh", "VM/startBwb.sh"):
            result = subprocess.run(
                ["bash", str(ROOT / filename), "/bin/false"],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                universal_newlines=True, timeout=10)
            self.assertEqual(result.returncode, 1, (filename, result.stdout))


if __name__ == "__main__":
    unittest.main()
