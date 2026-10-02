"""Workspace and runner regressions; no Qt, Docker daemon or root required.

Run with: python3 -m unittest discover -s tests -v
The Docker image smoke test additionally imports the complete production modules.
"""
import ast
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
import uuid
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]


def load_methods(path, class_name, names, namespace=None):
    """Exercise production methods without importing the legacy Qt GUI stack."""
    tree = ast.parse(path.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == class_name)
    body = [n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert len(body) == len(names)
    module = ast.Module(body=body)
    if "type_ignores" in module._fields:
        module.type_ignores = []
    ns = dict(os=os, sys=sys, tempfile=tempfile)
    ns.update(namespace or {})
    exec(compile(module, str(path), "exec"), ns)
    return type(class_name, (), {name: ns[name] for name in names})


class WorkspaceCases:
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="bwb-workspace-test-")
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        env = mock.patch.dict(os.environ, {}, clear=True)
        env.start()
        self.addCleanup(env.stop)
        no_shell = mock.patch("os.system", side_effect=AssertionError("No shell provisioning"))
        no_shell.start()
        self.addCleanup(no_shell.stop)
        quiet = mock.patch("sys.stderr", new=types.SimpleNamespace(write=lambda text: None))
        quiet.start()
        self.addCleanup(quiet.stop)
        cls = load_methods(ROOT / self.source / "DockerClient.py", "DockerClient",
                           ["shareHostPath", "findShareMountPoint", "findVolumeMappings"])
        self.client = cls()
        self.client.isContainer = True
        self.client.bwb_instance_id = "test-container"
        self.mounts = []
        self.client.cli = types.SimpleNamespace(containers=lambda: [
            {"Id": "test-container", "Mounts": self.mounts}])

    def mount(self, name, rw=True, source=None):
        dest = self.base / name
        dest.mkdir(parents=True, exist_ok=True)
        self.mounts.append({"Source": source or "/daemon-only/" + name,
                            "Destination": str(dest), "RW": rw, "Type": "bind"})
        return dest

    def select(self):
        self.client.findShareMountPoint(overwrite=True)
        return Path(os.environ["BWBSHARE"])

    def test_creates_local_directory_not_daemon_path(self):
        dest = self.mount("data")
        self.assertEqual(self.select(), dest / ".bwbshare")
        self.assertEqual(os.environ["BWBHOSTSHARE"], "/daemon-only/data/.bwbshare")
        self.assertEqual(list((dest / ".bwbshare").iterdir()), [])

    @unittest.skipIf(os.geteuid() == 0, "real permission check requires non-root")
    def test_precreated_child_under_unwritable_parent(self):
        dest = self.mount("drive")
        workspace = dest / ".bwbshare"
        workspace.mkdir()
        sentinel = workspace / "another-active-job"
        sentinel.write_text("keep")
        dest.chmod(0o555)
        self.addCleanup(dest.chmod, 0o755)
        self.assertFalse(os.access(str(dest), os.W_OK))
        self.assertEqual(self.select(), workspace)
        self.assertEqual(sentinel.read_text(), "keep")
        self.assertEqual(self.select(), workspace)
        self.assertEqual(sentinel.read_text(), "keep")

    @unittest.skipIf(os.geteuid() == 0, "real permission check requires non-root")
    def test_missing_child_under_unwritable_parent_falls_back(self):
        bad = self.mount("drive")
        good = self.mount("other")
        bad.chmod(0o555)
        self.addCleanup(bad.chmod, 0o755)
        self.assertEqual(self.select(), good / ".bwbshare")

    def test_restart_refreshes_daemon_mapping_and_preserves_files(self):
        self.mount("drive")
        workspace = self.select()
        (workspace / "active-job").write_text("keep")
        self.mounts[0]["Source"] = "/new-desktop-hash"
        self.assertEqual(self.select(), workspace)
        self.assertEqual(os.environ["BWBHOSTSHARE"], "/new-desktop-hash/.bwbshare")
        self.assertEqual((workspace / "active-job").read_text(), "keep")

    def test_explicit_workspace_derives_host_path(self):
        dest = self.mount("drive")
        requested = dest / "custom workspace"
        os.environ["BWBSHARE"] = str(requested)
        os.environ["BWBHOSTSHARE"] = "/stale-host-path"
        self.assertEqual(self.select(), requested)
        self.assertEqual(os.environ["BWBHOSTSHARE"], "/daemon-only/drive/custom workspace")

    def test_readonly_mount_skipped_but_kept_for_input_translation(self):
        bad = self.mount("input", rw=False)
        good = self.mount("output")
        self.assertEqual(self.select(), good / ".bwbshare")
        self.assertEqual(self.client.bwbMounts["/daemon-only/input"], str(bad))
        self.assertFalse((bad / ".bwbshare").exists())

    def test_most_specific_readonly_mount_blocks_explicit_workspace(self):
        parent = self.mount("parent")
        self.mount("parent/nested", rw=False)
        os.environ["BWBSHARE"] = str(parent / "nested" / ".bwbshare")
        with self.assertRaisesRegex(RuntimeError, "read-only"):
            self.select()

    def test_explicit_unmapped_workspace_fails_without_creating_it(self):
        self.mount("input")
        outside = self.base / "input-other" / ".bwbshare"
        os.environ["BWBSHARE"] = str(outside)
        with self.assertRaisesRegex(RuntimeError, "not inside a Docker mount"):
            self.select()
        self.assertFalse(outside.exists())

    def test_no_usable_mount_fails_without_private_tmp_fallback(self):
        self.mount("input", rw=False)
        with self.assertRaisesRegex(RuntimeError, "No writable shared Bwb workspace"):
            self.select()
        self.assertNotIn("BWBHOSTSHARE", os.environ)

    def test_x11_and_socket_filter_uses_destination_not_hashed_source(self):
        self.mounts.extend([
            {"Source": "/desktop/hash1", "Destination": "/tmp/.X11-unix", "RW": True, "Type": "bind"},
            {"Source": "/desktop/hash2", "Destination": "/var/run/docker.sock", "RW": True, "Type": "bind"},
        ])
        self.client.findVolumeMappings()
        self.assertEqual(self.client.bwbMounts, {})
        with mock.patch("os.path.isdir", return_value=True):
            self.assertEqual(self.client.shareHostPath("/tmp/.X11-unix/.bwbshare"),
                             "/desktop/hash1/.bwbshare")

    def test_file_mount_cannot_be_workspace(self):
        dest = self.base / "file"
        dest.write_text("not a directory")
        self.mounts.append({"Source": "/daemon/file", "Destination": str(dest),
                            "RW": True, "Type": "bind"})
        with self.assertRaisesRegex(RuntimeError, "not a directory"):
            self.select()

    def test_tmpfs_without_host_source_is_not_shared(self):
        dest = self.base / "private-tmp"
        dest.mkdir()
        self.mounts.append({"Destination": str(dest), "RW": True, "Type": "tmpfs"})
        good = self.mount("data")
        self.assertEqual(self.select(), good / ".bwbshare")

    def test_nested_writable_mount_uses_its_own_source(self):
        self.mount("parent")
        dest = self.mount("parent/nested", source="/nested-host")
        os.environ["BWBSHARE"] = str(dest / ".bwbshare")
        self.select()
        self.assertEqual(os.environ["BWBHOSTSHARE"], "/nested-host/.bwbshare")

    def test_native_explicit_workspace_uses_same_host_path(self):
        self.client.isContainer = False
        requested = self.base / "native"
        os.environ["BWBSHARE"] = str(requested)
        self.assertEqual(self.select(), requested)
        self.assertEqual(os.environ["BWBHOSTSHARE"], str(requested))

    def test_relative_explicit_workspace_rejected(self):
        os.environ["BWBSHARE"] = "relative"
        with self.assertRaisesRegex(RuntimeError, "absolute"):
            self.select()

    def test_root_explicit_workspace_rejected(self):
        os.environ["BWBSHARE"] = "/"
        with self.assertRaisesRegex(RuntimeError, "filesystem root"):
            self.select()


class DesktopWorkspaceTests(WorkspaceCases, unittest.TestCase):
    source = Path("coreutils")


class VMWorkspaceTests(WorkspaceCases, unittest.TestCase):
    source = Path("VM/coreutils")


@unittest.skipUnless(shutil.which("jq"), "runner requires jq")
class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="bwb-runner-test-")
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.share = self.base / "shared space"
        self.share.mkdir()
        self.output = self.base / "output.json"
        self.proc = "proc.test-" + uuid.uuid4().hex
        self.process_path = Path("/tmp") / self.proc
        self.addCleanup(shutil.rmtree, str(self.process_path), True)
        self.calls = self.base / "docker-calls"
        self.env = dict(os.environ, BWBSHARE=str(self.share), BWBHOSTSHARE=str(self.share),
                        NWORKERS="1", BWB_TEST_CALLS=str(self.calls),
                        PATH=str(ROOT / "tests/fixtures") + os.pathsep + os.environ["PATH"])

    def run_job(self, *commands):
        return subprocess.run(["bash", str(ROOT / "scripts/runDockerJob.sh"),
                               str(self.output), self.proc, str(self.base / "logs"),
                               *(commands or ("test-image success",))], env=self.env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              universal_newlines=True, timeout=20)

    def test_success_preserves_shared_root_and_other_jobs(self):
        sentinel = self.share / "other-job"
        sentinel.write_text("keep")
        result = self.run_job()
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual(json.loads(self.output.read_text()), [{"result": "ok"}])
        self.assertEqual(sentinel.read_text(), "keep")
        self.assertFalse((self.share / self.proc).exists())
        self.assertFalse((self.process_path / "errors").exists())

    def test_docker_failure_is_nonzero_and_has_real_error_code(self):
        result = self.run_job("test-image fail")
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue((self.process_path / "errors/job0.42").is_dir())
        self.assertFalse(self.output.exists())

    def test_multiple_workers_propagate_failure(self):
        self.env["NWORKERS"] = "3"
        result = self.run_job("test-image success", "test-image fail", "test-image success")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(self.calls.read_text().splitlines()), 3)
        self.assertFalse(self.output.exists())

    def test_multiple_workers_success(self):
        self.env["NWORKERS"] = "3"
        result = self.run_job(*(["test-image success"] * 6))
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual(len(json.loads(self.output.read_text())), 6)
        self.assertEqual(len(self.calls.read_text().splitlines()), 6)

    def test_missing_environment_never_launches_docker(self):
        del self.env["BWBSHARE"]
        result = self.run_job()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())

    def test_existing_job_directory_is_not_deleted(self):
        existing = self.share / self.proc
        existing.mkdir()
        (existing / "keep").write_text("keep")
        result = self.run_job()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())
        self.assertEqual((existing / "keep").read_text(), "keep")

    def test_file_instead_of_workspace_never_launches_docker(self):
        bad = self.base / "not-a-directory"
        bad.write_text("keep")
        self.env["BWBSHARE"] = str(bad)
        result = self.run_job()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())
        self.assertEqual(bad.read_text(), "keep")
        log = self.base / "logs" / self.proc / "logs/log0"
        self.assertIn("Cannot create shared workspace", log.read_text())

    @unittest.skipIf(os.geteuid() == 0, "real permission check requires non-root")
    def test_unwritable_workspace_never_launches_docker(self):
        self.share.chmod(0o555)
        self.addCleanup(self.share.chmod, 0o755)
        result = self.run_job()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.calls.exists())

    @unittest.skipIf(os.geteuid() == 0, "real permission check requires non-root")
    def test_precreated_workspace_under_unwritable_parent_runs(self):
        parent = self.base / "drive"
        parent.mkdir()
        workspace = parent / ".bwbshare"
        workspace.mkdir()
        parent.chmod(0o555)
        self.addCleanup(parent.chmod, 0o755)
        self.env.update(BWBSHARE=str(workspace), BWBHOSTSHARE=str(workspace))
        result = self.run_job()
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertTrue(workspace.is_dir())

    def test_bad_worker_count_never_launches_docker(self):
        self.env["NWORKERS"] = "0"
        self.assertNotEqual(self.run_job().returncode, 0)
        self.assertFalse(self.calls.exists())

    def test_traversal_identifier_rejected(self):
        self.proc = "../unsafe"
        self.assertNotEqual(self.run_job().returncode, 0)
        self.assertFalse(self.calls.exists())

    def test_metadata_with_spaces_and_json_types(self):
        result = self.run_job("test-image metadata")
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual(json.loads(self.output.read_text()),
                         [{"key with spaces": {"answer": 42}, "text": "hello world"}])


class FinishTests(unittest.TestCase):
    def test_failure_never_emits_outputs_or_success(self):
        for source in ("coreutils", "VM/coreutils"):
            cls = load_methods(ROOT / source / "BwBase.py", "OWBwBWidget", ["onRunFinished"],
                               {"Qt": types.SimpleNamespace(red="red", green="green")})
            for code, status in ((1, 0), (0, -1), (1, -1)):
                with self.subTest(source=source, code=code, status=status):
                    widget = cls()
                    widget.pConsole = mock.Mock()
                    widget.bgui = mock.Mock()
                    widget.repeat = True
                    widget.status = "running"
                    for name in ("reenableExec", "setStatusMessage", "initTriggers", "updateOutputs"):
                        setattr(widget, name, mock.Mock())
                    widget.onRunFinished(code, status)
                    self.assertEqual(widget.status, "error")
                    widget.updateOutputs.assert_not_called()
                    self.assertNotIn("Finished", [c[0][0] for c in widget.pConsole.writeMessage.call_args_list])

    def test_success_still_emits_outputs(self):
        for source in ("coreutils", "VM/coreutils"):
            cls = load_methods(ROOT / source / "BwBase.py", "OWBwBWidget", ["onRunFinished"],
                               {"Qt": types.SimpleNamespace(red="red", green="green")})
            widget = cls()
            widget.pConsole = mock.Mock()
            widget.bgui = mock.Mock()
            widget.repeat = False
            widget.status = "running"
            for name in ("reenableExec", "setStatusMessage", "initTriggers", "updateOutputs"):
                setattr(widget, name, mock.Mock())
            widget.onRunFinished(0, 0)
            self.assertEqual(widget.status, "finished")
            widget.updateOutputs.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
