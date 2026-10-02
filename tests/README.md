# Shared workspace regression tests

Run the dependency-light suite as a non-root user (Python 3.6+ and `jq`):

```sh
python3 -m unittest discover -s tests -v
bash -n scripts/runDockerJob.sh
```

Permission tests intentionally require a non-root user. The suite exercises the
production methods from both desktop and VM code without loading Qt, and runs
the real shell runner with a fake Docker executable. It never uses the daemon.

After building an image, test startup, full widget discovery, and nested Docker execution:

```sh
python3 tests/smoke_docker.py biodepot/bwb:latest__amd64
python3 tests/smoke_startup.py biodepot/bwb:latest__amd64
python3 tests/smoke_desktop.py biodepot/bwb:latest__amd64
```

This opt-in test needs Docker socket access. It starts disposable containers,
uses only temporary test data, and checks the runner's ability to use an existing
writable child under an unwritable parent, plus nested-container output and
failure propagation. This low-level runner scenario is not a valid GUI startup:
the GUI now requires a writable `/data` root as described below.

The startup test runs the actual GUI entry point twice, requires all 34 default
widgets to be discovered, and instantiates a standard widget while ensuring
other jobs' logs and scratch files survive. It then checks both container and GUI startup failures
with unwritable/read-only `/data`, unwritable logs or scratch space, missing
`/data` mappings, and a missing Docker socket. Failures must exit 1 with an
actionable message, before widget discovery. Tests never publish ports or run a
user workflow.

The desktop test additionally starts the image's unmodified default command,
including Xvfb, Fluxbox, and the real Bwb application. It uses a fresh registry
cache, checks for 34 widgets and a responsive internal web endpoint, and stops
only its own disposable container. No host ports or host X11 sockets are used.
For ARM64 under QEMU, use `--startup-timeout 300 --probe-timeout 120`; importing
the registry in the probe is also significantly slower under emulation.

## Workspace configuration

Bwb checks its workspace immediately at startup, before loading widgets. The
normal container launch requires a writable host directory mounted at `/data`.
The check performs real temporary directory/file writes in `/data`, `/data/.bwb`
(logs), and `/data/.bwbshare` (scratch). Existing job files are preserved. A
failure aborts startup with a nonzero exit and asks the user to restart with a
writable host directory mounted at `/data`, for example:

```sh
-v /path/to/writable-directory:/data
```

There is no silent startup fallback to X11 or an unrelated input mount. Read-only
mounts, files, and private tmpfs mounts cannot supply shared output space. Mount
discovery errors are reported separately from filesystem permission errors;
container identification accepts a verified short-ID hostname when mountinfo
does not expose the full ID.

To select a pre-provisioned directory, add these arguments to the existing Bwb
container launch command (create the host directory with suitable permissions):

```sh
-v /path/to/writable/bwb-share:/bwb-share -e BWBSHARE=/bwb-share
```

Use a WSL-accessible host path when launching from WSL. Do not hard-code Docker
Desktop's generated bind-mount hashes or set `BWBHOSTSHARE`: Bwb derives that
path from the current container's mounts. An explicit `BWBSHARE` overrides the
scratch location, but does not remove the requirement for writable `/data` and
`/data/.bwb`. Its own parent need not be writable if the scratch directory already
exists. An unusable explicit workspace aborts startup instead of silently moving
the workspace. Normal cleanup removes only the current job's directory.

The separate `/data/.bwb` log directory must also be writable when running the
GUI. This change does not grant Windows/DrvFS permissions or change ownership.

## Building the workspace hotfix

`Dockerfile.workspace-fix` applies the changed runtime modules, startup entry
points, and runner to the previous published multi-platform image, pinned by digest. This avoids
changing the legacy Python/Qt dependencies during a runtime hotfix.

The clean Dockerfile's Python 2 web layer pins `peewee==3.18.3`: the unbounded
requirement selected Python-3-only Peewee 4 and failed during compilation.
Both GHCR workflows build and run `python3 tests/smoke_web.py IMAGE` before
pushing. That test checks web imports, an in-memory SQLite round trip, and the
Flask landing page without mounting host directories or publishing ports.

Both Dockerfiles copy Docker CLI 29.8.2 and its Buildx plugin from the same
digest-pinned, multi-architecture official Docker image. This replaces the old
24.0.2/API 1.43 client without upgrading or reconfiguring the host daemon. The
CLI negotiates an API supported by the server; the integration test checks the
CLI connection separately from the Python Docker SDK and runs real child jobs.
In particular, Docker 29.0-29.2 normally require API 1.44+, while 29.3+ lowered
the minimum to 1.40. Inspect the actual server minimum rather than assuming all
29.x hosts behave identically.

```sh
docker build --platform linux/amd64 -f Dockerfile.workspace-fix \
  --build-arg VCS_REF="$(git rev-parse HEAD)" -t biodepot/bwb:workspace-fix__amd64 .
python3 tests/smoke_docker.py biodepot/bwb:workspace-fix__amd64
python3 tests/smoke_startup.py biodepot/bwb:workspace-fix__amd64
```

Build and test ARM64 separately using `linux/arm64` and an `__arm64` tag. Publish
the tested architecture tags first, then update the `latest` manifest using their
immutable digests. During a staged AMD64-first release, retain the previous
ARM64 digest until the ARM64 replacement passes its checks.
