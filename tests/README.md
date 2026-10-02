# Shared workspace regression tests

Run the dependency-light suite as a non-root user (Python 3.6+ and `jq`):

```sh
python3 -m unittest discover -s tests -v
bash -n scripts/runDockerJob.sh
```

Permission tests intentionally require a non-root user. The suite exercises the
production methods from both desktop and VM code without loading Qt, and runs
the real shell runner with a fake Docker executable. It never uses the daemon.

After building an image, test its actual modules and nested Docker execution:

```sh
python3 tests/smoke_docker.py biodepot/bwb:latest__amd64
```

This opt-in test needs Docker socket access. It starts disposable containers,
uses only temporary test data, and checks two Bwb starts with a writable shared
directory under an unwritable mount root. It also checks output collection and
propagation of a failing child container's exit code.

## Workspace configuration

Bwb preserves existing `.bwbshare` directories. Automatic selection tries a
dedicated share mount, `/data`, other directory mounts, and finally the mapped
X11 directory. Read-only mounts, files, and private tmpfs mounts cannot supply
shared output space.

To select a pre-provisioned directory, add these arguments to the existing Bwb
container launch command (create the host directory with suitable permissions):

```sh
-v /path/to/writable/bwb-share:/bwb-share -e BWBSHARE=/bwb-share
```

Use a WSL-accessible host path when launching from WSL. Do not hard-code Docker
Desktop's generated bind-mount hashes or set `BWBHOSTSHARE`: Bwb derives that
path from the current container's mounts. Provisioning only has to make the
shared directory writable; its parent need not be writable by Bwb. An explicitly
configured but unusable `BWBSHARE` produces an error instead of silently moving
the workspace. Normal cleanup removes only the current job's directory.

The separate `/data/.bwb` log directory must also be writable when running the
GUI. This change does not grant Windows/DrvFS permissions or change ownership.

## Building the workspace hotfix

`Dockerfile.workspace-fix` applies only the three changed desktop runtime files
to the previous published multi-platform image, pinned by digest. This avoids
changing the legacy Python/Qt dependencies: the original clean Dockerfile fails
while compiling an unpinned Python 2 `peewee` dependency.

```sh
docker build --platform linux/amd64 -f Dockerfile.workspace-fix \
  --build-arg VCS_REF="$(git rev-parse HEAD)" -t biodepot/bwb:workspace-fix__amd64 .
python3 tests/smoke_docker.py biodepot/bwb:workspace-fix__amd64
```

Build and test ARM64 separately using `linux/arm64` and an `__arm64` tag. Publish
the tested architecture tags first, then update the `latest` manifest using their
immutable digests. During a staged AMD64-first release, retain the previous
ARM64 digest until the ARM64 replacement passes its checks.
