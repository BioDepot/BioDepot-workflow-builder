"""Explicit startup preflight, separate from Orange's widget imports."""
import sys

_client = None


def get_docker_client():
    global _client
    if _client is None:
        from DockerClient import DockerClient
        _client = DockerClient("unix:///var/run/docker.sock", "local")
    return _client


def startup_preflight():
    """Return a checked client, or an actionable fatal startup error."""
    try:
        client = get_docker_client()
        client.checkStartupWorkspace()
        return client
    except Exception as exc:
        raise RuntimeError(
            "Bwb cannot start: {}\n\n"
            "Bwb requires a writable host directory mounted at /data. "
            "Restart Bwb with a different writable host directory mounted at /data, "
            "for example:\n"
            "  -v /path/to/writable-directory:/data\n"
            "The Bwb user must be able to create files there and in /data/.bwb "
            "and /data/.bwbshare (or the explicit BWBSHARE directory).\n"
            "If mount discovery failed, also check the Docker socket mapping "
            "and container identification. No workflow has been started.".format(exc))


def main():
    try:
        startup_preflight()
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
