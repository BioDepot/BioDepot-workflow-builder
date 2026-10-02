#!/usr/bin/env python3
"""Check the legacy web runtime in a disposable container without host mounts."""
import argparse
import subprocess


PROBE = r'''
import sys
import peewee
import flask_login
import gevent
import geventwebsocket
from lightop import app

assert sys.version_info[:2] == (2, 7), sys.version
assert peewee.__version__ == "3.18.3", peewee.__version__
database = peewee.SqliteDatabase(":memory:")
class Probe(peewee.Model):
    value = peewee.CharField()
    class Meta:
        database = database
database.connect()
database.create_tables([Probe])
Probe.create(value="ok")
assert Probe.get().value == "ok"
database.close()
response = app.test_client().get("/")
assert response.status_code == 200, response.status_code
assert b"redirect.html" in response.data, response.data
print("PASS: Python 2 web imports, Peewee SQLite round trip, and Flask route")
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image")
    args = parser.parse_args()
    subprocess.run([
        "docker", "run", "--rm", "-i", "--network", "none", "-w", "/tmp",
        "-e", "PYTHONPATH=/web", "-e", "PYTHONDONTWRITEBYTECODE=1",
        "--entrypoint", "python", args.image, "-"],
        input=PROBE, universal_newlines=True, check=True, timeout=120)


if __name__ == "__main__":
    main()
