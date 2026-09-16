#!/usr/bin/env python3
"""Dump the PETSc API description produced by PETSc's own getAPI.py as JSON.

Usage:  python3 getapi_dump.py PETSC_DIR OUTPUT.json[.gz]

Works with both layouts of getAPI.py:
  - PETSc <= 3.24: config/utils/getAPI.py, getAPI() run with cwd = PETSC_DIR, 9-tuple
  - PETSc >= 3.25: lib/petsc/bin/getAPI.py, getAPI(directory), 10-tuple (adds functiontypedefs)

Only the Python standard library is needed. The JSON is deterministic (sorted keys,
sets sorted) so two dumps of the same source tree are byte-identical and dumps of two
releases can be diffed.
"""
import gzip
import json
import os
import re
import sys


def petsc_version(petsc_dir):
    txt = open(os.path.join(petsc_dir, "include", "petscversion.h")).read()
    v = [re.search(r"#define\s+PETSC_VERSION_%s\s+(\d+)" % k, txt).group(1)
         for k in ("MAJOR", "MINOR", "SUBMINOR")]
    return ".".join(v)


def load_getapi(petsc_dir):
    for rel in ("lib/petsc/bin", "config/utils"):
        path = os.path.join(petsc_dir, rel, "getAPI.py")
        if os.path.isfile(path):
            sys.path.insert(0, os.path.dirname(path))
            import getAPI  # noqa: E402
            return getAPI, rel
    sys.exit("getAPI.py not found under %s (needs PETSc >= 3.23)" % petsc_dir)


def make_walk_deterministic():
    """getAPI.py records the first definition it meets when a function is defined in several files
    (e.g. device.c vs device.cxx) and walks the tree with os.walk/os.listdir, whose order depends on
    the filesystem. Sort them so two machines produce the same snapshot."""
    _walk, _listdir = os.walk, os.listdir

    def walk(top, *a, **k):
        for root, dirs, files in _walk(top, *a, **k):
            # MATLAB mex stubs re-declare PETSc functions with other signatures
            dirs[:] = sorted(d for d in dirs if d != "mex-scripts")
            files.sort()
            yield root, dirs, files

    os.walk = walk
    os.listdir = lambda *a, **k: sorted(_listdir(*a, **k))


def run_getapi(getAPI, layout, petsc_dir):
    make_walk_deterministic()
    # both layouts open include files relative to the cwd, so run inside PETSC_DIR
    cwd = os.getcwd()
    os.chdir(petsc_dir)
    try:
        if layout == "lib/petsc/bin":
            result = getAPI.getAPI(petsc_dir)
        else:
            result = getAPI.getAPI()
    finally:
        os.chdir(cwd)
    names9 = ["classes", "enums", "senums", "typedefs", "structs", "funcs",
              "includefiles", "mansecs", "submansecs"]
    names10 = ["classes", "enums", "senums", "typedefs", "functiontypedefs", "structs",
               "funcs", "includefiles", "mansecs", "submansecs"]
    names = names10 if len(result) == 10 else names9
    if len(result) != len(names):
        sys.exit("unexpected getAPI() return arity %d" % len(result))
    return dict(zip(names, result))


def to_json(obj):
    """Recursively convert getAPI's plain objects into JSON-able data."""
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    if isinstance(obj, dict):
        return {str(k): to_json(obj[k]) for k in sorted(obj, key=str)}
    if isinstance(obj, (set, frozenset)):
        return sorted(to_json(x) for x in obj)
    if isinstance(obj, (list, tuple)):
        return [to_json(x) for x in obj]
    if hasattr(obj, "__dict__"):
        d = {k: to_json(v) for k, v in vars(obj).items() if not k.startswith("_")}
        d["_class"] = type(obj).__name__
        return d
    return str(obj)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    petsc_dir = os.path.abspath(sys.argv[1])
    out = sys.argv[2]
    getAPI, layout = load_getapi(petsc_dir)
    api = run_getapi(getAPI, layout, petsc_dir)
    data = {
        "petsc_version": petsc_version(petsc_dir),
        "getapi_layout": layout,
    }
    data.update(to_json(api))
    text = json.dumps(data, indent=1, sort_keys=True) + "\n"
    if out.endswith(".gz"):
        # mtime=0 so the archive is reproducible byte for byte
        with open(out, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as f:
            f.write(text.encode())    # no file name or mtime in the header
    else:
        with open(out, "w") as f:
            f.write(text)
    nfun = len(data["funcs"]) + sum(len(c["functions"]) for c in data["classes"].values())
    print("PETSc %s (%s): %d classes, %d functions, %d enums, %d structs -> %s" % (
        data["petsc_version"], layout, len(data["classes"]), nfun,
        len(data["enums"]), len(data["structs"]), out))


if __name__ == "__main__":
    main()
