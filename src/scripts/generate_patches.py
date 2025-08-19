import argparse
import pathlib
import subprocess
import io
import re

RX = re.compile(r"diff --git a/([a-zA-Z\./]+) b/([a-zA-Z\./]+)")


parser = argparse.ArgumentParser()
parser.add_argument("--input", default=pathlib.Path("vendor", "jess"), type=pathlib.Path)
parser.add_argument("--base", default="testing")
parser.add_argument("--head", default="martin-update")
parser.add_argument("--output", default=pathlib.Path("patches"), type=pathlib.Path)
args = parser.parse_args()



proc = subprocess.run(["git", "checkout", args.head], cwd=args.input)
proc.check_returncode()

proc = subprocess.run(["git", "diff", args.base], cwd=args.input, capture_output=True)
proc.check_returncode()
patch = proc.stdout.decode()

proc = subprocess.run(["git", "checkout", args.base], cwd=args.input)
proc.check_returncode()



filename = None
buffer = []

for line in io.StringIO(patch):
    m = RX.match(line)
    if m:
        if filename is not None:
            with args.output.joinpath(f"{filename}.patch").open("w") as dst:
                dst.writelines(buffer)
                buffer.clear()
        path = pathlib.Path(m.group(1))
        filename = path.name
    buffer.append(line)

if filename is not None:
    with args.output.joinpath(f"{filename}.patch").open("w") as dst:
        dst.writelines(buffer)
        buffer.clear()

