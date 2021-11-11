import os
import tarfile
import docker
from io import BytesIO
import shutil

client = docker.from_env()

def exec_stream(container, command):
    exitcode, outstream = container.exec_run(
        command,
        tty=True,
        stdout=True,
        stream=True
    )

    print(f"=== {command} ===")
    for o in outstream:
        print(str(o, encoding="utf8"))
    print("====" + "="*len(command) + "====")


def copy_to(src, dst, preserve=True):
    name, dst = dst.split(':')
    container = client.containers.get(name)

    srcname = os.path.basename(src) if type(src) == str and preserve else os.path.basename(dst)
    bts = BytesIO()
    bts.seek(0)

    buf = BytesIO()
    buf.seek(0)
    buf.write(src.read())
    src.seek(0)

    tar = tarfile.open(fileobj=bts, mode='w')
    try:

        if type(src) == str:
            tar.add(src, arcname=srcname)
        else:
            src.seek(0)
            info = tarfile.TarInfo(name=srcname)
            info.size = len(src.read())
            src.seek(0)
            tar.addfile(info, fileobj=src)
    finally:
        tar.close()

    bts.seek(0)
    s = str(bts.read(), encoding="utf8")
    bts.seek(0)
    print("size of s:", len(s))
    container.put_archive(os.path.dirname(dst), bts.read())
    src = buf
    bts.close()
