import docker
from utils import *
import os
from base64 import b64encode
import json

client = docker.from_env()

def containers_list(client):
    print("\n=== CONTAINER LIST ===\n")
    for c in client.containers.list(all=True):
        print(c.name, c.id)
    print("======================")

#def generate(outstream, somebytes, debug):
#    s = ''
#    for o in outstream:
#        t = str(o, encoding="utf8").strip() + '\n'
#        s += t
#        if debug:
#            print("=== DEBUG OUTSTREAM GENERATE ===", t)
#        yield t
#    ret_value = {
#        'output.tar': b64encode(filebytes.read()).decode('utf-8'),
#        'text' : s
#    }
#    yield json.dumps(ret_value, sort_keys=False, indent=4, ensure_ascii=False)

def create_hobotnica(cntmtrx, anno, debug=False):
    hbt_container = client.containers.create(
        image="diffexprimage",
        tty=True
    )
    if debug:
        containers_list(client)
        print()
        print("=== HBT_CONTAINER NAME | CREATE ===")
        print(hbt_container.name)
        print("===================================")
        print()


    cntmtrx.seek(0)
    anno.seek(0)

    buf = BytesIO()
    buf.seek(0)
    buf.write(cntmtrx.read())
    buf.seek(0)
    copy_to(buf, str(hbt_container.name) + ':data/countmatrix.txt', preserve=False)

    buf = BytesIO()
    buf.seek(0)
    buf.write(anno.read())
    buf.seek(0)
    copy_to(buf, str(hbt_container.name) + ':data/annotation.txt', preserve=False)

    cntmtrx.seek(0)
    anno.seek(0)

    if debug:
        print("\nSTARTING CONTAINER...\n")

    hbt_container.start()
    return hbt_container



def call_hobotnica(hbt_container, debug=False):
    if debug:
        containers_list(client)
        print()
        print("=== HBT_CONTAINER NAME | CALL ===")
        print(hbt_container.name)
        hbt_container.start()
        print(hbt_container.status)
        print("=================================")
        print()
    

    if debug:
        exec_stream(hbt_container, "ls -lah")
        exec_stream(hbt_container, "ls -lah data")

    
    outstream = None
    exitcode, outstream = hbt_container.exec_run(
            'Rscript run.R data/countmatrix.txt data/annotation.txt output',
            tty=True,
            stdout=True,
            stream=True
    )


    if debug:
        print("\n=== Script is probably running... ===\n")
        containers_list(client)
        exec_stream(hbt_container, "ls -lah")
    s = ''
    for o in outstream:
        t = str(o, encoding="utf8").strip() + '\n'
        s += t
        if debug:
            print("=== DEBUG OUTSTREAM MAIN ===", t)
        yield t

    f = BytesIO()
    f.seek(0)
    if debug:
        print("\n=== PROBABLY COMPUTED ===\n")
        exec_stream(hbt_container, "ls -lah")
        containers_list(client)
    bits, stat = hbt_container.get_archive('output')

    print("\n=== STAT ===\n")
    print(stat)
    print("============")
    for chunk in bits:
        f.write(chunk)
    f.seek(0)
    ret_value = {
        'output.tar': b64encode(f.read()).decode('utf-8'),
        'text' : s
    }
    yield "=========="
    yield json.dumps(ret_value, sort_keys=False, indent=4, ensure_ascii=False)

    hbt_container.stop()
    hbt_container.remove()




