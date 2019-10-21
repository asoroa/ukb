#!/usr/bin/env python3

import subprocess
import os
import sys
from ukbprotocol import UkbSession

# derived variables

def start_ukb_server(port):
    '''Example function that runs UKB as a server.'''
    # Change the variables below
    ukb_path="/home/user/Ukb/Main/src" #path where ukb executables are
    ukb_kb_path="/home/user/lkbs/wn30g.bin" # path of KB
    ukb_dict_path="/home/user/lkbs/en30.lex" # path of dictionary
    # Additional arguments
    ukb_args = ['--ppr_w2w', '--dict_weight']

    ukb_server_exec = os.path.join(ukb_path, "ukb_wsd") # ukb executable
    ukb_server_exec_args = [ukb_server_exec, ["--port " + str(port) ] + ["--daemon"] + ["-K " + ukb_kb_path] + ['-D ' + ukb_dict_path] + ukb_args]
    return_code = subprocess.call(ukb_server_exec_args)
    if return_code != 0:
        print ('Could not start the UKB server : {}'.format(return_code))
        return false
    return true

def start_shutdown_server(port):
    '''Example function that shutdowns the UKB server'''
    ukb_server_exec_args = [ukb_server_exec, ["--port " + str(port) ] + ["--shutdown"]]
    subprocess.call(ukb_server_exec_args)

def test(fh, port):
    l = 1
    try:
        session = UkbSession(port)
        session.send("go")
        ctx = session.recv()
        print(ctx)
        for line in fh:
            session.send("ctx" + str(l))
            session.send(line)
            out = session.recv()
            print(out)
    except OSError as err:
        print("[E] " + str(err))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        fh = open(sys.argv[1], encoding='utf-8', errors='surrogateescape')
    else:
        fh = open(sys.stdin.fileno(), encoding='utf-8', errors='surrogateescape')
    test(fh, 6969)
