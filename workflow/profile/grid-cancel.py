#!/usr/bin/env python3
import sys
import htcondor


def print_and_exit(s):
    print(s)
    exit()

jobID, UUID, clusterID = sys.argv[1].split("_")

schedd = htcondor.Schedd()
schedd.act(htcondor.JobAction.Remove, f'ClusterId == "{clusterID}"')
