#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This is a test procedure. It generates weighted samplings of tau decay
# vertices.

import json
import os
import subprocess
import time

def run_sampling(model):
    """Run a sampling with RETRO
    """
    # Generate the configuration card
    outfile = { "1 / E": "events-E.json", "1 / E**2": "events-E2.json" }[model]
    s = 300E+03
    card = {
        "generator": {
            "theta": [90.0, 92.5],
            "energy": [model, [10**7.5, 10**10.5]],
            "position": [[-s, s], [-s, s], [0, 3E+03]]},

            "processor":  { "requested": 10000 },

        "logger": { "path": outfile },

        "topography": {
            "latitude": 43,
            "longitude": 87,
            "path": "flat/10"}}

    # Run RETRO
    with open("card.json", "wb+") as f:
        json.dump(card, f)
    p = subprocess.Popen("retro-run card.json", shell=True)
    p.communicate()
    os.remove("card.json")

if __name__ == "__main__":
    for model in ("1 / E", "1 / E**2"):
        print ""
        print "# Running `{:}` ...".format(model)
        t0 = time.time()
        run_sampling(model)
        print "  --> Done in {:.1f} s".format(time.time() - t0)
