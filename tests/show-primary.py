#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This is a test procedure. It plots the density of decaying taus as
# function of the tau energy. Use the `run-primary.py` script in order to
# generate the data.

import json
import numpy
import matplotlib.pyplot as plt

with open("events-primary.json", "rb") as f:
    data = json.load(f)
    energy, density = map(numpy.array, (data["energy"], data["density"]))

plt.style.use("plugins/ccin2p3/deps/mplstyle-l3/style/l3.mplstyle")
plt.figure()
plt.loglog(energy, density[:, 0], "k-")
plt.errorbar(energy, density[:, 0], yerr=density[:, 1],
             fmt="k.")
plt.xlabel(r"$\tau$ energy (GeV)")
plt.ylabel(r"density (GeV$^{-1}$ m$^{-3}$ sr$^{-1}$ s$^{-1}$)")
plt.show()
