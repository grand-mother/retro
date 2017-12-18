# -*- coding: utf-8 -*-
#
#  Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
#  Author: Valentin NIESS (niess@in2p3.fr)
#
#  Radio nEuTRino simulatiOn (RETRO).
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>

import json
import math
import os
import subprocess
import tempfile
# Custom imports.
import danton
from . import TAU_CTAU, TAU_MASS


class DantonError(Exception):
    """Custom exception for a DANTON error.
    """

    def __init__(self, message):
        super(DantonError, self).__init__(message)


def PrimarySampler(primary, generator, topography, topo_handle):
    """Closure for sampling the primary flux with DANTON.
    """
    infile, outfile = "", ""

    class ManageTemp:
        """Context for ensuring that temporary files are deleted at exit.
        """

        def __enter__(self): pass

        def __exit__(self, type, value, traceback):
            if os.path.exists(infile):
                os.remove(infile)
            if os.path.exists(outfile):
                os.remove(outfile)

    def get_tempfile():
        """Utility function for getting a reusable temporary file.
        """
        f = tempfile.NamedTemporaryFile(
            prefix="retro.", suffix=".json", dir=".", delete=False)
        name = f.name
        f.close()
        return name

    # Compute the energy range.
    emin, emax = float("inf"), 0.
    for _, opts in generator:
        e0, e1 = opts["energy"]
        if isinstance(e0, basestring):
            e0, e1 = e1
        if e0 < emin:
            emin = e0
        if e1 > emax:
            emax = e1
    emax *= 1E+03
    weight_factor = TAU_MASS / TAU_CTAU * (1. / emin - 1. / emax)

    # Configure for running DANTON.
    max_events = primary["events"]
    try:
        requested = primary["requested"]
    except KeyError:
        requested = -1
    with ManageTemp():
        infile, outfile = get_tempfile(), get_tempfile()
        particle = {"tau": None, "tau~": None}
        sampler = {
            "altitude": None,
            "elevation": None,
            "energy": None,
            "weight": particle}
        flux_model = ["power-law", {
            "energy": [emin, emax],
            "exponent": -2.,
            "weight": 1.}]
        card = {
            "events": max_events,
            "requested": requested,
            "output-file": outfile,
            "mode": "backward",
            "longitudinal": primary["longitudinal"],
            "decay": False,
            "particle-sampler": sampler,
            "primary-flux": {
                "nu_tau": flux_model,
                "nu_tau~": flux_model}}
        flat = topography["path"].startswith("flat")
        if flat:
            card["earth-model"] = {"sea": False}
        else:
            raise ValueError("non flat topography is not yet supported")
        run_cmd = "danton {:}".format(infile)

    def sample(pid, position, energy, direction):
        """Sample a bunch of primaries using DANTON.
        """
        with ManageTemp():
            # Configure the sampler.
            sampler["altitude"] = topo_handle.local_to_lla(position)[2]
            theta, _ = topo_handle.local_to_angular(position, direction)
            sampler["elevation"] = theta - 90.
            sampler["energy"] = energy
            if pid > 0.:
                particle["tau"], particle["tau~"] = 1., 0.
            else:
                particle["tau"], particle["tau~"] = 0., 1.

            # Dump the steering card.
            with open(infile, "wb+") as f:
                json.dump(card, f)

            # Run DANTON.
            p = subprocess.Popen(run_cmd, shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            _, stderr = p.communicate()
            if stderr:
                raise DantonError(stderr)

            # Compute the decay weight.
            p = math.sqrt((energy - TAU_MASS) * (energy + TAU_MASS))
            wd = weight_factor / p

            # Parse the result.
            if not os.path.exists(outfile):
                return [], 0
            n_events = 0
            primaries = []
            for event in danton.iter_event(outfile):
                n_events = event.id
                d = event.decay[0]
                primaries.append([event.weight * wd, event.primary.energy,
                                  d.generation, d.tau_i.position])
            os.remove(outfile)
            if len(primaries) < requested:
                n_events = max_events
            return primaries, n_events
    return sample
