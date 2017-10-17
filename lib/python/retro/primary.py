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
import os
import subprocess
import tempfile
# Custom imports.
import danton

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
            if os.path.exists(infile): os.remove(infile)
            if os.path.exists(outfile): os.remove(outfile)

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
        if e0 < emin: emin = e0
        if e1 > emax: emax = e1
    emax *= 1E+03

    # Configure for running DANTON.
    with ManageTemp():
        infile, outfile = get_tempfile(), get_tempfile()
        particle = { "tau" : None, "tau~" : None }
        sampler = {
            "altitude" : None,
            "elevation" : None,
            "energy" : None,
            "weight" : particle }
        flux_model = [ "power-law", {
            "energy" : [ emin, emax ],
            "exponent" : -2.,
            "weight" : 1. }]
        card = {
            "events" : primary["events"],
            "output-file" : outfile,
            "mode" : "backward",
            "longitudinal" : primary["longitudinal"],
            "decay" : False,
            "particle-sampler" : sampler,
            "primary-flux" : {
                "nu_tau" : flux_model,
                "nu_tau~" : flux_model }}
        flat = topography["path"].startswith("flat")
        if flat:
            card["earth-model"] = { "sea" : False }
        else:
            raise ValueError("non flat topography is not yet supported")
        run_cmd = "danton {:}".format(infile)

    def sample(pid, position, energy, direction):
        """Sample a bunch of primaries using DANTON.
        """
        with ManageTemp():
            # Configure the sampler.
            sampler["altitude"] = topo_handle.local_to_lla(position)[2]
            _, elevation = topo_handle.local_to_horizontal(position, direction)
            sampler["elevation"] = elevation
            sampler["energy"] = energy
            if pid > 0.: particle["tau"], particle["tau~"] = 1., 0.
            else: particle["tau"], particle["tau~"] = 0., 1.

            # Dump the steering card.
            with open(infile, "wb+") as f: json.dump(card, f)

            # Run DANTON.
            p = subprocess.Popen(run_cmd, shell=True, stdout=subprocess.PIPE,
              stderr=subprocess.PIPE)
            _, stderr = p.communicate()
            if stderr: raise DantonError(stderr)

            # Parse the result.
            if not os.path.exists(outfile): return []
            primaries = []
            for event in danton.iter_event(outfile):
                primaries.append([event.weight, event.primary.energy])
            os.remove(outfile)
            return primaries
    return sample
