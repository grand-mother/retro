#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
#  Author: Valentin NIESS (niess@in2p3.fr)
#
#  Radio nEuTRino simulatiOn (RETRO).
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>

import json
import math
import sys
import random
# Custom import(s)
from grand_tour import Topography

def run(generator, logger, topography, comment=None):
    """Generate some tau decay vertices according to the provided settings.
    """

    # Get a handle for the topography.
    topo = Topography(**topography)

    def PositionGenerator():
        """Closure for generating a tentative tau decay position.
        """
        x, y, z = generator["position"]
        weight = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0])
        def generate():
            return (random.uniform(*x), random.uniform(*y),
                    random.uniform(*z)), weight
        return generate
    generate_position = PositionGenerator()

    def DirectionGenerator():
        """Closure for generating a tentative tau direction before decay.
        """
        deg = math.pi / 180.
        c0, c1 = (math.sin(x * deg) for x in generator["elevation"])
        weight = (c1 - c0) * 2. * math.pi
        def generate():
            c = random.uniform(c0, c1)
            s = math.sqrt(1. - c * c)
            phi = random.uniform(0., 2. * math.pi)
            return (-s * math.cos(phi), -s*math.sin(phi), -c), weight
        return generate
    generate_direction = DirectionGenerator()

    def EnergyGenerator():
        """Closure for generating a tentative tau energy before decay.
        """
        e0 = generator["energy"][0]
        lnr = math.log(generator["energy"][1] / e0)
        def generate():
            e = e0 * math.exp(random.uniform(0., lnr))
            return e, e * lnr
        return generate
    generate_energy = EnergyGenerator()

    class EventLogger:
        """Encapsulation for logging events.
        """
        def __init__(self, path):
            event = {}
            self._path = path
            with open(path, "w+") as f: pass
            self._previous = -1

        def __call__(self, *args):
            event = { "previous" : self._previous, "tau_at_decay" : args }
            with open(self._path, "a") as f:
                self._previous = f.tell()
                json.dump(event, f)
                f.write("\n")
    log_event = EventLogger(**logger)

    def filter_event(energy, position, direction):
        """Raw event filter from first principles.
        """
        # First let us compute the decay length, assuming no energy loss.
        dl = energy * 4.89639E-05

        # Then let us compute the distance to the topography, propagating
        # backwards.
        dg = topo.distance(
          position, tuple([-c for c in direction]), limit=10. * dl)
        if dg is None: return 0.

        # As selection weight, let us consider the probability that no decay
        # occured on the path to rocks.
        return math.exp(-dg / dl)

    # Main loop over events.
    trials = 0
    while True:
        # Generate a tentative tau decay.
        trials += 1
        position, w0 = generate_position()
        if not topo.is_above(position): continue
        direction, w1 = generate_direction()
        energy, w2 = generate_energy()

        # Check if the generated direction is relevant considering the
        # generated position and energy.
        w3 = filter_event(energy, position, direction)
        if (w3 <= 0.) or (random.random() > w3) : continue
        weight = w0 * w1 * w2 / w3

        # TODO: decay the tau with danton, or alouette.

        # Log the decay.
        log_event(energy, position, direction, weight, trials)
        trials = 0

if __name__ == "__main__":
    # Read the settings from a JSON card.
    with open(sys.argv[1]) as f: settings = json.load(f)
    run(**settings)
