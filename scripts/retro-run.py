#!/usr/bin/env python
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
import sys
import random
# Custom import(s)
from retro.event import EventLogger
from retro.generator import Generator
from grand_tour import Topography

def run(generator, processor, logger, topography, comment=None):
    """Generate some tau decay vertices according to the provided settings.
    """

    # Get a handle for the topography.
    topo = Topography(**topography)

    # Instanciate a generator for the taus at decay.
    generate = Generator(**generator)

    # Instanciate and event logger.
    log_event = EventLogger(**logger)

    def filter_event(energy, position, direction):
        """Raw event filter from first principles.
        """
        # First let us compute the decay length, assuming no energy loss.
        dl = energy * 4.89639E-05

        # Then let us compute the distance to the topography, propagating
        # backwards.
        dg = topo.distance(
          position, [-c for c in direction], limit=10. * dl)
        if dg is None: return 0.

        # As selection weight, let us consider the probability that no decay
        # occured on the path to rocks.
        return math.exp(-dg / dl)

    # Main loop over events.
    requested, max_trials = processor["requested"], processor["trials"]
    trials, total_trials, done = 0, 0, 0
    while True:
        # Check the termination conditions.
        if requested and (done == requested): break
        if max_trials and (total_trials >= max_trials): break

        # Generate a tentative tau decay.
        trials += 1
        total_trials += 1
        position, w0 = generate.position()
        if not topo.is_above(position): continue
        direction, w1 = generate.direction()
        energy, w2 = generate.energy()

        # Check if the generated direction is relevant considering the
        # generated position and energy.
        w3 = filter_event(energy, position, direction)
        if (w3 <= 0.) or (random.random() > w3) : continue
        weight = w0 * w1 * w2 / w3

        # TODO: decay the tau with danton, or alouette.

        # Log the decay.
        log_event(energy, position, direction, weight, trials)
        trials = 0
        done += 1

if __name__ == "__main__":
    # Read the settings from a JSON card.
    with open(sys.argv[1]) as f: settings = json.load(f)
    run(**settings)
