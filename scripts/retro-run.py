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

import copy
import json
import math
import sys
import random
# Custom import(s)
from grand_tour import Topography
from retro import TAU_CTAU, TAU_MASS
from retro.event import EventLogger
from retro.generator import Generator
from retro.primary import PrimarySampler


def run(generator, processor, logger, topography, primary=None, antenna=None,
        preselector=None, comment=None):
    """Generate some tau decay vertices according to the provided settings.
    """

    # Get a handle for the topography.
    topo = Topography(**topography)

    # Instanciate a generator for the taus at decay. First let us vectorize
    # and flatten the options.
    if isinstance(generator, dict):
        generator = [[1., generator], ]
    for i, (w, opts) in enumerate(generator):
        if i > 0:
            new = copy.deepcopy(generator[i - 1][1])
            new.update(opts)
            generator[i][1] = new
    generate = Generator(generator, topo)

    # Instanciate an event logger.
    log_event = EventLogger(**logger)

    # Initialise the primary sampler.
    if primary:
        sample_primaries = PrimarySampler(primary, generator, topography, topo)
    else:
        def sample_primaries(pid, position, energy, direction): return []

    # Initialise the preselector
    if antenna is not None:
        from retro.preselector import Preselector
        if preselector is None:
            preselector = {}
        preselect = Preselector(topo, antenna, **preselector)

    def filter_vertex(energy, position, direction):
        """Vertex filter based on the tau decay length.
           TODO: check the grammage before as well.
        """
        # First let us compute the decay length, assuming no energy loss.
        dl = energy * TAU_CTAU / TAU_MASS

        # Then let us compute the distance to the topography, propagating
        # backwards.
        dg = topo.distance(
            position, [-c for c in direction], limit=10. * dl)
        if dg is None:
            return 0.

        # As selection weight, let us consider the probability that no decay
        # occured on the path to rocks.
        return math.exp(-dg / dl)

    # Infer the energy threshold for showers from the generation model.
    threshold = float("inf")
    for _, opts in generator:
        e0, e1 = opts["energy"]
        if isinstance(e0, basestring):
            e0, _ = e1
        if e0 < threshold:
            threshold = e0

    # Main loop over events.
    try:
        max_trials = processor["trials"]
    except KeyError:
        max_trials = None
    try:
        requested = processor["requested"]
    except KeyError:
        if max_trials is None:
            raise ValueError(
                "a requested or maximum number of events must be specified")
        requested = max_trials
    trials, neff, total_trials, done = 0, 0., 0, 0
    pid = 15
    while True:
        # Check the termination conditions.
        if requested and (done == requested):
            break
        if max_trials and (total_trials >= max_trials):
            break

        # Generate a tentative decay vertex.
        trials += 1
        total_trials += 1
        w0 = generate.model()
        position, w1 = generate.position()
        if not topo.is_above(position):
            continue
        direction, w2 = generate.direction(position)
        energy, w3 = generate.energy()

        # Check if the generated direction is relevant considering the
        # generated position and its energy.
        w4 = filter_vertex(energy, position, direction)
        neff += w4
        if (w4 == 0.) or (random.random() > w4):
            continue
        weight = w0 * w1 * w2 * w3 / w4

        # Generate a valid tau decay, i.e. with enough energy for the shower.
        while True:
            decay, state = generate.decay(pid, energy, direction)
            shower_energy = 0.
            for (pid_, momentum) in decay:
                aid = abs(pid_)
                if aid in (12, 13, 14, 16):
                    continue
                shower_energy += sum(m**2 for m in momentum)**0.5
            if shower_energy >= threshold:
                break
            trials += 1
            neff += w4

        # Preselect antennas that might detect the radio signal from the shower
        if antenna is not None:
            selection = preselect(shower_energy, position, direction)
            if len(selection) < 4:
                continue
        else:
            selection = None

        # Sample the primary flux.
        if primary and (primary["events"] > 0):
            primaries = sample_primaries(
                pid, position, energy, direction)
            if len(primaries[0]) == 0:
                continue
        else:
            primaries = ([], 0)

        # Build the tag.
        theta, phi = topo.local_to_angular(position, direction)
        if phi < 0.:
            phi += 360.
        tag = ("E.{:.0e}".format(energy * 1E+09).replace("+", ""),
               "X.{:.0f}".format(position[0]), "Y.{:.0f}".format(position[1]),
               "Z.{:.0f}".format(position[2]), "T.{:.0f}".format(theta),
               "P.{:.0f}".format(phi), "D.{:}".format(state))
        tag = "_".join(tag)

        # Log the event.
        log_event(tag=tag, tau_at_decay=(energy, position, direction),
                  decay=decay, primaries=primaries,
                  statistics=(weight, trials, neff), antennas=selection)
        trials, neff = 0, 0.
        done += 1


if __name__ == "__main__":
    # Read the settings from a JSON card.
    with open(sys.argv[1]) as f:
        settings = json.load(f)
    run(**settings)
