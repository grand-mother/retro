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
from retro import TAU_BR_MU
from retro.event import EventLogger
from retro.generator import Generator
from retro.primary import PrimarySampler
from retro.selector import Selector


def run(generator, processor, logger, topography, primary=None, setup=None,
        selector=None, comment=None):
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

    # Initialise the selector
    select = Selector(selector, topo, setup)

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
    trials, total_trials, done = 0, 0, 0
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
        w4 = select.vertex_weight(energy, position, direction)
        if (w4 == 0.) or (random.random() > w4):
            continue
        weight = w0 * w1 * w2 * w3 / w4

        # Generate a valid tau decay, i.e. not a muonic decay.
        shower_energy = 0.
        while shower_energy == 0.:
            decay, state = generate.decay(pid, energy, direction)
            for (pid_, momentum) in decay:
                aid = abs(pid_)
                if aid in (12, 13, 14, 16):
                    continue
                shower_energy += sum(m**2 for m in momentum)**0.5
        if shower_energy < threshold:
            continue
        weight *= 1. - TAU_BR_MU

        # Preselect antennas that might detect the radio signal from the shower
        selection = select.antennas(shower_energy, position, direction)
        if (selection is not None) and len(selection) < 4:
            continue

        # Sample the primary flux.
        if primary and (primary["events"] > 0):
            primaries, primary_trials = sample_primaries(
                pid, position, energy, direction)
            if len(primaries) == 0:
                continue
            for i, (wi, ei, _, _) in enumerate(primaries):
                primaries[i][0] = weight * wi * ei**2
        else:
            primaries, primary_trials = [], 0

        # Build the tag.
        lla = topo.local_to_lla(position)
        latitude, longitude, altitude = lla
        hz = topo.local_to_angular(position, direction)
        theta, phi = hz
        if phi < 0.:
            phi += 360.
        tag = ("E.{:.0e}".format(energy * 1E+09).replace("+", ""),
               "Z.{:.0f}".format(theta), "A.{:.0f}".format(phi),
               "La.{:.0f}".format(latitude), "Lo.{:.0f}".format(longitude),
               "H.{:.0f}".format(altitude), "D.{:}".format(state))
        tag = "_".join(tag)

        # Log the event.
        tau_at_decay = (weight, energy, position, direction, lla, hz)
        log_event(tag=tag, tau_at_decay=tau_at_decay, decay=decay,
                  primaries=primaries, statistics=(trials, primary_trials),
                  antennas=selection, origin=(topography["latitude"],
                                              topography["longitude"]))
        trials = 0
        done += 1


if __name__ == "__main__":
    # Read the settings from a JSON card.
    with open(sys.argv[1]) as f:
        settings = json.load(f)
    run(**settings)
