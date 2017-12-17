#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This is a test procedure. It computes the density of decaying taus as
# function of the tau energy. The Waxman Bachall bound is assumed for the
# primary neutrino flux.

import json
import time

import numpy

from retro.primary import PrimarySampler
from grand_tour import Topography


def Sampler(requested):
    """Encapsulation of a primary sampler
    """
    # Settings
    card = {
        "topography": {
            "latitude": 43,
            "longitude": 87,
            "path": "flat/10"},
        "generator": [[1., {"energy": [1E+06, 1E+12]}]],
        "primary": {
            "events": 100 * requested,
            "requested": requested,
            "longitudinal": False}}

    # Configure the sampler
    topo = Topography(**card["topography"])
    sample_primaries = PrimarySampler(card["primary"], card["generator"],
                                      card["topography"], topo)

    def sample(energy, theta, altitude=0.):
        position = (0., 0., altitude)
        direction = topo.angular_to_local(position, theta, 0.)
        return sample_primaries(15, position, energy, direction)
    return sample


if __name__ == "__main__":
    WBB = 2. / 3. * 1E-04
    theta = 91.5
    sample = Sampler(1000)
    energy = numpy.logspace(6., 10., 41)

    density = numpy.zeros((len(energy), 2))
    for i, ei in enumerate(energy):
        print "# Processing {:} / {:}".format(i + 1, len(energy))
        t0 = time.time()
        primaries = sample(ei, theta)
        n = primaries[1]
        w = numpy.array([primary[0] * WBB for primary in primaries[0]])
        mu = sum(w) / n
        sigma = numpy.sqrt((sum(w**2) / n - mu**2) / n)
        density[i, :] = (mu, sigma)
        print "  density = {:.3E} +- {:.3E} GeV^-1 m^-3 sr^-1 s^-1".format(
            mu, sigma)
        print "  --> Done in {:.1f} s".format(time.time() - t0)

    with open("events-primary.json", "wb+") as f:
        json.dump({"energy": energy.tolist(), "density": density.tolist()}, f)
