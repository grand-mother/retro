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
import types

from retro import TAU_CTAU, TAU_MASS


def _AntennaSelector(topography, antenna, check_xmax=True, shadowing=True):
    """Closure for selecting antennas at sight of the shower
    """
    # Import numpy only if the antenna selector is required
    import numpy

    # Load the antenna positions
    with open(antenna["position"], "rb") as f:
        ra = numpy.array(json.load(f))

    # Discretization accuracy for rays, in m
    deltar = 200.

    def select(self, shower_energy, position, direction):
        # Cone parameters
        gamma = numpy.deg2rad(3.)
        zcmin = 14E+03  # m
        zcmax = 165E+03 * shower_energy / 1E+09 + 55E+03  # m

        # Check if the shower crashes into a mountain early, before xmax
        if check_xmax:
            s = numpy.arange(0., zcmin + deltar, deltar)
            xs, ys, zs = [position[i] + direction[i] * s for i in xrange(3)]
            zg = [topography.ground_altitude(xi, yi) for xi, yi in zip(xs, ys)]
            if (zs <= zg).any():
                return []

        # Select the antenna(s) within the cone
        dr = ra - position
        zp = numpy.dot(dr, direction)
        rp2 = numpy.sum(dr**2, axis=1) - zp**2
        test_radius = rp2 <= ((zp - zcmin) * numpy.tan(gamma))**2
        test_edgemin = zp >= zcmin
        test_edgemax = zp <= zcmax
        index = numpy.nonzero(test_radius & test_edgemin & test_edgemax)[0]
        if len(index) == 0:
            return index

        if not shadowing:
            return ra[index, :].tolist()

        # Check for shadowing
        r0 = position + zcmin * numpy.array(direction)

        def check_shadowing(i):
            u = ra[i, :] - r0
            d = numpy.linalg.norm(u)
            u /= d
            s = numpy.arange(0., d, deltar)
            xj, yj, zj = [r0[j] + u[j] * s for j in xrange(3)]
            zg = [topography.ground_altitude(x, y) for x, y in zip(xj, yj)]
            if (zj <= zg).any():
                return False
            return True

        K = filter(check_shadowing, index)
        return ra[K, :].tolist()

    return select


def _VertexSelector(topography, limit):
    """Closure for selecting decay vertices
    """

    def select(self, energy, position, direction):
        """Vertex filter based on the tau decay length.
        """
        # First let us compute the decay length, assuming no energy loss.
        dl = energy * TAU_CTAU / TAU_MASS

        # Then let us compute the distance to the topography, propagating
        # backwards.
        dg = topography.distance(
            position, [-c for c in direction], limit=limit * dl)
        if dg is None:
            return math.exp(-limit)

        # As selection weight, let us consider the probability that no decay
        # occured on the path to rocks.
        return math.exp(-dg / dl)

    return select


class Selector(object):
    """Selector for tau events."""

    def __init__(self, selector, topo_handle, antenna):
        # Overwrite the selection methods, if relevant
        if (selector is None) or not selector:
            return

        if ((antenna is not None) and
                ((antenna not in selector.keys()) or selector["antenna"])):
            if ((antenna not in selector.keys()) or
                    not isinstance(selector["antenna"], dict)):
                opts = {}
            else:
                opts = selector["antenna"]
            self.antennas = types.MethodType(
                _AntennaSelector(topo_handle, antenna, **opts), self)
        if ("vertex" in selector.keys()) and (selector["vertex"]["limit"] > 0.):
            self.vertex_weight = types.MethodType(
                _VertexSelector(topo_handle, selector["vertex"]["limit"]), self)

    def antennas(self, shower_energy, position, direction):
        """Select antennas at sight of the decay
        """
        return None

    def vertex_weight(self, energy, position, direction):
        """Compute the vertex selection weight
        """
        return 1.
