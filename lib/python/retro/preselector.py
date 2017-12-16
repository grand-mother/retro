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
import numpy


def Preselector(topography, antenna, check_xmax=True, shadowing=True):
    """Closure for preselecting antennas at sight of the shower
    """
    # Load the antenna positions
    with open(antenna["position"], "rb") as f:
        ra = numpy.array(json.load(f))

    # Discretization accuracy for rays, in m
    deltar = 200.

    def preselect(shower_energy, position, direction):
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

    return preselect
