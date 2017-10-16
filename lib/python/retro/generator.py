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

import math
import random
import types

def _PositionGenerator(position):
    """Closure for generating a tentative tau decay position.
    """
    x, y, z = position
    weight = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0])
    def generate(self):
        return (random.uniform(*x), random.uniform(*y),
                random.uniform(*z)), weight
    return generate

def _DirectionGenerator(elevation):
    """Closure for generating a tentative tau direction before decay.
    """
    deg = math.pi / 180.
    c0, c1 = (math.sin(x * deg) for x in elevation)
    weight = (c1 - c0) * 2. * math.pi
    def generate(self):
        c = random.uniform(c0, c1)
        s = math.sqrt(1. - c * c)
        phi = random.uniform(0., 2. * math.pi)
        return (-s * math.cos(phi), -s*math.sin(phi), -c), weight
    return generate

def _EnergyGenerator(energy):
    """Closure for generating a tentative tau energy before decay.
    """
    e0 = energy[0]
    lnr = math.log(energy[1] / e0)
    def generate(self):
        e = e0 * math.exp(random.uniform(0., lnr))
        return e, e * lnr
    return generate

class Generator(object):
    """Generator for tau vertices, before decay."""
    def __init__(self, position, elevation, energy):
        # Implement the generation routines.
        self.position = types.MethodType(_PositionGenerator(position), self)
        self.direction = types.MethodType(_DirectionGenerator(elevation), self)
        self.energy = types.MethodType(_EnergyGenerator(energy), self)

    def position(self):
        """Generate a tentative tau decay position.
        """
        pass

    def direction(self):
        """Generate a tentative tau direction before decay.
        """
        pass

    def energy(self):
        """Generate a tentative tau energy before decay..
        """
        pass
