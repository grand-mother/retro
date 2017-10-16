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

import collections
from ctypes import *
import math
import random
import types

def _PositionGenerator(position, topo_handle):
    """Closure for generating a tentative tau decay position.
    """
    x, y, z = position
    weight = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0])
    def generate():
        xi, yi = random.uniform(*x), random.uniform(*y)
        zi = random.uniform(*z) + topo_handle.ground_altitude(xi, yi)
        return (xi, yi, zi), weight
    return generate

def _DirectionGenerator(elevation, topo_handle):
    """Closure for generating a tentative tau direction before decay.
    """
    deg = math.pi / 180.
    if isinstance(elevation[0], basestring): model, range_ = elevation
    else: model, range_ = "uniform", elevation
    c0, c1 = (math.sin(x * deg) for x in range_)
    if model == "uniform":
        weight = (c1 - c0) * 2. * math.pi
        def generate(position):
            c = random.uniform(c0, c1)
            elevation = -math.asin(c) / deg
            azimuth = random.uniform(-180, 180.)
            direction = topo_handle.horizontal_to_local(
              position, azimuth, elevation)
            return direction, weight
    elif model == "linear":
        if c0 < 0.: raise ValueError("invalid range for linear pdf")
        a, b = c0**2, (c1**2 - c0**2)
        weight = b * math.pi
        def generate(position):
            while True:
                c = (a + b * random.random())**0.5
                if c > 0.: break
            elevation = -math.asin(c) / deg
            azimuth = random.uniform(-180, 180.)
            direction = topo_handle.horizontal_to_local(
              position, azimuth, elevation)
            return direction, weight / c
    else:
        raise ValueError("invalid generation model for the direction")
    return generate

def _EnergyGenerator(energy):
    """Closure for generating a tentative tau energy before decay.
    """
    e0 = energy[0]
    lnr = math.log(energy[1] / e0)
    def generate():
        e = e0 * math.exp(random.uniform(0., lnr))
        return e, e * lnr
    return generate

def _ModelGenerator(pdf):
    """Closure for picking a generation model.
    """
    s, cdf, weight = 0., [], []
    for p in pdf:
        s += p
        cdf.append(p)
        weight.append(1. / p)

    def generate(self):
        u = random.random()
        for i, ci in enumerate(cdf):
            if u <= ci: break
        self._current = self._models[i]
        return weight[i]
    return generate

# Handle for the ALOUETTE C library.
_alouette = None

class AlouetteError(Exception):
    """Custom exception for an ALOUETTE library error.
    """
    def __init__(self, rc):
        message = _alouette.alouette_strerror(rc)
        super(AlouetteError, self).__init__(
          "alouette error ({:})".format(message))
        self.rc = rc

def _DecayGenerator():
    """Closure for decaying a tau.
    """
    global _alouette
    if _alouette is None:
        # Load the C library.
        lib = cdll.LoadLibrary("libalouette.so")

        # Prototypes.
        lib.alouette_initialise.argtypes = (c_int, POINTER(c_int))
        lib.alouette_strerror.argtypes = (c_int,)
        lib.alouette_strerror.restype = c_char_p
        lib.alouette_decay.argtypes = (
          c_int, POINTER(c_double), POINTER(c_double))
        lib.alouette_product.argtypes = (POINTER(c_int), POINTER(c_double))

        # Initialise with a random seed.
        seed = c_int(random.randint(0, 2**32 - 1))
        rc = lib.alouette_initialise(1, byref(seed))
        if rc != 0: raise AlouetteError(rc)

        _alouette = lib

    def generate(self, pid, energy, direction):
        m = 1.77682
        p = math.sqrt((energy - m) * (energy + m))
        momentum = (3 * c_double)(*[x * p for x in direction])
        polarisation = (3 * c_double)(*direction)
        rc = _alouette.alouette_decay(pid, momentum, polarisation)
        if rc != 0: raise AlouetteError(rc)
        pid = c_int(0)
        products = []
        while True:
            rc = _alouette.alouette_product(byref(pid), momentum)
            if rc != 0: break
            products.append([pid.value, [x for x in momentum]])
        return products
    return generate

class Generator(object):
    """Generator for tau vertices, before decay."""
    def __init__(self, generator, topo_handle):
        # Build the generation routines.
        total = sum(g[0] for g in generator)
        models, pdf = [], []
        GenerationModel = collections.namedtuple("GenerationModel",
          ("position", "direction", "energy"))
        for weight, opts in generator:
            pdf.append(weight / total)
            models.append(GenerationModel(
                _PositionGenerator(opts["position"], topo_handle),
                _DirectionGenerator(opts["elevation"], topo_handle),
                _EnergyGenerator(opts["energy"])))
        self._models = models
        self._current = models[0]

        # Overwrite the pick and decay methods.
        self.model = types.MethodType(_ModelGenerator(pdf), self)
        self.decay = types.MethodType(_DecayGenerator(), self)

    def model(self):
        """Randomly select a generation model.
        """
        pass

    def position(self):
        """Generate a tentative tau decay position.
        """
        return self._current.position()

    def direction(self, position):
        """Generate a tentative tau direction before decay.
        """
        return self._current.direction(position)

    def energy(self):
        """Generate a tentative tau energy before decay..
        """
        return self._current.energy()

    def decay(self, pid, energy, direction):
        """Generate a tentative tau decay.
        """
        pass
