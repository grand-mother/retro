#!/usr/bin/env python
import numpy
import pylab
from retro.event import EventIterator
from grand_tour import Topography

topo = Topography(42.928056, 86.741667, "share/topography", 25)
elevation = []
for event in EventIterator("events-ulastai.json"):
    energy, position, direction = event["tau_at_decay"]
    az, el = topo.local_to_horizontal(position, direction)
    elevation.append(el)

p, x = numpy.histogram(elevation, 40, density=True)
x = 0.5 * (x[1:] + x[:-1])
pylab.figure()
pylab.plot(x, p, "ko-")
pylab.show()
