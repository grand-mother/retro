#!/usr/bin/env python
import json
import math
import sys
import random
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
        c1, c0 = (math.sin(x * deg) for x in generator["elevation"])
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

    def EventLogger():
        """Closure for logging events.
        """
        event = {}
        path = logger["path"]
        with open(path, "w+") as f: pass
        def log(*args):
            event["tau_at_decay"] = args
            with open(path, "a") as f:
                json.dump(event, f)
                f.write("\n")
        return log
    log_event = EventLogger()

    # Main loop over events.
    trials = 0
    while True:
        # Generate a tentative tau decay.
        trials += 1
        position, w0 = generate_position()
        if not topo.is_above(position): continue
        direction, w1 = generate_direction()
        energy, w2 = generate_energy()

        # TODO: check if the generated direction is relevant considering the
        # generated position and energy.
        weight = w0 * w1 * w2

        # TODO: decay the tau with danton, or alouette.

        # Log the decay.
        log_event(energy, position, direction, weight, trials)
        trials = 0

if __name__ == "__main__":
    with open(sys.argv[1]) as f: settings = json.load(f)
    run(**settings)
