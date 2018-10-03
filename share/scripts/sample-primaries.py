#!/usr/bin/env python
import sys
sys.path.append("lib/python")
from retro.event import EventIterator

import argparse
import ctypes as C
import json
import signal


def DantonSampler(primary_exponent):
    """Closure providing a DANTON primary sampler"""
    lib = C.cdll.LoadLibrary("deps/danton/lib/libdanton.so")

    # Configure the error handling
    lib.danton_error_count.argtypes = (C.c_void_p,)
    lib.danton_error_count.restype = C.c_int
    lib.danton_error_pop.argtypes = (C.c_void_p,)
    lib.danton_error_pop.restype = C.c_char_p

    def Define(function, arguments = None, result = None, contextual=False):
        """Export a DANTON library function"""
        f = getattr(lib, function)
        if arguments is not None:
            f.argtypes = arguments
        if result is not None:
            f.restype = result

        def call(*args):
            """Encapsulation of the function with error check"""
            r = f(*args)
            if contextual:
                c = context
            else:
                c = None
            if lib.danton_error_count(c) > 0:
                raise RuntimeError(lib.danton_error_pop(c))
            return r

        # Export the encapsulated function
        globals()[function] = call

    # Initialise the library
    Define("danton_initialise",
        arguments = (C.c_char_p, C.c_char_p, C.c_char_p, C.c_void_p,
                     C.c_void_p))

    danton_initialise(None, None, None, None, None)

    # Utility functions for index conversion
    Define("danton_particle_index", arguments = (C.c_int,), result = C.c_uint)
    Define("danton_particle_pdg", arguments = (C.c_uint,), result = C.c_int)

    # Configure the geometry
    Define("danton_earth_model",
        arguments = (C.c_char_p, C.c_char_p, C.c_int, C.c_char_p, C.c_double,
                     C.POINTER(C.c_int)))

    sea = C.c_int(0)
    danton_earth_model("WGS84", "share/SRTMGL1", 0, "Rock", 2.65E+03,
        C.byref(sea))

    # Set the particle sampler
    class Sampler(C.Structure):
        _fields_ = (("latitude", C.c_double), ("longitude", C.c_double),
            ("altitude", 2 * C.c_double), ("azimuth", 2 * C.c_double),
            ("elevation", 2 * C.c_double), ("energy", 2 * C.c_double),
            ("weight", 8 * C.c_double))

    Define("danton_sampler_create", result = C.POINTER(Sampler))
    Define("danton_sampler_update", arguments = (C.POINTER(Sampler),),
        contextual=True)

    sampler = danton_sampler_create()
    weight = (8 * C.c_double)(*(8 * [0.,]))
    weight[danton_particle_index(15)] = 1.
    sampler.contents.weight = weight

    # Set the primary model
    class Primary(C.Structure):
        _fields_ = (("flux", C.c_void_p), ("energy", 2 * C.c_double))

    Define("danton_powerlaw_create",
        arguments = (C.c_double, C.c_double, C.c_double, C.c_double),
        result = C.POINTER(Primary))

    primary = danton_powerlaw_create(1E+06, 1E+14, primary_exponent, 1.)

    # Set the event recorder
    class State(C.Structure):
        _fields_ = (("pid", C.c_int), ("energy", C.c_double), ("position", 3 *
        C.c_double), ("direction", 3 * C.c_double))

    class Event(C.Structure):
        _fields_ = (("id", C.c_long), ("weight", C.c_double),
            ("primary", C.POINTER(State)), ("generation", C.c_int),
            ("vertex", C.POINTER(State)), ("final", C.POINTER(State)),
            ("n_products", C.c_int), ("product", C.c_void_p))

    danton_event_cb = C.CFUNCTYPE(
        C.c_int, C.c_void_p, C.c_void_p, C.POINTER(Event))

    class Recorder(C.Structure):
        _fields_ = (("record_event", danton_event_cb),
            ("record_grammage", C.c_void_p), ("events", C.c_long))

    primaries_container = []

    @danton_event_cb
    def record_event(context, recorder, event):
        r = C.cast(recorder, C.POINTER(Recorder)).contents
        weight = event.contents.weight / r.events
        primary = event.contents.primary.contents
        primaries_container.append((weight, primary.energy))
        return 0

    recorder = Recorder(record_event, None, 0)

    # Set the simulation context
    class Context(C.Structure):
        _fields_ = (("mode", C.c_uint), ("longitudinal", C.c_int),
            ("decay", C.c_int), ("primary", 6 * C.POINTER(Primary)),
            ("sampler", C.POINTER(Sampler)), ("recorder", C.POINTER(Recorder)),
            ("run_action", C.c_void_p))

    Define("danton_context_create", result = C.POINTER(Context))

    context = danton_context_create()
    context.contents.mode = 0
    context.contents.longitudinal = 1
    context.contents.decay = 0
    primaries = (6 * C.POINTER(Primary))(*(6 * [None,]))
    primaries[danton_particle_index(16)] = primary
    context.contents.primary = primaries
    context.contents.sampler = sampler
    context.contents.recorder = C.pointer(recorder)

    Define("danton_run", arguments = (C.POINTER(Context), C.c_long, C.c_long),
        contextual=True)

    def sample(events, energy, latitude, longitude, altitude, azimuth,
        elevation):
        # Configure the sampler
        s = sampler.contents
        s.energy[0] = s.energy[1] = C.c_double(energy)
        s.latitude = C.c_double(latitude)
        s.longitude = C.c_double(longitude)
        s.altitude[0] = s.altitude[1] = C.c_double(altitude)
        s.azimuth[0] = s.azimuth[1] = C.c_double(azimuth)
        s.elevation[0] = s.elevation[1] = C.c_double(elevation)
        danton_sampler_update(sampler)

        # Configure the recorder
        recorder.events = events

        # Run DANTON
        del primaries_container[:]
        danton_run(context, 1000 * events, events)
        return primaries_container

    return sample

# Get the primary sampler
sample_primaries = DantonSampler(primary_exponent=-2.)


if __name__ == "__main__":
    # Unpack the arguments
    parser = argparse.ArgumentParser(
        description="Backward sample primary neutrinos given tau events.")
    parser.add_argument("files", metavar="FILE", type = str, nargs = "+",
        help = "a file containing RETRO tau events")
    parser.add_argument("--flatten", dest="flatten", action="store_const",
                    const = True, default = False,
                    help = "flatten the result (default: order by tau event)")
    parser.add_argument("-e", "--events", dest="events", action="store",
                    type = int, default = 100,
                    help = "the requested number of neutrinos per tau event "
                           "(default: 100)")
    parser.add_argument("-s", "--skip", dest="skip", action="store",
                    type = int, default = -1,
                    help = "ignore the first SKIP events. By default no "
                           "event is ignored")
    parser.add_argument("-m", "--max", dest="max", action="store",
                    type = int, default = 0,
                    help = "the maximum number of events to process. By "
                           "default all events are processed")
    parser.add_argument("-o", "--output", dest="output", action="store",
                    type = str, default = None,
                    help = "specify an output file for the sampled neutrinos")
    args = parser.parse_args()

    # Configure the the signal handler for Ctrl+C
    def signal_handler(sig, frame):
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)

    # Set the output stream
    if args.output is not None:
        with open(args.output, "w+") as stream:
            pass
    else:
        stream = sys.stdout

    # Loop over events
    done = 0
    for filename in args.files:
        for i, event in enumerate(EventIterator(filename)):
            if i < args.skip:
                continue

            # Unpack the tau event
            tau = event["tau_at_decay"]
            energy = tau[1]
            latitude, longitude, altitude = tau[4]
            azimuth, elevation = -tau[5][1], 90. - tau[5][0]
            if azimuth < 0: azimuth += 360

            # Sample the primaries
            primaries = sample_primaries(args.events, energy, latitude,
                longitude, altitude, azimuth, elevation)

            # Dump the result
            if args.output is not None:
                stream = open(args.output, "a")
            else:
                stream = sys.stdout

            if args.flatten:
                for primary in primaries:
                    stream.write("{:.5E} {:.5E} {:.5E}\n".format(
                        primary[0], primary[1], energy))
            else:
                stream.write(
                    json.dumps((filename, i, energy, primaries)))
                stream.write("\n")
            if args.output is not None:
                stream.close()

            # Check for termination
            done += 1
            if done == args.max:
                break
