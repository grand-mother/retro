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

class EventLogger:
    """Encapsulation for logging events.
    """
    def __init__(self, path):
        event = {}
        self._path = path
        with open(path, "w+") as f: pass
        self._previous = -1

    def __call__(self, *args):
        event = { "previous" : self._previous, "tau_at_decay" : args }
        with open(self._path, "a") as f:
            self._previous = f.tell()
            json.dump(event, f)
            f.write("\n")
