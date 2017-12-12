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
    """Encapsulation for logging events in JSON format.
    """

    def __init__(self, path):
        event = {}
        self._path = path
        with open(path, "w+") as f:
            pass
        self._previous = -1

    def __call__(self, **kwargs):
        event = {"previous": self._previous}
        event.update(kwargs)
        with open(self._path, "a") as f:
            self._previous = f.tell()
            json.dump(event, f)
            f.write("\n")


class EventIterator:
    """Iterator over events stored in a JSON file."""

    def __init__(self, path):
        f = open(path)
        self._file = f
        self._prev = -1
        self._current = None

    def __del__(self):
        if self._file:
            self._file.close()
            self._file = None

    def __iter__(self):
        self.rewind()
        return self

    def next(self):
        """Pop the next event from file."""
        try:
            event = json.loads(self._file.readline())
        except ValueError:
            raise StopIteration()
        self._prev = event["previous"]
        self._current = event
        return event

    def previous(self):
        """Pop the previous event from file."""
        if self._prev >= 0:
            self._file.seek(self._prev, 0)
        else:
            self._file.seek(0, 0)
        return self.next()

    def rewind(self):
        """Rewind the file."""
        self._file.seek(0, 0)
        self._prev = -1
        self._current = None
