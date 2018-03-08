/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * Radio nEuTRino simulatiOn (RETRO).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

#include "random.h"
#include "danton.h"

static struct danton_context * context;

void random_initialise(struct danton_context * context_) { context = context_; }

double random_uniform01(void) { return danton_get_uniform01(context); }
