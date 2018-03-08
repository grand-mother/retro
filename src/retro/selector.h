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

#ifndef RETRO_SELECTOR_H_
#define RETRO_SELECTOR_H_

struct retro_selector;
typedef double selector_function_t(const struct retro_selector * selector,
    double energy, const double * position, const double * direction);

struct roar_handler;
struct gt_topography;

struct retro_selector {
        struct roar_handler * handler;
        struct gt_topography * topography;

        double vertex_limit;
        selector_function_t * vertex;

        int setup_cone;
        int setup_xmax;
        int setup_shadowing;
        int setup_n;
        int setup_size;
        double * setup_data;
        double ** setup_selection;
        selector_function_t * setup;
};

struct retro_card;
void selector_initialise(struct retro_selector * selector,
    const struct retro_card * card, struct roar_handler * handler,
    struct gt_topography * topography);

void selector_finalise(struct retro_selector * selector);

#endif
