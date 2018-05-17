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
typedef double selector_function_t(struct retro_selector * selector,
    double energy, const double * position, const double * direction);

struct roar_handler;
struct gt_topography;

enum setup_cone_model {
        SETUP_CONE_MODEL_NONE = 0,
        SETUP_CONE_MODEL_3DEG,
        SETUP_CONE_MODEL_AGRESSIVE
};

enum setup_array_model {
        SETUP_ARRAY_MODEL_NONE = 0,
        SETUP_ARRAY_MODEL_FILE,
        SETUP_ARRAY_MODEL_WORLD_WIDE
};

struct retro_selector {
        struct roar_handler * handler;
        struct gt_topography * topography;

        double vertex_limit;
        selector_function_t * vertex;

        enum setup_cone_model setup_cone;
        int setup_xmax;
        int setup_shadowing;
        enum setup_array_model array;
        double setup_ww_step;
        double setup_ww_height;
        int setup_n;
        int setup_size;
        double * setup_data;
        double ** setup_selection;
        selector_function_t * setup;
        selector_function_t * setup_gamma;
};

struct retro_card;
void selector_initialise(struct retro_selector * selector,
    const struct retro_card * card, struct roar_handler * handler,
    struct gt_topography * topography);

void selector_finalise(struct retro_selector * selector);

#endif
