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

#ifndef RETRO_GENERATOR_H_
#define RETRO_GENERATOR_H_

/* Random generator models */
enum retro_generator_mode {
        RETRO_GENERATOR_MODE_UNIFORM = 0,
        RETRO_GENERATOR_MODE_LINEAR,
        RETRO_GENERATOR_MODE_1_OVER_E,
        RETRO_GENERATOR_MODE_1_OVER_E2,
        RETRO_GENERATOR_MODE_LOCAL,
        RETRO_GENERATOR_MODE_GEODETIC,
};

struct retro_generator;
typedef double generator_function_t(
    struct retro_generator * generator, double * value);

struct gt_topography;
struct roar_handler;

struct retro_product {
        int pid;
        double momentum[3];
};

struct retro_generator {
        struct roar_handler * handler;
        struct gt_topography * topography;

        enum retro_generator_mode position_mode;
        double position_weight;
        double position_parameter[5][2];
        generator_function_t * position;

        double direction_weight;
        double direction_parameter[3];
        generator_function_t * direction;

        double energy_weight;
        double energy_parameter[2];
        generator_function_t * energy;

        long long decay_state;
        int decay_pid;
        double decay_momentum[3];
        double decay_polarisation[3];
        int decay_product_n;
#define RETRO_MAX_DECAY_PRODUCTS 10
        struct retro_product decay_product[RETRO_MAX_DECAY_PRODUCTS];
        generator_function_t * decay;
};

struct retro_card;
void generator_initialise(struct retro_generator * generator,
    struct retro_card * card, struct roar_handler * handler,
    struct gt_topography * topography);

void generator_decay_initialise(struct retro_generator * generator, int pid,
    double energy, double * direction);

#endif
