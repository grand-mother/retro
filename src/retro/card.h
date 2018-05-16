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

#ifndef RETRO_CARD_H_
#define RETRO_CARD_H_

#include "generator.h"

/* Data structure for a RETRO card */
struct retro_card {
        struct jsmn_tea * tea;
        const char * path;

        enum retro_generator_mode generator_position_mode;
        enum retro_generator_mode generator_theta_mode;
        enum retro_generator_mode generator_energy_mode;
        double generator_theta[2];
        double generator_energy[2];
        double generator_position[3][2];

        long long processor_requested;
        long long processor_trials;

        double selector_vertex_limit;
        int selector_setup_cone;
        int selector_setup_xmax;
        int selector_setup_shadowing;

        char * logger_path;

        double topography_latitude;
        double topography_longitude;
        double topography_density;
        char * topography_path;
        int topography_stack_size;

        int primary_events;
        int primary_requested;
        int primary_longitudinal;

        char * setup_path;
};

/* Initialise a RETRO card with default values */
void retro_card_initialise(struct retro_card * card);

/* Properly destroy a RETRO card */
void retro_card_destroy(struct retro_card * card);

/* Update a RETRO card from a JSON file */
struct roar_handler;
void retro_card_update(
    struct retro_card * card, struct roar_handler * handler, char * path);

#endif
