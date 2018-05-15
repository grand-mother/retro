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

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "jsmn-tea.h"
#include "retro/card.h"
#include "retro/selector.h"

/* Callback for catching errors. */
static int catch_error(
    struct roar_handler * handler, roar_function_t * referent, int code)
{
        return EXIT_SUCCESS;
}

/* Initialise a RETRO card with default values */
void retro_card_initialise(struct retro_card * card)
{
        memset(card, 0x0, sizeof(*card));
        card->tea = NULL;
        card->path = NULL;

        card->generator_theta_mode = RETRO_GENERATOR_MODE_UNIFORM;
        card->generator_energy_mode = RETRO_GENERATOR_MODE_1_OVER_E;
        card->generator_theta[0] = 80.;
        card->generator_theta[1] = 100.;
        card->generator_energy[0] = 1E-06;
        card->generator_energy[1] = 1E-12;
        card->processor_requested = -1;
        card->selector_setup_cone = SETUP_CONE_MODEL_3DEG;
        card->selector_setup_xmax = 1;
        card->selector_setup_shadowing = 1;
        card->logger_path = NULL;
        card->topography_density = 2.65E+03;
        card->topography_path = NULL;
        card->topography_stack_size = 1;
        card->primary_requested = -1;
        card->primary_longitudinal = 1;
        card->setup_path = NULL;
}

/* Properly destroy a RETRO card */
void retro_card_destroy(struct retro_card * card)
{
        jsmn_tea_destroy(&card->tea);

        free(card->logger_path);
        card->logger_path = NULL;
        free(card->topography_path);
        card->topography_path = NULL;
        free(card->setup_path);
        card->setup_path = NULL;
}

static void raise_error_key(struct retro_card * card, const char * key)
{
        ROAR_ERRNO_FORMAT(card->tea->handler, &retro_card_update, EINVAL,
            "[%s #%d] invalid key `%s`", card->path, card->tea->index, key);
}

static void raise_error_memory(struct retro_card * card)
{
        ROAR_ERRNO(card->tea->handler, &retro_card_update, ENOMEM);
}

static void raise_error_mode(struct retro_card * card, const char * mode)
{
        ROAR_ERRNO_FORMAT(card->tea->handler, &retro_card_update, EINVAL,
            "[%s #%d] unknown mode `%s`", card->path, card->tea->index, mode);
}

static void raise_error_size(struct retro_card * card, int size)
{
        ROAR_ERRNO_FORMAT(card->tea->handler, &retro_card_update, EINVAL,
            "[%s #%d] invalid array size `%d`", card->path, card->tea->index,
            size);
}

static void parse_double2(struct retro_card * card, double * value)
{
        int size;
        jsmn_tea_next_array(card->tea, &size);
        if (size != 2) raise_error_size(card, size);
        jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_DOUBLE, value);
        jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_DOUBLE, value + 1);
}

static void parse_string(struct retro_card * card, char ** string)
{
        char * s;
        jsmn_tea_next_string(card->tea, 0, &s);
        free(*string);
        if (s == NULL) {
                *string = NULL;
                return;
        }
        const int n = strlen(s) + 1;
        *string = malloc(n);
        if (*string == NULL) raise_error_memory(card);
        memcpy(*string, s, n);
}

static void parse_mode_double2(
    struct retro_card * card, char ** mode, double * value)
{
        int size;
        jsmn_tea_next_array(card->tea, &size);
        if (size != 2) raise_error_size(card, size);
        *mode = NULL;
        card->tea->handler->pre = catch_error;
        int jrc = jsmn_tea_next_string(card->tea, 0, mode);
        card->tea->handler->pre = NULL;
        if (jrc == JSMN_SUCCESS) {
                parse_double2(card, value);
        } else {
                jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_DOUBLE, value);
                jsmn_tea_next_number(
                    card->tea, JSMN_TEA_TYPE_DOUBLE, value + 1);
        }
}

#define LOOP_OVER_KEYS_BEGIN                                                   \
        int i;                                                                 \
        for (jsmn_tea_next_object(card->tea, &i); i; i--) {                    \
                char * tag;                                                    \
                jsmn_tea_next_string(card->tea, 1, &tag);

#define LOOP_OVER_KEYS_END                                                     \
        else { raise_error_key(card, tag); }                                   \
        }

static void update_generator(struct retro_card * card)
{
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "energy") == 0) {
                char * mode;
                parse_mode_double2(card, &mode, card->generator_energy);
                if (mode != NULL) {
                        if (strcmp(mode, "uniform") == 0) {
                                card->generator_energy_mode =
                                    RETRO_GENERATOR_MODE_UNIFORM;
                        } else if (strcmp(mode, "1 / E") == 0) {
                                card->generator_energy_mode =
                                    RETRO_GENERATOR_MODE_1_OVER_E;
                        } else if (strcmp(mode, "1 / E**2") == 0) {
                                card->generator_energy_mode =
                                    RETRO_GENERATOR_MODE_1_OVER_E2;
                        } else {
                                raise_error_mode(card, mode);
                        }
                }
        } else if (strcmp(tag, "position") == 0) {
                int size;
                jsmn_tea_next_array(card->tea, &size);
                if (size != 3) raise_error_size(card, size);
                for (; size > 0; size--)
                        parse_double2(
                            card, &card->generator_position[3 - size][0]);
        } else if (strcmp(tag, "theta") == 0) {
                char * mode;
                parse_mode_double2(card, &mode, card->generator_theta);
                if (mode != NULL) {
                        if (strcmp(mode, "uniform") == 0) {
                                card->generator_theta_mode =
                                    RETRO_GENERATOR_MODE_UNIFORM;
                        } else if (strcmp(mode, "linear") == 0) {
                                card->generator_theta_mode =
                                    RETRO_GENERATOR_MODE_LINEAR;
                        } else {
                                raise_error_mode(card, mode);
                        }
                }
        }
        LOOP_OVER_KEYS_END
}

static void update_processor(struct retro_card * card)
{
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "requested") == 0) {
                jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_LONG_LONG,
                    &card->processor_requested);
        } else if (strcmp(tag, "trials") == 0) {
                card->tea->handler->pre = catch_error;
                int jrc = jsmn_tea_next_null(card->tea);
                card->tea->handler->pre = NULL;
                if (jrc == JSMN_SUCCESS) {
                        card->processor_trials = -1;
                } else {
                        jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_LONG_LONG,
                            &card->processor_trials);
                }
        }
        LOOP_OVER_KEYS_END
}

static void update_selector_vertex(struct retro_card * card)
{
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "limit") == 0) {
                jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_DOUBLE,
                    &card->selector_vertex_limit);
        }
        LOOP_OVER_KEYS_END
}

static void update_selector_setup(struct retro_card * card)
{
        card->tea->handler->pre = catch_error;
        int enable;
        int jrc = jsmn_tea_next_bool(card->tea, &enable);
        card->tea->handler->pre = NULL;
        if (jrc == JSMN_SUCCESS) {
                if (enable) {
                        card->selector_setup_cone = SETUP_CONE_MODEL_3DEG;
                        card->selector_setup_xmax = 1;
                        card->selector_setup_shadowing = 1;
                } else {
                        card->selector_setup_cone = SETUP_CONE_MODEL_NONE;
                        card->selector_setup_xmax = 0;
                        card->selector_setup_shadowing = 0;
                }
                return;
        }

        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "cone") == 0) {
                card->tea->handler->pre = catch_error;
                int jrc =
                    jsmn_tea_next_bool(card->tea, &card->selector_setup_cone);
                card->tea->handler->pre = NULL;
                if (jrc != JSMN_SUCCESS) {
                        char * s;
                        jsmn_tea_next_string(card->tea, 0, &s);
                        if (strcmp(s, "3deg") == 0) {
                                card->selector_setup_cone =
                                    SETUP_CONE_MODEL_3DEG;
                        } else if (strcmp(s, "agressive") == 0) {
                                card->selector_setup_cone =
                                    SETUP_CONE_MODEL_AGRESSIVE;
                        } else {
                                raise_error_mode(card, s);
                        }
                }
        } else if (strcmp(tag, "xmax") == 0) {
                jsmn_tea_next_bool(card->tea, &card->selector_setup_xmax);
        } else if (strcmp(tag, "shadowing") == 0) {
                jsmn_tea_next_bool(card->tea, &card->selector_setup_shadowing);
        }
        LOOP_OVER_KEYS_END
}

static void update_selector(struct retro_card * card)
{
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "vertex") == 0) {
                update_selector_vertex(card);
        } else if (strcmp(tag, "setup") == 0) {
                update_selector_setup(card);
        }
        LOOP_OVER_KEYS_END
}

static void update_logger(struct retro_card * card)
{
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "path") == 0) {
                parse_string(card, &card->logger_path);
        }
        LOOP_OVER_KEYS_END
}

static void update_primary(struct retro_card * card)
{
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "requested") == 0) {
                jsmn_tea_next_number(
                    card->tea, JSMN_TEA_TYPE_INT, &card->primary_requested);
        } else if (strcmp(tag, "events") == 0) {
                card->tea->handler->pre = catch_error;
                int jrc = jsmn_tea_next_null(card->tea);
                card->tea->handler->pre = NULL;
                if (jrc == JSMN_SUCCESS) {
                        card->primary_events = -1;
                } else {
                        jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_INT,
                            &card->primary_events);
                }
        } else if (strcmp(tag, "longitudinal") == 0) {
                jsmn_tea_next_bool(card->tea, &card->primary_longitudinal);
        }
        LOOP_OVER_KEYS_END
}

static void update_topography(struct retro_card * card)
{
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "latitude") == 0) {
                jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_DOUBLE,
                    &card->topography_latitude);
        } else if (strcmp(tag, "longitude") == 0) {
                jsmn_tea_next_number(card->tea, JSMN_TEA_TYPE_DOUBLE,
                    &card->topography_longitude);
        } else if (strcmp(tag, "density") == 0) {
                jsmn_tea_next_number(
                    card->tea, JSMN_TEA_TYPE_DOUBLE, &card->topography_density);
        } else if (strcmp(tag, "path") == 0) {
                parse_string(card, &card->topography_path);
        } else if (strcmp(tag, "stack_size") == 0) {
                jsmn_tea_next_number(
                    card->tea, JSMN_TEA_TYPE_INT, &card->topography_stack_size);
        }
        LOOP_OVER_KEYS_END
}

static void update_setup(struct retro_card * card)
{
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "path") == 0) { parse_string(card, &card->setup_path); }
        LOOP_OVER_KEYS_END
}

/* Update a RETRO card from a JSON file */
void retro_card_update(
    struct retro_card * card, struct roar_handler * handler, char * path)
{
        card->tea = jsmn_tea_create(path, JSMN_TEA_MODE_LOAD, handler);
        card->path = path;

        /* Loop over the fields. */
        LOOP_OVER_KEYS_BEGIN
        if (strcmp(tag, "comment") == 0) {
                jsmn_tea_next_string(card->tea, 0, NULL);
        } else if (strcmp(tag, "generator") == 0) {
                update_generator(card);
        } else if (strcmp(tag, "processor") == 0) {
                update_processor(card);
        } else if (strcmp(tag, "selector") == 0) {
                update_selector(card);
        } else if (strcmp(tag, "logger") == 0) {
                update_logger(card);
        } else if (strcmp(tag, "primary") == 0) {
                update_primary(card);
        } else if (strcmp(tag, "topography") == 0) {
                update_topography(card);
        } else if (strcmp(tag, "setup") == 0) {
                update_setup(card);
        }
        LOOP_OVER_KEYS_END

        jsmn_tea_destroy(&card->tea);
}
