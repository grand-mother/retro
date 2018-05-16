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
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "alouette.h"
#include "grand-tour.h"
#include "roar.h"

#include "retro/card.h"
#include "retro/constant.h"
#include "retro/generator.h"
#include "retro/random.h"

static double generate_position(
    struct retro_generator * generator, double * position)
{
        int i;
        for (i = 0; i < 3; i++) {
                position[i] = generator->position_parameter[i][0] +
                    generator->position_parameter[i][1] * random_uniform01();
        }

        if (generator->position_mode == RETRO_GENERATOR_MODE_GEODETIC) {
                int rc;
                double lla[3];
                if ((rc = gt_to_lla(generator->topography, position, lla)) !=
                    EXIT_SUCCESS) {
                        ROAR_ERRWP_MESSAGE(generator->handler,
                            &generate_position, -1, "turtle error",
                            turtle_strerror(rc));
                }
                if ((lla[0] < generator->position_parameter[3][0]) ||
                    (lla[0] > generator->position_parameter[3][1]) ||
                    (lla[1] < generator->position_parameter[4][0]) ||
                    (lla[1] > generator->position_parameter[4][1]))
                        return -1.;
        }

        int rc;
        double zg;
        if ((rc = gt_ground_altitude(
                 generator->topography, position, 0, &zg)) != EXIT_SUCCESS) {
                ROAR_ERRWP_MESSAGE(generator->handler, &generate_position, -1,
                    "turtle error", turtle_strerror(rc));
        }
        position[2] += zg;

        return generator->position_weight;
}

static double generate_direction_uniform(
    struct retro_generator * generator, double * angle)
{
        const double c = generator->direction_parameter[0] +
            generator->direction_parameter[1] * random_uniform01();
        angle[0] = acos(c) * 180. / M_PI;
        angle[1] = -180. + 360. * random_uniform01();
        return generator->direction_weight;
}

static double generate_direction_linear(
    struct retro_generator * generator, double * angle)
{
        double c;
        for (;;) {
                c = generator->direction_parameter[0] +
                    generator->direction_parameter[1] * random_uniform01();
                if (c > 0.) break;
        }
        c = sqrt(c);
        angle[0] = acos(c * generator->direction_parameter[2]) * 180. / M_PI;
        angle[1] = -180. + 360. * random_uniform01();
        return generator->direction_weight / c;
}

static double generate_energy_uniform(
    struct retro_generator * generator, double * energy)
{
        *energy = generator->energy_parameter[0] +
            generator->energy_parameter[1] * random_uniform01();
        return generator->energy_weight;
}

static double generate_energy_1_over_E(
    struct retro_generator * generator, double * energy)
{
        const double e = generator->energy_parameter[0] *
            exp(generator->energy_parameter[1] * random_uniform01());
        *energy = e;
        return generator->energy_weight * e;
}

static double generate_energy_1_over_E2(
    struct retro_generator * generator, double * energy)
{
        const double e = generator->energy_parameter[0] /
            (1. - generator->energy_parameter[1] * random_uniform01());
        *energy = e;
        return generator->energy_weight * e * e;
}

void generator_decay_initialise(struct retro_generator * generator, int pid,
    double energy, double * direction)
{
        static int initialised = 0;
        if (!initialised) {
                unsigned int state[3] = { 0, 0, 0 };
                state[0] = (unsigned int)(random_uniform01() * 900000000);
                int rc;
                if ((rc = alouette_initialise(1, state)) != EXIT_SUCCESS) {
                        ROAR_ERRWP_MESSAGE(generator->handler,
                            &generator_decay_initialise, -1, "alouette error",
                            alouette_strerror(rc));
                }
                initialised = 1;
        }

        generator->decay_pid = pid;
        const double p = sqrt((energy + TAU_MASS) * (energy - TAU_MASS));
        int i;
        for (i = 0; i < 3; i++) {
                generator->decay_momentum[i] = direction[i] * p;
                generator->decay_polarisation[i] = direction[i];
        }
}

static double generate_decay(
    struct retro_generator * generator, double * shower_energy)
{
        unsigned int state[3];
        double esh = -1.;
        while (esh <= 0.) {
                alouette_random_state(state);
                int rc;
                if ((rc = alouette_decay(generator->decay_pid,
                         generator->decay_momentum,
                         generator->decay_polarisation)) != EXIT_SUCCESS) {
                        ROAR_ERRWP_MESSAGE(generator->handler, &generate_decay,
                            -1, "alouette error", alouette_strerror(rc));
                }

                int pid, i = 0;
                double p[3];
                while (alouette_product(&pid, p) == EXIT_SUCCESS) {
                        generator->decay_product[i].pid = pid;
                        memcpy(
                            generator->decay_product[i].momentum, p, sizeof(p));
                        i++;
                        if (i >= RETRO_MAX_DECAY_PRODUCTS) {
                                ROAR_ERRWP_MESSAGE(generator->handler,
                                    &generate_decay, -1, "retro error",
                                    "maximum number of decay products "
                                    "exceeded");
                        }

                        int aid = abs(pid);
                        if (((aid >= 12) && (aid <= 14)) || (aid == 16))
                                continue;
                        esh += sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                }
                generator->decay_product_n = i;
        }

        generator->decay_state =
            (state[2] * 1000000000 + state[1]) * 1000000000 + state[0];
        *shower_energy = esh;
        return 1. - TAU_BR_MU;
}

static void update_bounding_box(struct retro_generator * generator,
    double latitude, double longitude, double bounding_box[2][2])
{
        double lla[3] = { latitude, longitude, 0. }, local[3];
        int rc;
        if ((rc = gt_from_lla(generator->topography, lla, local)) !=
            EXIT_SUCCESS) {
                ROAR_ERRWP_MESSAGE(generator->handler, &generator_initialise,
                    -1, "turtle error", turtle_strerror(rc));
        }

        int i;
        for (i = 0; i < 2; i++) {
                if (local[i] < bounding_box[i][0])
                        bounding_box[i][0] = local[i];
                else if (local[i] > bounding_box[i][1])
                        bounding_box[i][1] = local[i];
        }
}

static double angle2cos(double angle)
{
        const double c = cos(angle * M_PI / 180.);
        return (fabs(c) < FLT_EPSILON) ? 0. : c;
}

void generator_initialise(struct retro_generator * generator,
    struct retro_card * card, struct roar_handler * handler,
    struct gt_topography * topography)
{
        generator->handler = handler;
        generator->topography = topography;

        generator->position_mode = card->generator_position_mode;
        generator->position_weight = 1;
        if (generator->position_mode == RETRO_GENERATOR_MODE_GEODETIC) {
                /* Compute the bounding box parameters */
                double bounding_box[2][2] = { { DBL_MAX, -DBL_MAX },
                        { DBL_MAX, -DBL_MAX } };
                update_bounding_box(generator,
                    card->topography_latitude + card->generator_position[0][0],
                    card->topography_longitude + card->generator_position[1][0],
                    bounding_box);
                update_bounding_box(generator,
                    card->topography_latitude + card->generator_position[0][1],
                    card->topography_longitude + card->generator_position[1][0],
                    bounding_box);
                update_bounding_box(generator,
                    card->topography_latitude + card->generator_position[0][1],
                    card->topography_longitude + card->generator_position[1][1],
                    bounding_box);
                update_bounding_box(generator,
                    card->topography_latitude + card->generator_position[0][0],
                    card->topography_longitude + card->generator_position[1][1],
                    bounding_box);
                update_bounding_box(generator,
                    card->topography_latitude + card->generator_position[0][0],
                    card->topography_longitude, bounding_box);
                update_bounding_box(generator,
                    card->topography_latitude + card->generator_position[0][1],
                    card->topography_longitude, bounding_box);
                update_bounding_box(generator, card->topography_latitude,
                    card->topography_longitude + card->generator_position[1][0],
                    bounding_box);
                update_bounding_box(generator, card->topography_latitude,
                    card->topography_longitude + card->generator_position[1][1],
                    bounding_box);

                /* Set the geodetic box */
                int j;
                for (j = 0; j < 2; j++) {
                        generator->position_parameter[3][j] =
                            card->topography_latitude +
                            card->generator_position[0][j];
                        generator->position_parameter[4][j] =
                            card->topography_longitude +
                            card->generator_position[1][j];
                }

                /* Copy the local box */
                int i;
                for (i = 0; i < 2; i++) {
                        for (j = 0; j < 2; j++) {
                                card->generator_position[i][j] =
                                    bounding_box[i][j];
                        }
                }
        }

        int i;
        for (i = 0; i < 3; i++) {
                const double d = card->generator_position[i][1] -
                    card->generator_position[i][0];
                generator->position_parameter[i][0] =
                    card->generator_position[i][0];
                generator->position_parameter[i][1] = d;
                generator->position_weight *= d;
        }
        generator->position = &generate_position;

        const double c0 = angle2cos(card->generator_theta[1]);
        const double c1 = angle2cos(card->generator_theta[0]);
        if (card->generator_theta_mode == RETRO_GENERATOR_MODE_UNIFORM) {
                const double dc = c1 - c0;
                generator->direction_parameter[0] = c0;
                generator->direction_parameter[1] = dc;
                generator->direction_weight = 2. * M_PI * dc;
                generator->direction = &generate_direction_uniform;
        } else if (card->generator_theta_mode == RETRO_GENERATOR_MODE_LINEAR) {
                if (c1 * c0 < 0.) {
                        ROAR_ERRNO_MESSAGE(handler, &generator_initialise,
                            EINVAL,
                            "invalid theta values for linear generator");
                }
                const double c02 = c0 * c0;
                generator->direction_parameter[0] = c02;
                const double dc2 = c1 * c1 - c02;
                generator->direction_parameter[1] = dc2;
                generator->direction_parameter[2] = (c1 < 0.) ? -1. : 1.;
                generator->direction_weight = fabs(dc2) * M_PI;
                generator->direction = &generate_direction_linear;
        } else {
                ROAR_ERRNO_MESSAGE(handler, &generator_initialise, EINVAL,
                    "invalid generator mode for the theta angle");
        }

        if (card->generator_energy_mode == RETRO_GENERATOR_MODE_UNIFORM) {
                const double de =
                    card->generator_energy[1] - card->generator_energy[0];
                generator->energy_parameter[0] = card->generator_energy[0];
                generator->energy_parameter[1] = de;
                generator->energy_weight = de;
                generator->energy = &generate_energy_uniform;
        } else if (card->generator_energy_mode ==
            RETRO_GENERATOR_MODE_1_OVER_E) {
                const double lne =
                    log(card->generator_energy[1] / card->generator_energy[0]);
                generator->energy_parameter[0] = card->generator_energy[0];
                generator->energy_parameter[1] = lne;
                generator->energy_weight = lne;
                generator->energy = &generate_energy_1_over_E;
        } else if (card->generator_energy_mode ==
            RETRO_GENERATOR_MODE_1_OVER_E2) {
                const double r =
                    1. - card->generator_energy[0] / card->generator_energy[1];
                generator->energy_parameter[0] = card->generator_energy[0];
                generator->energy_parameter[1] = r;
                generator->energy_weight = r / card->generator_energy[0];
                generator->energy = &generate_energy_1_over_E2;
        } else {
                ROAR_ERRNO_MESSAGE(handler, &generator_initialise, EINVAL,
                    "invalid generator mode for the energy");
        }

        generator->decay = &generate_decay;
}
