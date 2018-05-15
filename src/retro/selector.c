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

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "danton.h"
#include "grand-tour.h"
#include "roar.h"

#include "retro/card.h"
#include "retro/constant.h"
#include "retro/selector.h"

static double select_none(const struct retro_selector * selector, double energy,
    const double * position, const double * direction)
{
        return -1.;
}

static double select_vertex(const struct retro_selector * selector,
    double energy, const double * position, const double * direction)
{
        /* First let us compute the decay length, assuming no energy loss */
        const double dl = energy * TAU_CTAU / TAU_MASS;

        /* Then let us compute the distance to the topography, propagating
         * backwards */
        double dg;
        double u[3] = { -direction[0], -direction[1], -direction[2] };
        int rc;
        if ((rc = gt_ground_distance(selector->topography, position, u,
                 selector->vertex_limit * dl, &dg)) != EXIT_SUCCESS) {
                ROAR_ERRWP_MESSAGE(selector->handler, &select_vertex, -1,
                    "turtle error", turtle_strerror(rc));
        }
        if (dg < 0.) return exp(-selector->vertex_limit);

        /* As selection weight, let us consider the probability that no decay
         * occured on the path to rocks */
        return exp(-dg / dl);
}

static double select_setup(const struct retro_selector * selector,
    double energy, const double * position, const double * direction);

/* Check if the given ray intersects the topography */
static int topography_intersect(const struct retro_selector * selector,
    const double r0[3], const double u[3], double smax)
{
#define RAY_STEP 200.
        struct turtle_datum * datum = danton_get_datum();
        const double ds = RAY_STEP;
        double r[3] = { r0[0], r0[1], r0[2] };
        double s;
        for (s = 0.; s < smax; s += ds) {
                /* Check the topography */
                double latitude, longitude, altitude;
                enum turtle_return rc;
                if ((rc = turtle_datum_geodetic(datum, r, &latitude, &longitude,
                         &altitude)) != TURTLE_RETURN_SUCCESS) {
                        ROAR_ERRWP_MESSAGE(selector->handler, &select_setup, -1,
                            "turtle error", turtle_strerror(rc));
                }
                double zg;
                if (selector->topography->flat) {
                        zg = 0.;
                } else {
                        if ((rc = turtle_datum_elevation(datum, latitude,
                                 longitude, &zg)) != TURTLE_RETURN_SUCCESS) {
                                ROAR_ERRWP_MESSAGE(selector->handler,
                                    &select_setup, -1, "turtle error",
                                    turtle_strerror(rc));
                        }
                }
                if (altitude < zg) return EXIT_FAILURE;

                /* Update the track position. */
                r[0] += ds * u[0];
                r[1] += ds * u[1];
                r[2] += ds * u[2];
        }
        return EXIT_SUCCESS;
}

static double select_setup(const struct retro_selector * selector,
    double energy, const double * position, const double * direction)
{
        /* Cone parameters */
        double tan_gamma = 0.;
        if (selector->setup_cone) {
                const double gamma = selector->setup_gamma(
                    selector, energy, position, direction);
                tan_gamma = tan(gamma * M_PI / 180.);
        }
        const double zcmin = 14E+03;
        const double zcmax = 165E+03 * energy / 1E+09 + 55E+03;

        /* Convert to ECEF */
        double r0[3], u[3];
        gt_to_ecef(selector->topography, position, 0, r0);
        gt_to_ecef(selector->topography, direction, 1, u);

        /* Check if the shower crashes into a mountain early, before xmax */
        if (selector->setup_xmax) {
                if (topography_intersect(selector, r0, u, zcmin + RAY_STEP) ==
                    EXIT_FAILURE)
                        return -1.;
        }
        if (!selector->setup_cone) return 100.;

        /* Check the antennas */
        const double r1[3] = { r0[0] + zcmin * u[0], r0[1] + zcmin * u[1],
                r0[2] + zcmin * u[2] };

        selector->setup_selection[0] = NULL;
        int i, triggers = 0;
        double * ra;
        for (i = 0, ra = selector->setup_data; i < selector->setup_n;
             i++, ra += selector->setup_size) {
                /* Check that the antenna is inside the cone */
                const double dx = ra[0] - r0[0];
                const double dy = ra[1] - r0[1];
                const double dz = ra[2] - r0[2];

                const double zp = u[0] * dx + u[1] * dy + u[2] * dz;
                if ((zp < zcmin) || (zp > zcmax)) continue;
                const double rp2 = dx * dx + dy * dy + dz * dz - zp * zp;
                const double rho = (zp - zcmin) * tan_gamma;
                if (rp2 > rho * rho) continue;

                /* Check for any shadowing by the topography. */
                if (selector->setup_shadowing) {
                        double n[3] = { dx - zcmin * u[0], dy - zcmin * u[1],
                                dz - zcmin * u[2] };
                        const double s2 =
                            n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
                        if (s2 <= FLT_EPSILON) continue;
                        const double smax = sqrt(s2);
                        if (smax <= RAY_STEP) continue;
                        const double d = 1. / smax;
                        n[0] *= d;
                        n[1] *= d;
                        n[2] *= d;
                        if (topography_intersect(selector, r1, n,
                                smax - RAY_STEP) == EXIT_FAILURE)
                                continue;
                }

                selector->setup_selection[triggers] = ra;
                triggers++;
        }
        selector->setup_selection[triggers] = NULL;

        return triggers;
}

static double cone_model_3deg(const struct retro_selector * selector,
    double energy, const double * position, const double * direction)
{
        return 3.;
}

static double cone_model_agressive(const struct retro_selector * selector,
    double energy, const double * position, const double * direction)
{
        return 0.47 * log(energy / 1E+08) + 0.9;
}

void selector_initialise(struct retro_selector * selector,
    const struct retro_card * card, struct roar_handler * handler,
    struct gt_topography * topography)
{
        selector->handler = handler;
        selector->topography = topography;

        /* Vertex selector */
        if (card->selector_vertex_limit > 0.) {
                selector->vertex_limit = card->selector_vertex_limit;
                selector->vertex = &select_vertex;
        } else {
                selector->vertex_limit = 0.;
                selector->vertex = &select_none;
        }

        /* Setup selector */
        int setup = (card->selector_setup_cone || card->selector_setup_xmax ||
                        card->selector_setup_shadowing) &&
            (card->setup_path != NULL);
        if (!setup) {
                selector->setup = NULL;
                return;
        }

        selector->setup_cone = card->selector_setup_cone;
        selector->setup_xmax = card->selector_setup_xmax;
        selector->setup_shadowing = card->selector_setup_shadowing;
        selector->setup = &select_setup;

        /* Parse the setup data */
        FILE * fid = fopen(card->setup_path, "rb");
        if (fid == NULL) {
                ROAR_ERRNO_MESSAGE(selector->handler, &selector_initialise,
                    errno, card->setup_path);
        }

        int n_antennas = -1, size = 0;
        char buffer[2048];
        const int nmax = sizeof(buffer) / sizeof(*buffer);
        while (!feof(fid)) {
                const int n = fread(buffer, sizeof(*buffer), nmax, fid);
                int i;
                char * c;
                for (i = 0, c = buffer; i < n; i++, c++) {
                        if (*c == '[')
                                n_antennas++;
                        else if (*c == ',')
                                size++;
                }
        }
        if (n_antennas <= 0) {
                ROAR_ERRNO_FORMAT(selector->handler, &selector_initialise,
                    EINVAL, "empty setup file `%s`", card->setup_path);
        }
        size = (size + 1 - n_antennas) / n_antennas + 1;

        selector->setup_n = n_antennas;
        selector->setup_size = size;
        selector->setup_data =
            malloc(size * n_antennas * sizeof(*selector->setup_data));
        if (selector->setup_data == NULL) {
                ROAR_ERRNO_MESSAGE(selector->handler, &selector_initialise,
                    ENOMEM, card->setup_path);
        }
        if (selector->setup_cone) {
                selector->setup_selection = malloc(
                    (n_antennas + 1) * sizeof(*selector->setup_selection));
                if (selector->setup_selection == NULL) {
                        ROAR_ERRNO_MESSAGE(selector->handler,
                            &selector_initialise, ENOMEM, card->setup_path);
                }
                selector->setup_gamma =
                    (selector->setup_cone == SETUP_CONE_MODEL_3DEG) ?
                    &cone_model_3deg :
                    &cone_model_agressive;
        } else {
                selector->setup_selection = NULL;
                selector->setup_gamma = NULL;
        }

        rewind(fid);
        int nread = nmax - 1;
        double * d = selector->setup_data;
        char * c = buffer;
        while (!feof(fid)) {
                const int nr = fread(c, sizeof(*c), nread, fid);
                if (nr <= 0) break;
                c[nr] = 0x0;
                c = buffer;
                char * last = c;
                for (; *c != 0x0; c++) {
                        if (isdigit(*c) || (*c == '-')) {
                                last = c;
                                *d++ = strtod(c, &c);
                                if (*c == '-') c++;
                                if (*c == 0x0) break;
                        }
                }
                d--;
                const int m = nmax - 1 - (last - buffer);
                memmove(buffer, last, m);
                c = buffer + m;
                nread = nmax - 1 - m;
        }

        /* Convert the antenna positions to ECEF */
        int i;
        for (i = 0, d = selector->setup_data; i < n_antennas; i++, d += size) {
                double ecef[3];
                gt_to_ecef(topography, d, 0, ecef);
                memcpy(d, ecef, sizeof(ecef));
        }

        fclose(fid);
}

void selector_finalise(struct retro_selector * selector)
{
        free(selector->setup_data);
        selector->setup_data = NULL;
        free(selector->setup_selection);
        selector->setup_selection = NULL;
}
