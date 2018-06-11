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

static double select_none(struct retro_selector * selector, double energy,
    const double * position, const double * direction)
{
        return -1.;
}

static double select_vertex(struct retro_selector * selector, double energy,
    const double * position, const double * direction)
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

static double select_setup(struct retro_selector * selector, double energy,
    const double * position, const double * direction);

static double topography_altitude(struct retro_selector * selector,
    struct turtle_datum * datum, double latitude, double longitude)
{
        if (selector->topography->flat) return 0.;

        double zg;
        const enum turtle_return rc =
            turtle_datum_elevation(datum, latitude, longitude, &zg);
        if (rc == TURTLE_RETURN_PATH_ERROR) {
                return -1.;
        } else if (rc != TURTLE_RETURN_SUCCESS) {
                ROAR_ERRWP_MESSAGE(selector->handler, &select_setup, -1,
                    "turtle error", turtle_strerror(rc));
        }
        return zg;
}

/* Check if the given ray intersects the topography */
static int topography_intersect(struct retro_selector * selector,
    const double r0[3], const double u[3], double smin, double smax,
    double * lla)
{
#define RAY_STEP 10.
        struct turtle_datum * datum = danton_get_datum();
        double r[3] = { r0[0] + smin * u[0], r0[1] + smin * u[1],
                r0[2] + smin * u[2] };
        double s = smin;
        for (;;) {
                /* Get the geodetic coordinates */
                double latitude, longitude, altitude;
                enum turtle_return rc;
                if ((rc = turtle_datum_geodetic(datum, r, &latitude, &longitude,
                         &altitude)) != TURTLE_RETURN_SUCCESS) {
                        ROAR_ERRWP_MESSAGE(selector->handler, &select_setup, -1,
                            "turtle error", turtle_strerror(rc));
                }

                /* Check the altitude */
                double zg;
                if (altitude <= 0.) {
                        if (lla != NULL) {
                                lla[0] = latitude;
                                lla[1] = longitude;
                        }
                        return EXIT_SUCCESS;
                } else if (altitude <= 8850.) {
                        /* Check the topography */
                        zg = topography_altitude(
                            selector, datum, latitude, longitude);
                        if (altitude < zg) {
                                if (lla != NULL) {
                                        lla[0] = latitude;
                                        lla[1] = longitude;
                                }
                                return EXIT_SUCCESS;
                        }
                } else {
                        zg = 8850.;
                }

                /* Update the track position. */
                double ds = 0.5 * fabs(altitude - zg);
                if (ds < RAY_STEP) ds = RAY_STEP;
                s += ds;
                if (s > smax) break;

                r[0] += ds * u[0];
                r[1] += ds * u[1];
                r[2] += ds * u[2];
        }
        return EXIT_FAILURE;
}

static int check_antenna(struct retro_selector * selector, const double * ra,
    const double * u, const double * r0, const double * r1, double zcmin,
    double zcmax, double tan_gamma)
{
        /* Check that the antenna is inside the cone */
        const double dx = ra[0] - r0[0];
        const double dy = ra[1] - r0[1];
        const double dz = ra[2] - r0[2];

        const double zp = u[0] * dx + u[1] * dy + u[2] * dz;
        if ((zp < zcmin) || (zp > zcmax)) return EXIT_FAILURE;
        const double rp2 = dx * dx + dy * dy + dz * dz - zp * zp;
        const double rho = zp * tan_gamma;
        if (rp2 > rho * rho) return EXIT_FAILURE;

        /* Check for any shadowing by the topography. */
        if (selector->setup_shadowing) {
                double n[3] = { dx, dy, dz };
                const double s2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
                if (s2 <= FLT_EPSILON) return EXIT_FAILURE;
                const double smax = sqrt(s2);
                if (smax <= RAY_STEP) return EXIT_FAILURE;
                const double d = 1. / smax;
                n[0] *= d;
                n[1] *= d;
                n[2] *= d;
#define SAFETY_DIST 50.
                if (topography_intersect(selector, r0, n, 0.,
                        smax - SAFETY_DIST, NULL) == EXIT_SUCCESS)
                        return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
}

static double select_setup(struct retro_selector * selector, double energy,
    const double * position, const double * direction)
{
        /* Cone parameters */
        double gamma = 0.;
        if (selector->setup_cone) {
                gamma = selector->setup_gamma(
                            selector, energy, position, direction) *
                    M_PI / 180.;
        }
        const double zcmin = 14E+03;
        const double zcmax = 165E+03 * energy / 1E+09 + 55E+03;

        /* Convert to ECEF */
        double r0[3], u[3];
        gt_to_ecef(selector->topography, position, 0, r0);
        gt_to_ecef(selector->topography, direction, 1, u);

        /* Check if the shower crashes into a mountain early, before xmax */
        if (selector->setup_xmax) {
                if (topography_intersect(selector, r0, u, 0, zcmin + RAY_STEP,
                        NULL) == EXIT_SUCCESS)
                        return -1.;
        }
        if (!selector->setup_cone) return 100.;

        /* Select antennas */
        const double r1[3] = { r0[0] + zcmin * u[0], r0[1] + zcmin * u[1],
                r0[2] + zcmin * u[2] };
        const double tan_gamma = tan(gamma);
        int triggers = 0;

        if (selector->array == SETUP_ARRAY_MODEL_WORLD_WIDE) {
                /* Generate the antenna array. 1st let us compute a bounding
                 * box for the cone footprint.
                 */
                const double n0nrm = 1. / sqrt(u[0] * u[0] + u[1] * u[1]);
                const double n0[3] = { -u[1] * n0nrm, u[0] * n0nrm, 0. };
                double n1[3] = { n0[1] * u[2], -n0[0] * u[2],
                        n0[0] * u[1] - n0[1] * u[0] };
                const double n1nrm =
                    1. / sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
                int i;
                for (i = 0; i < 3; i++) n1[i] *= n1nrm;

                double bounding_box[2][2] = { { DBL_MAX, -DBL_MAX },
                        { DBL_MAX, -DBL_MAX } };
                int k;
                for (k = 10; k >= 0; k--) {
                        const double cg = cos(gamma * k * 0.1);
                        const double sg = sin(gamma * k * 0.1);

                        const int n = 60;
                        const double dphi = 2. * M_PI / n;
                        double phi = 0.;
                        for (i = 0; i < n; i++, phi += dphi) {
                                const double c0 = sg * cos(phi);
                                const double c1 = sg * sin(phi);
                                const double v[3] = { c0 * n0[0] + c1 * n1[0] +
                                            cg * u[0],
                                        c0 * n0[1] + c1 * n1[1] + cg * u[1],
                                        c0 * n0[2] + c1 * n1[2] + cg * u[2] };
                                double lla[2];
                                if (topography_intersect(selector, r0, v, zcmin,
                                        zcmax, lla) == EXIT_FAILURE)
                                        continue;
                                int j;
                                for (j = 0; j < 2; j++) {
                                        if (lla[j] < bounding_box[j][0])
                                                bounding_box[j][0] = lla[j];
                                        if (lla[j] > bounding_box[j][1])
                                                bounding_box[j][1] = lla[j];
                                }
                        }
                        if (k == 10) {
                                /* If the cone envelop doesn't intersect the
                                 * topography it is useless to look for inner
                                 * intersections.
                                 */
                                if ((bounding_box[0][0] == DBL_MAX) ||
                                    (bounding_box[1][0] == DBL_MAX))
                                        return 0;
                        }
                }
                if ((bounding_box[0][0] == DBL_MAX) ||
                    (bounding_box[1][0] == DBL_MAX))
                        return 0.;

                /* Second, let us loop over antennas in the bounding box */
                struct turtle_datum * datum = danton_get_datum();
                const double rE = 6367444.65;
                const double dlat =
                    selector->setup_ww_step * 180. / (rE * M_PI);
                bounding_box[0][0] = ((int)(bounding_box[0][0] / dlat)) * dlat;
                bounding_box[0][1] =
                    ((int)(bounding_box[0][1] / dlat) + 1) * dlat;
                double latitude;
                for (latitude = bounding_box[0][0];
                     latitude <= bounding_box[0][1]; latitude += dlat) {
                        const double dlong = selector->setup_ww_step * 180. /
                            (rE * sin(latitude * M_PI / 180.) * M_PI);
                        const double longmin =
                            ((int)(bounding_box[1][0] / dlong)) * dlong;
                        const double longmax =
                            ((int)(bounding_box[1][1] / dlong) + 1) * dlong;
                        double longitude;
                        for (longitude = longmin; longitude <= longmax;
                             longitude += dlong) {
                                double altitude = topography_altitude(
                                    selector, datum, latitude, longitude);
                                if (altitude < 0.) continue;
                                altitude += selector->setup_ww_height;
                                double ra[3];
                                enum turtle_return rc;
                                if ((rc = turtle_datum_ecef(datum, latitude,
                                         longitude, altitude, ra)) !=
                                    TURTLE_RETURN_SUCCESS) {
                                        ROAR_ERRWP_MESSAGE(selector->handler,
                                            &select_setup, -1, "turtle error",
                                            turtle_strerror(rc));
                                }
                                if (check_antenna(selector, ra, u, r0, r1,
                                        zcmin, zcmax,
                                        tan_gamma) == EXIT_FAILURE)
                                        continue;

                                /* Add the antenna to the candidates */
                                triggers++;
                                if (triggers > selector->setup_n) {
/* Re-allocate memory for the candidate
 * antennas
 */
#ifndef PAGESIZE
#define PAGESIZE 4096
#endif
                                        const int size = selector->setup_size *
                                            sizeof(*selector->setup_data);
                                        const int n_pages = selector->setup_n /
                                            (PAGESIZE / size);
                                        void * tmp =
                                            realloc(selector->setup_data,
                                                (n_pages + 1) * PAGESIZE);
                                        if (tmp == NULL) {
                                                ROAR_ERRNO(selector->handler,
                                                    &select_setup, ENOMEM);
                                        }
                                        selector->setup_data = tmp;
                                        selector->setup_n += PAGESIZE / size;
                                }
                                double antenna[5] = { 0 }, normal[3] = { 0 };
                                gt_from_ecef(
                                    selector->topography, ra, 0, antenna);
                                gt_ground_normal(selector->topography, antenna,
                                    0, 200., normal, antenna + 3);
                                memcpy(selector->setup_data +
                                        selector->setup_size * (triggers - 1),
                                    antenna, sizeof(antenna));
                        }
                }
        } else {
                /* Check the antennas */
                selector->setup_selection[0] = NULL;
                int i;
                double * ra;
                for (i = 0, ra = selector->setup_data; i < selector->setup_n;
                     i++, ra += selector->setup_size) {
                        if (check_antenna(selector, ra, u, r0, r1, zcmin, zcmax,
                                tan_gamma) == EXIT_FAILURE)
                                continue;
                        selector->setup_selection[triggers] = ra;
                        triggers++;
                }
                selector->setup_selection[triggers] = NULL;
        }

        return triggers;
}

static double cone_model_3deg(struct retro_selector * selector, double energy,
    const double * position, const double * direction)
{
        return 3.;
}

static double cone_model_agressive(struct retro_selector * selector,
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
        selector->array = SETUP_ARRAY_MODEL_NONE;
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

        if (selector->setup_cone) {
                selector->setup_gamma =
                    (selector->setup_cone == SETUP_CONE_MODEL_3DEG) ?
                    &cone_model_3deg :
                    &cone_model_agressive;
        } else {
                selector->setup_gamma = NULL;
        }

        if ((strlen(card->setup_path) >= 5) &&
            (strncmp(card->setup_path, "ww://", 5) == 0)) {
                /* Use a world wide array */
                selector->array = SETUP_ARRAY_MODEL_WORLD_WIDE;
                if (sscanf(card->setup_path, "ww://%lf/%lf",
                        &selector->setup_ww_step,
                        &selector->setup_ww_height) != 2) {
                        ROAR_ERRNO_MESSAGE(selector->handler,
                            &selector_initialise, EINVAL, card->setup_path);
                }
                selector->setup_n = 0;
                selector->setup_size = 5;
                selector->setup_data = NULL;
                selector->setup_selection = NULL;
                return;
        }

        /* Parse the array data from a file */
        selector->array = SETUP_ARRAY_MODEL_FILE;
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
        } else {
                selector->setup_selection = NULL;
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
