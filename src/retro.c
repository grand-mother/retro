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

/* Standard library includes. */
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The DANTON API */
#include "danton.h"
#include "danton/primary/powerlaw.h"

/* The GRAND_TOUR/TURTLE API */
#include "grand-tour.h"
#include "turtle.h"

/* The RETRO components */
#include "retro/card.h"
#include "retro/constant.h"
#include "retro/generator.h"
#include "retro/random.h"
#include "retro/selector.h"

/* The error handling with ROAR */
#include "roar.h"

/* Error handler. */
static struct roar_handler handler = { NULL, NULL, NULL, NULL };

/* Handle for the simulation context. */
static struct danton_context * context = NULL;

/* Storage for the GRAND topography */
static struct gt_topography topography;

/* Configuration card */
static struct retro_card card;

/* Data storage for selection functions */
static struct retro_selector selector;

/* Storage for the Monte-Carlo statistics */
struct {
        long long done;
        long long trials;
} statistics = { 0, 0 };

/* Storage for primary events */
struct primary_data {
        double weight;
        double energy;
        int generation;
        int medium;
        double local[3];
        double lla[3];
};

static double primary_weight = 0.;
static int primary_n = 0;
static int primary_trials = 0;
static struct primary_data * pdata = NULL;

/* Finalise and exit to the OS. */
static int gracefully_exit(int rc)
{
        free(pdata);
        pdata = NULL;
        retro_card_destroy(&card);
        selector_finalise(&selector);
        if (context != NULL) {
                int i;
                for (i = 0; i < DANTON_PARTICLE_N_NU; i++)
                        danton_destroy((void **)&context->primary[i]);
                danton_destroy((void **)&context->sampler);
                danton_context_destroy(&context);
        }
        danton_finalise();
        exit(rc);
}

/* Show help and exit. */
static void exit_with_help(int rc)
{
        /* This is the help string formated as it should output on the screen.
         * Therefore it follows a specific formating rule.
         */
        // clang-format off
        fprintf(stderr,
"Usage: retro [DATACARD_1.JSON] ... [DATACARD_N.JSON]\n"
"Simulate decaying tau neutrinos for a radio-detector\n"
"\n"
"Data card(s):\n"
"Syntax and examples available from https://github.com/grand-mother/retro.\n"
"\n"
"Exit status:\n"
" %d  if OK,\n"
" %d  if an error occurred.\n"
"\n"
"License: GNU LGPLv3\n"
"Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC.\n"
"Author: Valentin NIESS (niess@in2p3.fr)\n"
"\n", EXIT_SUCCESS, EXIT_FAILURE);
        // clang-format on
        gracefully_exit(rc);
}

/* Post-error callback. */
int handle_post_error(
    struct roar_handler * handler, roar_function_t * referent, int code)
{
        /* Finalise and exit to the OS. */
        return gracefully_exit(EXIT_FAILURE);
}

/* Callback for recording an event. */
static int record_event(struct danton_context * context,
    struct danton_recorder * recorder, const struct danton_event * event)
{
        primary_trials = event->id + 1;
        struct primary_data * p = pdata + primary_n++;
        const double etau = event->final->energy;
        p->weight = event->weight * primary_weight /
            sqrt((etau - TAU_MASS) * (etau + TAU_MASS));
        p->energy = event->primary->energy;
        p->generation = event->generation;

        memcpy(p->local, event->vertex->position, sizeof(p->local));
        int rc;
        if ((rc = gt_to_lla(&topography, event->vertex->position, p->lla)) !=
            EXIT_SUCCESS) {
                ROAR_ERRWP_MESSAGE(&handler, &record_event, -1, "turtle error",
                    turtle_strerror(rc));
        }

        int above;
        if ((rc = gt_ground_above(&topography, event->vertex->position,
                 &above)) != EXIT_SUCCESS) {
                if (rc == TURTLE_RETURN_PATH_ERROR) {
                        above = p->lla[2] >= 0.;
                } else {
                        ROAR_ERRWP_MESSAGE(&handler, &record_event, -1,
                            "turtle error", turtle_strerror(rc));
                }
        }
        p->medium = above;

        return EXIT_SUCCESS;
}

int main(int narg, char * argv[]);

static void raise_danton(void)
{
        ROAR_ERRWP_MESSAGE(
            &handler, &main, -1, "danton error", danton_error_pop(NULL));
}

static void raise_turtle(int rc)
{
        ROAR_ERRWP_MESSAGE(
            &handler, &main, -1, "turtle error", turtle_strerror(rc));
}

/* The main executable */
int main(int narg, char * argv[])
{
        if (narg <= 1) exit_with_help(EXIT_SUCCESS);

        /* Configure the error handler */
        errno = 0;
        handler.stream = stderr;
        handler.post = &handle_post_error;

        /* Parse the configuration card(s) */
        retro_card_initialise(&card);
        while (*(++argv) != NULL) retro_card_update(&card, &handler, *argv);

        if ((card.primary_requested > 0) || (card.primary_events > 0)) {
                int n = (card.primary_requested > 0) ? card.primary_requested :
                                                       card.primary_events;
                pdata = malloc(n * sizeof(*pdata));
                if (pdata == NULL) { ROAR_ERRNO(&handler, &main, ENOMEM); }
        }

        /* Reset the output file */
        FILE * fd = fopen(card.logger_path, "w+");
        if (fd == NULL) {
                ROAR_ERRNO_MESSAGE(&handler, &main, 0, card.logger_path);
        }
        fclose(fd);

        /* Location of the previous event in the output file. */
        long previous[2] = { -1, 0 };

        /* Initialise DANTON */
        if (danton_initialise(NULL, NULL, NULL, NULL, NULL) != EXIT_SUCCESS)
            raise_danton();

        /* Configure the topography */
        int sea = 0;
        if (danton_earth_model("WGS84", card.topography_path,
                card.topography_stack_size, "Rock", card.topography_density,
                &sea) != EXIT_SUCCESS)
                raise_danton();

        struct turtle_datum * datum = danton_get_datum();
        gt_initialise(&topography, card.topography_latitude,
            card.topography_longitude, card.topography_path,
            card.topography_stack_size, NULL, datum);

        /* Create a DANTON simulation context */
        context = danton_context_create();
        if (context == NULL) raise_danton();

        /* Initialise the random engine. */
        random_initialise(context);

        /* Configure the simulation */
        context->mode = DANTON_MODE_BACKWARD;
        context->longitudinal = card.primary_longitudinal;
        context->decay = 0;

        struct danton_primary * primary =
            (struct danton_primary *)danton_powerlaw_create(
                card.generator_energy[0], card.generator_energy[1] * 1E+03, -2.,
                1.);
        if (primary == NULL) raise_danton();
        context->primary[5] = primary;
        primary_weight = TAU_MASS / TAU_CTAU *
            (1. / card.generator_energy[0] - 1E-03 / card.generator_energy[1]);

        struct danton_sampler * sampler = danton_sampler_create();
        if (sampler == NULL) raise_danton();
        context->sampler = sampler;
        sampler->weight[7] = 1.;

        struct danton_recorder recorder = { &record_event, NULL };
        context->recorder = &recorder;

        struct retro_generator generator;
        generator_initialise(&generator, &card, &handler, &topography);
        selector_initialise(&selector, &card, &handler, &topography);

        /* Container for the tau data */
        struct {
                double energy;
                double position[3];
                double direction[3];
        } tau_at_decay;

        /* Monte-Carlo main loop */
        long long trials = 0;
        for (;;) {
                /* Check the termination conditions */
                if ((card.processor_requested > 0) &&
                    (statistics.done >= card.processor_requested))
                        break;
                if ((card.processor_trials > 0) &&
                    (statistics.trials >= card.processor_trials))
                        break;

                /* Generate a tentative decay vertex */
                statistics.trials += 1;
                trials += 1;
                double w;
                if ((w = generator.position(
                         &generator, tau_at_decay.position)) <= 0.)
                        continue;
                double angle[2];
                w *= generator.direction(&generator, angle);
                int rc;
                if ((rc = gt_from_angular(&topography, tau_at_decay.position,
                         angle, tau_at_decay.direction)) != EXIT_SUCCESS)
                        raise_turtle(rc);

                w *= generator.energy(&generator, &tau_at_decay.energy);

                /* Check if the generated direction is relevant considering the
                 * generated position and its energy
                 */
                const double psel =
                    selector.vertex(&selector, tau_at_decay.energy,
                        tau_at_decay.position, tau_at_decay.direction);
                if ((psel <= 0.) || (random_uniform01() > psel)) continue;
                w /= psel;

                /* Generate a valid tau decay, i.e. not a muonic decay */
                generator_decay_initialise(&generator, 15, tau_at_decay.energy,
                    tau_at_decay.direction);
                double shower_energy;
                w *= generator.decay(&generator, &shower_energy);
                if (shower_energy < card.generator_energy[0]) continue;

                /* Preselect antennas that might detect the radio signal from
                 * the shower */
                int triggers = 0;
                if (selector.setup != NULL) {
                        if ((triggers = selector.setup(&selector, shower_energy,
                                 tau_at_decay.position,
                                 tau_at_decay.direction)) < 4.)
                                continue;
                }

                /* Sample the primary flux */
                double lla[3];
                gt_to_lla(&topography, tau_at_decay.position, lla);

                primary_trials = 0;
                if ((card.primary_requested > 0) || (card.primary_events > 0)) {
                        sampler->latitude = lla[0];
                        sampler->longitude = lla[1];
                        sampler->altitude[0] = sampler->altitude[1] = lla[2];
                        sampler->azimuth[0] = sampler->azimuth[1] = -angle[1];
                        sampler->elevation[0] = sampler->elevation[1] =
                            90. - angle[0];
                        sampler->energy[0] = sampler->energy[1] =
                            tau_at_decay.energy;

                        if (danton_sampler_update(context->sampler) !=
                            EXIT_SUCCESS) {
                                ROAR_ERRWP_MESSAGE(&handler, &main, -1,
                                    "danton error", danton_error_pop(NULL));
                        }

                        primary_n = 0;
                        if (danton_run(context, card.primary_events,
                                card.primary_requested) != EXIT_SUCCESS)
                                ROAR_ERRWP_MESSAGE(&handler, &main, -1,
                                    "danton error", danton_error_pop(context));

                        if (primary_n <= 0) continue;
                        if ((card.primary_requested <= 0) ||
                            (primary_n < card.primary_requested)) {
                                primary_trials = card.primary_events;
                        }
                }

                /* Build the event tag */
                if (angle[1] < 0.) angle[1] += 360.;
                char subtag[16], tag[4096];
                sprintf(subtag, "%.0e", tau_at_decay.energy * 1E+09);
                memmove(subtag + 2, subtag + 3, 3 * sizeof(char));
                sprintf(tag,
                    "E.%s_Z.%.0lf_A.%.0lf_La.%.0lf_Lo.%.0lf_H.%.0lf_D.%lld",
                    subtag, angle[0], angle[1], lla[0], lla[1], lla[2],
                    generator.decay_state);

                /* Dump the event */
                FILE * fd = fopen(card.logger_path, "a+");
                if (fd == NULL) {
                        ROAR_ERRNO_MESSAGE(
                            &handler, &main, 0, card.logger_path);
                }
                fprintf(fd, "{\"tag\" : \"%s\", ", tag);
                fprintf(fd,
                    "\"tau_at_decay\" : [%.5E, %.5E, [%.3lf, %.3lf, "
                    "%.3lf], [%.5E, %.5E, %.5E], [%.8lf, %.8lf, "
                    "%.3lf], [%.3lf, %.3lf]], ",
                    w, tau_at_decay.energy, tau_at_decay.position[0],
                    tau_at_decay.position[1], tau_at_decay.position[2],
                    tau_at_decay.direction[0], tau_at_decay.direction[1],
                    tau_at_decay.direction[2], lla[0], lla[1], lla[2], angle[0],
                    angle[1]);
                fputs("\"decay\" : [", fd);
                struct retro_product * dp;
                int i;
                for (i = 0, dp = generator.decay_product;
                     i < generator.decay_product_n; i++, dp++) {
                        if (i > 0) fputs(", ", fd);
                        fprintf(fd, "[%d, [%.5E, %.5E, %.5E]]", dp->pid,
                            dp->momentum[0], dp->momentum[1], dp->momentum[2]);
                }
                fputs("], \"primaries\" : [", fd);
                struct primary_data * p;
                for (i = 0, p = pdata; i < primary_n; i++, p++) {
                        p->weight *= w * p->energy * p->energy;
                        if (i > 0) fputs(", ", fd);
                        fprintf(fd,
                            "[%.5E, %.5E, %d, %d, [%.3lf, %.3lf, "
                            "%.3lf], [%.3lf, %.3lf, %.3lf]]",
                            p->weight, p->energy, p->generation, p->medium,
                            p->local[0], p->local[1], p->local[2], p->lla[0],
                            p->lla[1], p->lla[2]);
                }
                fprintf(fd, "], \"statistics\" : [%lld, %d], ", trials,
                    primary_trials);
                fputs("\"antennas\" : [", fd);
                double ** s = selector.setup_selection;
                double * data = selector.setup_data;
                for (i = 0; i < triggers; i++) {
                        double *antenna, *local;
                        double tmp[3];
                        if (selector.array == SETUP_ARRAY_MODEL_FILE) {
                                antenna = *s++;
                                gt_from_ecef(&topography, antenna, 0, tmp);
                                local = tmp;
                        } else {
                                local = antenna = data;
                                data += selector.setup_size;
                        }
                        if (i > 0) fputs(", ", fd);
                        fprintf(fd, "[%.3lf, %.3lf, %.3lf", local[0], local[1],
                            local[2]);
                        int j;
                        for (j = 3; j < selector.setup_size; j++)
                                fprintf(fd, ", %.3lf", antenna[j]);
                        fputs("]", fd);
                }
                fprintf(fd, "], \"origin\" : [%.8lf, %.8lf], ",
                    card.topography_latitude, card.topography_longitude);
                fprintf(fd, "\"previous\" : %ld}\n", previous[0]);
                previous[0] = previous[1];
                previous[1] = ftell(fd);
                fclose(fd);

                trials = 0;
                statistics.done += 1;
        }

        gracefully_exit(EXIT_SUCCESS);
}
