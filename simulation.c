/* Functions pertinent to the outer simulation steps */

#include <math.h>

#include "lbm.h"

void timestep(const param_t params, const accel_area_t accel_area,
    speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    accelerate_flow(params,accel_area,cells,obstacles);
    propagate(params,cells,tmp_cells);
    rebound(params,cells,tmp_cells,obstacles);
    collision(params,cells,tmp_cells,obstacles);
}

void accelerate_flow(const param_t params, const accel_area_t accel_area,
    speed_t* cells, int* obstacles)
{
    int ii,jj;     /* generic counters */
    double w1,w2;  /* weighting factors */

    /* compute weighting factors */
    w1 = params.density * params.accel / 9.0;
    w2 = params.density * params.accel / 36.0;

    if (accel_area.col_or_row == ACCEL_COLUMN)
    {
        jj = accel_area.idx;

        for (ii = 0; ii < params.ny; ii++)
        {
            /* if the cell is not occupied and
            ** we don't send a density negative */
            if (!obstacles[ii*params.nx + jj] &&
            (cells[ii*params.nx + jj].speeds[4] - w1) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[7] - w2) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[8] - w2) > 0.0 )
            {
                /* increase 'north-side' densities */
                cells[ii*params.nx + jj].speeds[2] += w1;
                cells[ii*params.nx + jj].speeds[5] += w2;
                cells[ii*params.nx + jj].speeds[6] += w2;
                /* decrease 'south-side' densities */
                cells[ii*params.nx + jj].speeds[4] -= w1;
                cells[ii*params.nx + jj].speeds[7] -= w2;
                cells[ii*params.nx + jj].speeds[8] -= w2;
            }
        }
    }
    else
    {
        ii = accel_area.idx;

        for (jj = 0; jj < params.nx; jj++)
        {
            /* if the cell is not occupied and
            ** we don't send a density negative */
            if (!obstacles[ii*params.nx + jj] &&
            (cells[ii*params.nx + jj].speeds[3] - w1) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[6] - w2) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[7] - w2) > 0.0 )
            {
                /* increase 'east-side' densities */
                cells[ii*params.nx + jj].speeds[1] += w1;
                cells[ii*params.nx + jj].speeds[5] += w2;
                cells[ii*params.nx + jj].speeds[8] += w2;
                /* decrease 'west-side' densities */
                cells[ii*params.nx + jj].speeds[3] -= w1;
                cells[ii*params.nx + jj].speeds[6] -= w2;
                cells[ii*params.nx + jj].speeds[7] -= w2;
            }
        }
    }
}

void propagate(const param_t params, speed_t* cells, speed_t* tmp_cells)
{
    int ii,jj;            /* generic counters */

    /* loop over _all_ cells */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            int x_e,x_w,y_n,y_s;  /* indices of neighbouring cells */
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
            y_n = (ii + 1) % params.ny;
            x_e = (jj + 1) % params.nx;
            y_s = (ii == 0) ? (ii + params.ny - 1) : (ii - 1);
            x_w = (jj == 0) ? (jj + params.nx - 1) : (jj - 1);
            /* propagate densities to neighbouring cells, following
            ** appropriate directions of travel and writing into
            ** scratch space grid */
            tmp_cells[ii *params.nx + jj].speeds[0]  = cells[ii*params.nx + jj].speeds[0]; /* central cell, */
                                                     /* no movement   */
            tmp_cells[ii *params.nx + x_e].speeds[1] = cells[ii*params.nx + jj].speeds[1]; /* east */
            tmp_cells[y_n*params.nx + jj].speeds[2]  = cells[ii*params.nx + jj].speeds[2]; /* north */
            tmp_cells[ii *params.nx + x_w].speeds[3] = cells[ii*params.nx + jj].speeds[3]; /* west */
            tmp_cells[y_s*params.nx + jj].speeds[4]  = cells[ii*params.nx + jj].speeds[4]; /* south */
            tmp_cells[y_n*params.nx + x_e].speeds[5] = cells[ii*params.nx + jj].speeds[5]; /* north-east */
            tmp_cells[y_n*params.nx + x_w].speeds[6] = cells[ii*params.nx + jj].speeds[6]; /* north-west */
            tmp_cells[y_s*params.nx + x_w].speeds[7] = cells[ii*params.nx + jj].speeds[7]; /* south-west */
            tmp_cells[y_s*params.nx + x_e].speeds[8] = cells[ii*params.nx + jj].speeds[8]; /* south-east */
        }
    }
}

void rebound(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,jj;  /* generic counters */

    /* loop over the cells in the grid */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* if the cell contains an obstacle */
            if (obstacles[ii*params.nx + jj])
            {
                /* called after propagate, so taking values from scratch space
                ** mirroring, and writing into main grid */
                cells[ii*params.nx + jj].speeds[1] = tmp_cells[ii*params.nx + jj].speeds[3];
                cells[ii*params.nx + jj].speeds[2] = tmp_cells[ii*params.nx + jj].speeds[4];
                cells[ii*params.nx + jj].speeds[3] = tmp_cells[ii*params.nx + jj].speeds[1];
                cells[ii*params.nx + jj].speeds[4] = tmp_cells[ii*params.nx + jj].speeds[2];
                cells[ii*params.nx + jj].speeds[5] = tmp_cells[ii*params.nx + jj].speeds[7];
                cells[ii*params.nx + jj].speeds[6] = tmp_cells[ii*params.nx + jj].speeds[8];
                cells[ii*params.nx + jj].speeds[7] = tmp_cells[ii*params.nx + jj].speeds[5];
                cells[ii*params.nx + jj].speeds[8] = tmp_cells[ii*params.nx + jj].speeds[6];
            }
        }
    }
}

void collision(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,jj,kk;                 /* generic counters */
    const double c_sq = 1.0/3.0;  /* square of speed of sound */
    const double w0 = 4.0/9.0;    /* weighting factor */
    const double w1 = 1.0/9.0;    /* weighting factor */
    const double w2 = 1.0/36.0;   /* weighting factor */

    double u_x,u_y;               /* av. velocities in x and y directions */
    double u_sq;                  /* squared velocity */
    double local_density;         /* sum of densities in a particular cell */
    double u[NSPEEDS];            /* directional velocities */
    double d_equ[NSPEEDS];        /* equilibrium densities */

    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* don't consider occupied cells */
            if (!obstacles[ii*params.nx + jj])
            {
                /* compute local density total */
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += tmp_cells[ii*params.nx + jj].speeds[kk];
                }

                /* compute x velocity component */
                u_x = (tmp_cells[ii*params.nx + jj].speeds[1] +
                        tmp_cells[ii*params.nx + jj].speeds[5] +
                        tmp_cells[ii*params.nx + jj].speeds[8]
                    - (tmp_cells[ii*params.nx + jj].speeds[3] +
                        tmp_cells[ii*params.nx + jj].speeds[6] +
                        tmp_cells[ii*params.nx + jj].speeds[7]))
                    / local_density;

                /* compute y velocity component */
                u_y = (tmp_cells[ii*params.nx + jj].speeds[2] +
                        tmp_cells[ii*params.nx + jj].speeds[5] +
                        tmp_cells[ii*params.nx + jj].speeds[6]
                    - (tmp_cells[ii*params.nx + jj].speeds[4] +
                        tmp_cells[ii*params.nx + jj].speeds[7] +
                        tmp_cells[ii*params.nx + jj].speeds[8]))
                    / local_density;

                /* velocity squared */
                u_sq = u_x * u_x + u_y * u_y;

                /* directional velocity components */
                u[1] =   u_x;        /* east */
                u[2] =         u_y;  /* north */
                u[3] = - u_x;        /* west */
                u[4] =       - u_y;  /* south */
                u[5] =   u_x + u_y;  /* north-east */
                u[6] = - u_x + u_y;  /* north-west */
                u[7] = - u_x - u_y;  /* south-west */
                u[8] =   u_x - u_y;  /* south-east */

                /* equilibrium densities */
                /* zero velocity density: weight w0 */
                d_equ[0] = w0 * local_density * (1.0 - u_sq / (2.0 * c_sq));
                /* axis speeds: weight w1 */
                d_equ[1] = w1 * local_density * (1.0 + u[1] / c_sq
                    + (u[1] * u[1]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[2] = w1 * local_density * (1.0 + u[2] / c_sq
                    + (u[2] * u[2]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[3] = w1 * local_density * (1.0 + u[3] / c_sq
                    + (u[3] * u[3]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[4] = w1 * local_density * (1.0 + u[4] / c_sq
                    + (u[4] * u[4]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                /* diagonal speeds: weight w2 */
                d_equ[5] = w2 * local_density * (1.0 + u[5] / c_sq
                    + (u[5] * u[5]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[6] = w2 * local_density * (1.0 + u[6] / c_sq
                    + (u[6] * u[6]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[7] = w2 * local_density * (1.0 + u[7] / c_sq
                    + (u[7] * u[7]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[8] = w2 * local_density * (1.0 + u[8] / c_sq
                    + (u[8] * u[8]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));

                /* relaxation step */
                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    cells[ii*params.nx + jj].speeds[kk] = 
                        (tmp_cells[ii*params.nx + jj].speeds[kk] + params.omega * 
                        (d_equ[kk] - tmp_cells[ii*params.nx + jj].speeds[kk]));
                }
            }
        }
    }
}

double av_velocity(const param_t params, speed_t* cells, int* obstacles)
{
    int    ii,jj,kk;       /* generic counters */
    int    tot_cells = 0;  /* no. of cells used in calculation */
    double tot_u;          /* accumulated magnitudes of velocity for each cell */

    double local_density;  /* total density in cell */
    double u_x;            /* x-component of velocity for current cell */
    double u_y;            /* y-component of velocity for current cell */

    /* initialise */
    tot_u = 0.0;

    /* loop over all non-blocked cells */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* ignore occupied cells */
            if (!obstacles[ii*params.nx + jj])
            {
                /* local density total */
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += cells[ii*params.nx + jj].speeds[kk];
                }

                /* x-component of velocity */
                u_x = (cells[ii*params.nx + jj].speeds[1] +
                        cells[ii*params.nx + jj].speeds[5] +
                        cells[ii*params.nx + jj].speeds[8]
                    - (cells[ii*params.nx + jj].speeds[3] +
                        cells[ii*params.nx + jj].speeds[6] +
                        cells[ii*params.nx + jj].speeds[7])) /
                    local_density;

                /* compute y velocity component */
                u_y = (cells[ii*params.nx + jj].speeds[2] +
                        cells[ii*params.nx + jj].speeds[5] +
                        cells[ii*params.nx + jj].speeds[6]
                    - (cells[ii*params.nx + jj].speeds[4] +
                        cells[ii*params.nx + jj].speeds[7] +
                        cells[ii*params.nx + jj].speeds[8])) /
                    local_density;

                /* accumulate the norm of x- and y- velocity components */
                tot_u += sqrt(u_x*u_x + u_y*u_y);
                /* increase counter of inspected cells */
                ++tot_cells;
            }
        }
    }

    return tot_u / (double)tot_cells;
}

