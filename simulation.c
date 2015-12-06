/* Functions pertinent to the outer simulation steps */

#include <math.h>
#include <stdlib.h> //for malloc
#include <stdio.h>  //for printf
#include <string.h> //for memcpy

#include "lbm.h"
#include "mpi.h"

void swap(speed_t** one, speed_t** two) {
  speed_t* temp = *one;
  *one = *two;
  *two = temp;
}

void timestep(const param_t params, const accel_area_t accel_area,
    speed_t* cells, speed_t* tmp_cells, int* obstacles)
{

    accelerate_flow(params,accel_area,cells,obstacles);
    propagate(params,cells,tmp_cells);
    int ii,jj,kk;
    for (ii = 0; ii < params.loc_ny; ii++)
    {
      for (jj = 0; jj < params.loc_nx; jj++)
      {
        for(kk = 0; kk < NSPEEDS; kk++) {
          if(tmp_cells[ii*params.loc_nx + jj].speeds[kk] == -1.0) {
            printf("MINUS ONE IN TMP AT %d %d %d\n", ii, jj, kk);
          }
        }
      }
    }
    rebound(params,cells,tmp_cells,obstacles);
    collision(params,cells,tmp_cells,obstacles);
    //halo_exchange(params,cells,tmp_cells);
    //swap(&cells,&tmp_cells);
    //if(params.rank == MASTER) printf("\n");
}

//TODO: set indices may not be in local range??
void accelerate_flow(const param_t params, const accel_area_t accel_area,
    speed_t* cells, int* obstacles)
{
    int ii,jj;     // generic counters
    float w1,w2;  // weighting factors

    // compute weighting factors
    w1 = params.density * params.accel / 9.0;
    w2 = params.density * params.accel / 36.0;

    if (accel_area.col_or_row == ACCEL_COLUMN)
    {
        jj = accel_area.idx;

        for (ii = 0; ii < params.loc_ny; ii++)
        {
            // if the cell is not occupied and
            // we don't send a density negative
            if (!obstacles[ii*params.nx + (jj + (params.rank * params.loc_nx))] &&
            (cells[ii*params.loc_nx + jj].speeds[4] - w1) > 0.0 &&
            (cells[ii*params.loc_nx + jj].speeds[7] - w2) > 0.0 &&
            (cells[ii*params.loc_nx + jj].speeds[8] - w2) > 0.0 )
            {
                // increase 'north-side' densities
                cells[ii*params.loc_nx + jj].speeds[2] += w1;
                cells[ii*params.loc_nx + jj].speeds[5] += w2;
                cells[ii*params.loc_nx + jj].speeds[6] += w2;
                // decrease 'south-side' densities
                cells[ii*params.loc_nx + jj].speeds[4] -= w1;
                cells[ii*params.loc_nx + jj].speeds[7] -= w2;
                cells[ii*params.loc_nx + jj].speeds[8] -= w2;
            }
        }
    }
    else
    {
        ii = accel_area.idx;

        for (jj = 0; jj < params.loc_nx; jj++)
        {
            // if the cell is not occupied and
            // we don't send a density negative
            if (!obstacles[ii*params.nx + (jj + (params.rank * params.loc_nx))] &&
            (cells[ii*params.loc_nx + jj].speeds[3] - w1) > 0.0 &&
            (cells[ii*params.loc_nx + jj].speeds[6] - w2) > 0.0 &&
            (cells[ii*params.loc_nx + jj].speeds[7] - w2) > 0.0 )
            {
                // increase 'east-side' densities
                cells[ii*params.loc_nx + jj].speeds[1] += w1;
                cells[ii*params.loc_nx + jj].speeds[5] += w2;
                cells[ii*params.loc_nx + jj].speeds[8] += w2;
                // decrease 'west-side' densities
                cells[ii*params.loc_nx + jj].speeds[3] -= w1;
                cells[ii*params.loc_nx + jj].speeds[6] -= w2;
                cells[ii*params.loc_nx + jj].speeds[7] -= w2;
            }
        }
    }
}

void propagate(const param_t params, speed_t* cells, speed_t* tmp_cells)
{
    int ii,jj;  // generic counters
    float* sendbuf_l = (float*)malloc(sizeof(float) * params.loc_ny * 3);
    float* sendbuf_r = (float*)malloc(sizeof(float) * params.loc_ny * 3);
    float* recvbuf_l = (float*)malloc(sizeof(float) * params.loc_ny * 3);
    float* recvbuf_r = (float*)malloc(sizeof(float) * params.loc_ny * 3);

    // loop over local cells
    for (ii = 0; ii < params.loc_ny; ii++)
    {
        for (jj = 0; jj < params.loc_nx; jj++)
        {
            int x_e,x_w,y_n,y_s;  // indices of neighbouring cells
            // determine indices of axis-direction neighbours
            // which will be in the halos if 'out of bounds' in the x-direction
            // or wrap around if 'out of bounds' in the y-direction
            y_n = (ii + 1) % params.loc_ny;
            x_e = ((jj) + 1);
            y_s = (ii == 0) ? (ii + params.loc_ny - 1) : (ii - 1);
            x_w = ((jj) - 1);
            /* propagate densities to neighbouring cells, following
            ** appropriate directions of travel and writing into
            ** scratch space grid */
            if(jj == 0)
            {
              tmp_cells[ii *params.loc_nx + (jj)].speeds[0] = cells[ii*params.loc_nx + (jj)].speeds[0]; // central cell
              tmp_cells[ii *params.loc_nx + x_e   ].speeds[1] = cells[ii*params.loc_nx + (jj)].speeds[1]; // east
              tmp_cells[y_n*params.loc_nx + (jj)].speeds[2] = cells[ii*params.loc_nx + (jj)].speeds[2]; // north
              tmp_cells[y_s*params.loc_nx + (jj)].speeds[4] = cells[ii*params.loc_nx + (jj)].speeds[4]; // south
              tmp_cells[y_n*params.loc_nx + x_e   ].speeds[5] = cells[ii*params.loc_nx + (jj)].speeds[5]; // north-east
              tmp_cells[y_s*params.loc_nx + x_e   ].speeds[8] = cells[ii*params.loc_nx + (jj)].speeds[8]; // south-east

              sendbuf_l[y_n * 3 + 0] = cells[ii*params.loc_nx + (jj)].speeds[6]; // north-west
              sendbuf_l[ii  * 3 + 1] = cells[ii*params.loc_nx + (jj)].speeds[3]; // west
              sendbuf_l[y_s * 3 + 2] = cells[ii*params.loc_nx + (jj)].speeds[7]; // south-west
            } else if(jj == params.loc_nx - 1) {
              tmp_cells[ii *params.loc_nx + (jj)].speeds[0] = cells[ii*params.loc_nx + (jj)].speeds[0]; // central cell
              tmp_cells[y_n*params.loc_nx + (jj)].speeds[2] = cells[ii*params.loc_nx + (jj)].speeds[2]; // north
              tmp_cells[ii *params.loc_nx + x_w   ].speeds[3] = cells[ii*params.loc_nx + (jj)].speeds[3]; // west
              tmp_cells[y_s*params.loc_nx + (jj)].speeds[4] = cells[ii*params.loc_nx + (jj)].speeds[4]; // south
              tmp_cells[y_n*params.loc_nx + x_w   ].speeds[6] = cells[ii*params.loc_nx + (jj)].speeds[6]; // north-west
              tmp_cells[y_s*params.loc_nx + x_w   ].speeds[7] = cells[ii*params.loc_nx + (jj)].speeds[7]; // south-west

              sendbuf_r[y_n * 3 + 0] = cells[ii*params.loc_nx + (jj)].speeds[5]; // north-east
              sendbuf_r[ii  * 3 + 1] = cells[ii*params.loc_nx + (jj)].speeds[1]; // east
              sendbuf_r[y_s * 3 + 2] = cells[ii*params.loc_nx + (jj)].speeds[8]; // south-east
            } else {
              tmp_cells[ii *params.loc_nx + (jj)].speeds[0] = cells[ii*params.loc_nx + (jj)].speeds[0]; // central cell
              tmp_cells[ii *params.loc_nx + x_e   ].speeds[1] = cells[ii*params.loc_nx + (jj)].speeds[1]; // east
              tmp_cells[y_n*params.loc_nx + (jj)].speeds[2] = cells[ii*params.loc_nx + (jj)].speeds[2]; // north
              tmp_cells[ii *params.loc_nx + x_w   ].speeds[3] = cells[ii*params.loc_nx + (jj)].speeds[3]; // west
              tmp_cells[y_s*params.loc_nx + (jj)].speeds[4] = cells[ii*params.loc_nx + (jj)].speeds[4]; // south
              tmp_cells[y_n*params.loc_nx + x_e   ].speeds[5] = cells[ii*params.loc_nx + (jj)].speeds[5]; // north-east
              tmp_cells[y_n*params.loc_nx + x_w   ].speeds[6] = cells[ii*params.loc_nx + (jj)].speeds[6]; // north-west
              tmp_cells[y_s*params.loc_nx + x_w   ].speeds[7] = cells[ii*params.loc_nx + (jj)].speeds[7]; // south-west
              tmp_cells[y_s*params.loc_nx + x_e   ].speeds[8] = cells[ii*params.loc_nx + (jj)].speeds[8]; // south-east
            }
            //} else {
            //  tmp_cells[ii *params.loc_nx + (jj)].speeds[0] = cells[ii*params.loc_nx + (jj)].speeds[0]; // central cell
            //  tmp_cells[y_n*params.loc_nx + (jj)].speeds[2] = cells[ii*params.loc_nx + (jj)].speeds[2]; // north
            //  tmp_cells[y_s*params.loc_nx + (jj)].speeds[4] = cells[ii*params.loc_nx + (jj)].speeds[4]; // south
            //}
        }
    }

    int tag = 0;           // scope for adding extra information to a message
    MPI_Status status;     // struct used by MPI_Recv

    // send to the left, receive from right
    MPI_Sendrecv(sendbuf_l, params.loc_ny * 3, MPI_FLOAT, params.left, tag,
                 recvbuf_r, params.loc_ny * 3, MPI_FLOAT, params.right,tag, MPI_COMM_WORLD, &status);

    for(ii = 0; ii < params.loc_ny; ii++)
    {
      tmp_cells[ii*params.loc_nx + params.loc_nx - 1].speeds[6] = recvbuf_r[ii * 3 + 0];
      tmp_cells[ii*params.loc_nx + params.loc_nx - 1].speeds[3] = recvbuf_r[ii * 3 + 1];
      tmp_cells[ii*params.loc_nx + params.loc_nx - 1].speeds[7] = recvbuf_r[ii * 3 + 2];
    }

    // send to the right, receive from left
    MPI_Sendrecv(sendbuf_r, params.loc_ny * 3, MPI_FLOAT, params.right, tag,
                 recvbuf_l, params.loc_ny * 3, MPI_FLOAT, params.left,  tag, MPI_COMM_WORLD, &status);
    for(ii = 0; ii < params.loc_ny; ii++)
    {
      tmp_cells[ii*params.loc_nx + 0].speeds[5] = recvbuf_l[ii * 3 + 0];
      tmp_cells[ii*params.loc_nx + 0].speeds[1] = recvbuf_l[ii * 3 + 1];
      tmp_cells[ii*params.loc_nx + 0].speeds[8] = recvbuf_l[ii * 3 + 2];
    }


    // Communicate changes in tmp_cells to other ranks
/*
    int tag = 0;           // scope for adding extra information to a message
    MPI_Status status;     // struct used by MPI_Recv

    // send to the left, receive from right
    //if(params.rank == 2) printf("LEFT HALO :: LEFT VALUES (TO SEND):\n");
    for(ii = 0; ii < params.loc_ny; ii++)
    {
      //send to the left the left side of the left halo
      memcpy(params.sendbuf[ii].speeds, tmp_cells[ii*(params.loc_nx+2) + 0].speeds, sizeof(float)*NSPEEDS);
      //if(params.rank == 2) printf("%f %f %f\n", params.sendbuf[ii].speeds[6], params.sendbuf[ii].speeds[3], params.sendbuf[ii].speeds[7]);
    }
    for(ii = 0; ii < params.loc_ny; ii++) {
      MPI_Sendrecv(params.sendbuf[ii].speeds, NSPEEDS, MPI_FLOAT, params.left, tag,
                   params.recvbuf[ii].speeds, NSPEEDS, MPI_FLOAT, params.right,tag, MPI_COMM_WORLD, &status);
    }
    //if(params.rank == 2) printf("RIGHT HALO :: LEFT VALUES (RECEIVED):\n");
    for(ii = 0; ii < params.loc_ny; ii++)
    {
      //memcpy(cells[ii*(params.loc_nx+2) + params.loc_nx + 1].speeds, params.recvbuf[ii].speeds, sizeof(float)*NSPEEDS);
      //receive on the right the left side of the right halo
      //if(params.rank == 2) printf("%f %f %f\n", params.recvbuf[ii].speeds[6], params.recvbuf[ii].speeds[3], params.recvbuf[ii].speeds[7]);
      tmp_cells[ii*(params.loc_nx+2) + params.loc_nx + 1].speeds[6] = params.recvbuf[ii].speeds[6];
      tmp_cells[ii*(params.loc_nx+2) + params.loc_nx + 1].speeds[3] = params.recvbuf[ii].speeds[3];
      tmp_cells[ii*(params.loc_nx+2) + params.loc_nx + 1].speeds[7] = params.recvbuf[ii].speeds[7];
    }

    // send to the right, receive from left
    //if(params.rank == 2) printf("RIGHT HALO :: RIGHT VALUES (TO SEND):\n");
    for(ii = 0; ii < params.loc_ny; ii++)
    {
      // send to the right the right side of the right halo
      memcpy(params.sendbuf[ii].speeds, tmp_cells[ii*(params.loc_nx+2) + params.loc_nx + 1].speeds, sizeof(float)*NSPEEDS);
      //if(params.rank == 2) printf("%f %f %f\n", params.sendbuf[ii].speeds[5], params.sendbuf[ii].speeds[1], params.sendbuf[ii].speeds[8]);
    }
    for(ii = 0; ii < params.loc_ny; ii++) {
      MPI_Sendrecv(params.sendbuf[ii].speeds, NSPEEDS, MPI_FLOAT, params.right, tag,
                   params.recvbuf[ii].speeds, NSPEEDS, MPI_FLOAT, params.left,  tag, MPI_COMM_WORLD, &status);
    }
    //if(params.rank == 2) printf("LEFT HALO :: RIGHT VALUES (RECEIVED):\n");
    for(ii = 0; ii < params.loc_ny; ii++)
    {
      //memcpy(cells[ii*(params.loc_nx+2) + 0].speeds, params.recvbuf[ii].speeds, sizeof(float)*NSPEEDS);
      // receive on the left the right side of the left halo
      //if(params.rank == 2) printf("%f %f %f\n", params.recvbuf[ii].speeds[5], params.recvbuf[ii].speeds[1], params.recvbuf[ii].speeds[8]);
      tmp_cells[ii*(params.loc_nx+2) + 0].speeds[5] = params.recvbuf[ii].speeds[5];
      tmp_cells[ii*(params.loc_nx+2) + 0].speeds[1] = params.recvbuf[ii].speeds[1];
      tmp_cells[ii*(params.loc_nx+2) + 0].speeds[8] = params.recvbuf[ii].speeds[8];
    }
    /*
    int kk;
    for(ii = 0; ii < params.loc_ny; ii++)
    {
      for(kk = 0; kk < NSPEEDS; kk++) {
        if(params.rank == 2) printf("%f ", tmp_cells[ii*(params.loc_nx+2)+0].speeds[kk]);
      }
      if(params.rank == 2) printf("\t\t");
      for(kk = 0; kk < NSPEEDS; kk++) {
        if(params.rank == 2) printf("%f ", tmp_cells[ii*(params.loc_nx+2)+params.loc_nx+1].speeds[kk]);
      }
      if(params.rank == 2) printf("\n");
    }
    if(params.rank == 2) printf("\n");*/

}

void rebound(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,jj;  /* generic counters */

    /* loop over the cells in the grid */
    for (ii = 0; ii < params.loc_ny; ii++)
    {
        for (jj = 0; jj < params.loc_nx; jj++)
        {
            /* if the cell contains an obstacle */
            //int obs_jj = jj-1;
            //int obs_idx = (obs_jj < 0 && params.rank == MASTER) ? params.nx-1 : (ii*params.nx + (obs_jj + (params.rank * params.loc_nx))) % params.nx;
            if (obstacles[ii*params.nx + (jj + (params.rank * params.loc_nx))])//ii*params.nx + (obs_jj + (params.rank * params.loc_nx))])
            {
                /* called after propagate, so taking values from scratch space
                ** mirroring, and writing into main grid */
                cells[ii*params.loc_nx + (jj)].speeds[1] = tmp_cells[ii*params.loc_nx + (jj)].speeds[3];
                cells[ii*params.loc_nx + (jj)].speeds[2] = tmp_cells[ii*params.loc_nx + (jj)].speeds[4];
                cells[ii*params.loc_nx + (jj)].speeds[3] = tmp_cells[ii*params.loc_nx + (jj)].speeds[1];
                cells[ii*params.loc_nx + (jj)].speeds[4] = tmp_cells[ii*params.loc_nx + (jj)].speeds[2];
                cells[ii*params.loc_nx + (jj)].speeds[5] = tmp_cells[ii*params.loc_nx + (jj)].speeds[7];
                cells[ii*params.loc_nx + (jj)].speeds[6] = tmp_cells[ii*params.loc_nx + (jj)].speeds[8];
                cells[ii*params.loc_nx + (jj)].speeds[7] = tmp_cells[ii*params.loc_nx + (jj)].speeds[5];
                cells[ii*params.loc_nx + (jj)].speeds[8] = tmp_cells[ii*params.loc_nx + (jj)].speeds[6];
            }
        }
    }
}

void collision(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,jj,kk;                 /* generic counters */
    const float c_sq = 1.0/3.0;  /* square of speed of sound */
    const float w0 = 4.0/9.0;    /* weighting factor */
    const float w1 = 1.0/9.0;    /* weighting factor */
    const float w2 = 1.0/36.0;   /* weighting factor */

    float u_x,u_y;               /* av. velocities in x and y directions */
    float u_sq;                  /* squared velocity */
    float local_density;         /* sum of densities in a particular cell */
    float u[NSPEEDS];            /* directional velocities */
    float d_equ[NSPEEDS];        /* equilibrium densities */

    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    for (ii = 0; ii < params.loc_ny; ii++)
    {
        for (jj = 0; jj < params.loc_nx; jj++)
        {
            /* don't consider occupied cells */
            if (!obstacles[ii*params.nx + (jj + (params.rank * params.loc_nx))])
            {
                /* compute local density total */
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += tmp_cells[ii*params.loc_nx + jj].speeds[kk];
                }

                /* compute x velocity component */
                u_x = (tmp_cells[ii*params.loc_nx + jj].speeds[1] +
                        tmp_cells[ii*params.loc_nx + jj].speeds[5] +
                        tmp_cells[ii*params.loc_nx + jj].speeds[8]
                    - (tmp_cells[ii*params.loc_nx + jj].speeds[3] +
                        tmp_cells[ii*params.loc_nx + jj].speeds[6] +
                        tmp_cells[ii*params.loc_nx + jj].speeds[7]))
                    / local_density;

                /* compute y velocity component */
                u_y = (tmp_cells[ii*params.loc_nx + jj].speeds[2] +
                        tmp_cells[ii*params.loc_nx + jj].speeds[5] +
                        tmp_cells[ii*params.loc_nx + jj].speeds[6]
                    - (tmp_cells[ii*params.loc_nx + jj].speeds[4] +
                        tmp_cells[ii*params.loc_nx + jj].speeds[7] +
                        tmp_cells[ii*params.loc_nx + jj].speeds[8]))
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
                    cells[ii*params.loc_nx + jj].speeds[kk] =
                        (tmp_cells[ii*params.loc_nx + jj].speeds[kk] + params.omega *
                        (d_equ[kk] - tmp_cells[ii*params.loc_nx + jj].speeds[kk]));
                }
            }
        }
    }
}

float av_velocity(const param_t params, speed_t* cells, int* obstacles)
{
    int    ii,jj,kk;       // generic counters
    int    tot_cells = 0;  // no. of cells used in calculation
    float tot_u;          // accumulated magnitudes of velocity for each cell

    float local_density;  // total density in cell
    float u_x;            // x-component of velocity for current cell
    float u_y;            // y-component of velocity for current cell

    // initialise
    tot_u = 0.0;

    // loop over all non-blocked cells
    for (ii = 0; ii < params.loc_ny; ii++)
    {
        for (jj = 0; jj < params.loc_nx; jj++)
        {
            // ignore occupied cells
            if (!obstacles[ii*params.nx + (jj + (params.rank * params.loc_nx))])
            {
                // local density total
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += cells[ii*params.loc_nx + jj].speeds[kk];
                }

                // x-component of velocity
                u_x = (cells[ii*params.loc_nx + jj].speeds[1] +
                        cells[ii*params.loc_nx + jj].speeds[5] +
                        cells[ii*params.loc_nx + jj].speeds[8]
                    - (cells[ii*params.loc_nx + jj].speeds[3] +
                        cells[ii*params.loc_nx + jj].speeds[6] +
                        cells[ii*params.loc_nx + jj].speeds[7])) /
                    local_density;

                // compute y velocity component
                u_y = (cells[ii*params.loc_nx + jj].speeds[2] +
                        cells[ii*params.loc_nx + jj].speeds[5] +
                        cells[ii*params.loc_nx + jj].speeds[6]
                    - (cells[ii*params.loc_nx + jj].speeds[4] +
                        cells[ii*params.loc_nx + jj].speeds[7] +
                        cells[ii*params.loc_nx + jj].speeds[8])) /
                    local_density;

                // accumulate the norm of x- and y- velocity components
                tot_u += sqrt(u_x*u_x + u_y*u_y);
                // increase counter of inspected cells
                ++tot_cells;
            }
        }
    }
    return tot_u / (float)tot_cells;
}
