/* Functions pertinent to the outer simulation steps */

#include <math.h>
#include <stdlib.h> //for malloc
#include <stdio.h>  //for printf
#include <string.h> //for memcpy

#include <omp.h>
#include "lbm.h"
#include "mpi.h"

float timestep(const param_t params, const accel_area_t accel_area,
    speed_t* cells, speed_t* tmp_cells, speed_t* tmp_tmp_cells, char* obstacles)
{
    accelerate_flow(params,accel_area,cells,obstacles);
    return d2q9bgk(params,cells,tmp_cells,tmp_tmp_cells,obstacles);
}

void accelerate_flow(const param_t params, const accel_area_t accel_area,
    speed_t* cells, char* obstacles)
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
            if (!obstacles[get_global_y_coord(params, params.rank, ii) * params.nx + jj] &&
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
        int global_ii = get_global_y_coord(params, params.rank, 0);
        /*int global_ii = 0;
        int kk;
        for(kk = 0; kk < params.rank; kk++)
        {
          global_ii += params.loc_nys[kk];
        }*/
        if(ii > global_ii && ii < global_ii + params.loc_nys[params.rank])
        {
          for (jj = 0; jj < params.loc_nx; jj++)
          {
              // if the cell is not occupied and
              // we don't send a density negative
              if (!obstacles[ii * params.nx + jj] &&
              (cells[(ii - global_ii)*params.loc_nx + jj].speeds[3] - w1) > 0.0 &&
              (cells[(ii - global_ii)*params.loc_nx + jj].speeds[6] - w2) > 0.0 &&
              (cells[(ii - global_ii)*params.loc_nx + jj].speeds[7] - w2) > 0.0 )
              {
                  // increase 'east-side' densities
                  cells[(ii - global_ii)*params.loc_nx + jj].speeds[1] += w1;
                  cells[(ii - global_ii)*params.loc_nx + jj].speeds[5] += w2;
                  cells[(ii - global_ii)*params.loc_nx + jj].speeds[8] += w2;
                  // decrease 'west-side' densities
                  cells[(ii - global_ii)*params.loc_nx + jj].speeds[3] -= w1;
                  cells[(ii - global_ii)*params.loc_nx + jj].speeds[6] -= w2;
                  cells[(ii - global_ii)*params.loc_nx + jj].speeds[7] -= w2;
              }
          }
        }
    }
}

float d2q9bgk(const param_t params, speed_t* cells, speed_t* tmp_cells, speed_t* tmp_tmp_cells, char* obstacles)
{
  int ii,jj,kk,yy,ll;  // generic counters
  int x_e,x_w,y_n,y_s;  // indices of neighbouring cells

  float u_x,u_y;               /* av. velocities in x and y directions */
  float u_sq;                  /* squared velocity */
  float local_density;         /* sum of densities in a particular cell */
  float u[NSPEEDS];            /* directional velocities */
  float d_equ[NSPEEDS];        /* equilibrium densities */

  int tot_cells = 0;  // no. of cells used in calculation
  float tot_u = 0.0;          // accumulated magnitudes of velocity for each cell

  // Propagate the halos, and fill the buffers with the data to send
  for (jj = 0; jj < params.loc_nx; jj++)
  {
    // Bottom row
    ii = 0;
    y_n = (0 + 1);
    x_e = (jj + 1) % params.loc_nx;
    y_s = (0 - 1);
    x_w = (jj == 0) ? (jj + params.loc_nx - 1) : (jj - 1);

    tmp_cells[ii *params.loc_nx + jj ].speeds[0] = cells[ii *params.loc_nx + jj ].speeds[0]; // central cell
    tmp_cells[ii *params.loc_nx + jj ].speeds[1] = cells[ii *params.loc_nx + x_w].speeds[1]; // east
    tmp_cells[ii *params.loc_nx + jj ].speeds[4] = cells[y_n*params.loc_nx + jj ].speeds[4]; // south
    tmp_cells[ii *params.loc_nx + jj ].speeds[3] = cells[ii *params.loc_nx + x_e].speeds[3]; // west
    tmp_cells[ii *params.loc_nx + jj ].speeds[8] = cells[y_n*params.loc_nx + x_w].speeds[8]; // south east
    tmp_cells[ii *params.loc_nx + jj ].speeds[7] = cells[y_n*params.loc_nx + x_e].speeds[7]; // south west

    params.sendbuf_d[x_w * 3 + 0] = cells[ii*params.loc_nx + jj].speeds[7]; // south-west
    params.sendbuf_d[jj  * 3 + 1] = cells[ii*params.loc_nx + jj].speeds[4]; // south
    params.sendbuf_d[x_e * 3 + 2] = cells[ii*params.loc_nx + jj].speeds[8]; // south-east

    // Top row
    ii = params.loc_ny - 1;
    y_n = ((params.loc_ny - 1) + 1);
    x_e = (jj + 1) % params.loc_nx;
    y_s = ((params.loc_ny - 1) - 1);
    x_w = (jj == 0) ? (jj + params.loc_nx - 1) : (jj - 1);

    tmp_cells[ii *params.loc_nx + jj ].speeds[0] = cells[ii *params.loc_nx + jj ].speeds[0]; // central cell
    tmp_cells[ii *params.loc_nx + jj ].speeds[1] = cells[ii *params.loc_nx + x_w].speeds[1]; // east
    tmp_cells[ii *params.loc_nx + jj ].speeds[2] = cells[y_s*params.loc_nx + jj ].speeds[2]; // north
    tmp_cells[ii *params.loc_nx + jj ].speeds[3] = cells[ii *params.loc_nx + x_e].speeds[3]; // west
    tmp_cells[ii *params.loc_nx + jj ].speeds[5] = cells[y_s*params.loc_nx + x_w].speeds[5]; // north-east
    tmp_cells[ii *params.loc_nx + jj ].speeds[6] = cells[y_s*params.loc_nx + x_e].speeds[6]; // north-west

    params.sendbuf_u[x_w * 3 + 0] = cells[ii *params.loc_nx + jj].speeds[6]; // north-west
    params.sendbuf_u[jj  * 3 + 1] = cells[ii *params.loc_nx + jj].speeds[2]; // north
    params.sendbuf_u[x_e * 3 + 2] = cells[ii *params.loc_nx + jj].speeds[5]; // north-east
  }

  int tag = 0;           // scope for adding extra information to a message
  MPI_Status status;     // struct used by MPI_Recv
  MPI_Request request_d;   // request struct used in non-blocking comms calls
  MPI_Request request_u;   // request struct used in non-blocking comms calls

  // asynchronously send and receive the halo data above and below
  MPI_Isend(params.sendbuf_d, params.loc_nx * 3, MPI_FLOAT, params.down, tag, MPI_COMM_WORLD, &request_d);
  MPI_Isend(params.sendbuf_u, params.loc_nx * 3, MPI_FLOAT, params.up,   tag, MPI_COMM_WORLD, &request_u);
  MPI_Irecv(params.recvbuf_u, params.loc_nx * 3, MPI_FLOAT, params.up,   tag, MPI_COMM_WORLD, &request_d);
  MPI_Irecv(params.recvbuf_d, params.loc_nx * 3, MPI_FLOAT, params.down, tag, MPI_COMM_WORLD, &request_u);

  #pragma omp parallel default(none) shared(cells,tmp_cells,tmp_tmp_cells,obstacles,tot_u,tot_cells) private(ii,jj,kk,yy,ll,u_x,u_y,u_sq,local_density,u,d_equ,status,request_u,request_d) firstprivate(x_w,x_e,y_s,y_n)
  {
    #pragma omp for reduction(+:tot_u,tot_cells) schedule(static)
    // loop over all cells, but start *after* the first halo so as to do computation that
    // doesnt need halo data first. By the time this is done, hopefully wont have to wait
    // long in MPI_Wait for the async receive to have completed!
    for (yy = 1; yy < params.loc_ny - 1; yy++)
    {
      ii = yy;

      for (jj = 0; jj < params.loc_nx; jj++)
      {
        loop_body(params, cells, tmp_cells, tmp_tmp_cells, obstacles, ii, jj, &tot_cells, &tot_u);
      }
    }
  }

  for (yy = params.loc_ny - 1; yy < params.loc_ny + 1; yy++)
  {
    ii = yy;
    // If this is the last halo, wait to receive the data
    if(ii == params.loc_ny - 1) {
      MPI_Wait(&request_d, &status);
      for(ll = 0; ll < params.loc_nx; ll++)
      {
        tmp_cells[(params.loc_ny - 1) * params.loc_nx + ll].speeds[7] = params.recvbuf_u[ll * 3 + 0];
        tmp_cells[(params.loc_ny - 1) * params.loc_nx + ll].speeds[4] = params.recvbuf_u[ll * 3 + 1];
        tmp_cells[(params.loc_ny - 1) * params.loc_nx + ll].speeds[8] = params.recvbuf_u[ll * 3 + 2];
      }
    }
    // if this is the row "after" the last halo, we are on the *first* halo (which we originally missed)
    // So now we can wait until we have received the data
    else if(ii == params.loc_ny) {
      ii = 0;
      MPI_Wait(&request_u, &status);
      for(ll = 0; ll < params.loc_nx; ll++)
      {
        tmp_cells[ll].speeds[6] = params.recvbuf_d[ll * 3 + 0];
        tmp_cells[ll].speeds[2] = params.recvbuf_d[ll * 3 + 1];
        tmp_cells[ll].speeds[5] = params.recvbuf_d[ll * 3 + 2];
      }
    }

    for (jj = 0; jj < params.loc_nx; jj++)
    {
      loop_body(params, cells, tmp_cells, tmp_tmp_cells, obstacles, ii, jj, &tot_cells, &tot_u);
    }
  }

  // Reduction
  int total_cells = 0;
  MPI_Reduce(&tot_cells, &total_cells, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
  float total_u = 0.0;
  MPI_Reduce(&tot_u, &total_u, 1, MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);

  return total_u / (float)total_cells;
}


inline void loop_body(const param_t params, speed_t* cells, speed_t* tmp_cells, speed_t* tmp_tmp_cells, char* obstacles, int ii, int jj, int* tot_cells, float* tot_u)
{
  int kk;
  int y_n, x_e, y_s, x_w;

  float u_x,u_y;               /* av. velocities in x and y directions */
  float u_sq;                  /* squared velocity */
  float local_density;         /* sum of densities in a particular cell */
  float u[NSPEEDS];            /* directional velocities */
  float d_equ[NSPEEDS];        /* equilibrium densities */

  // Don't propagate the halos, they have already been done!
  if(ii != 0 && ii != params.loc_ny - 1) {
    // determine indices of axis-direction neighbours
    // which will be in the halos if 'out of bounds' in the x-direction
    // or wrap around if 'out of bounds' in the y-direction
    y_n = (ii + 1);
    x_e = (jj + 1) % params.loc_nx;
    y_s = (ii - 1);
    x_w = (jj == 0) ? (jj + params.loc_nx - 1) : (jj - 1);
    /* propagate densities to neighbouring cells, following
    ** appropriate directions of travel and writing into
    ** scratch space grid */
    tmp_cells[ii *params.loc_nx + jj ].speeds[0] = cells[ii *params.loc_nx + jj ].speeds[0]; // central cell
    tmp_cells[ii *params.loc_nx + jj ].speeds[1] = cells[ii *params.loc_nx + x_w].speeds[1]; // east
    tmp_cells[ii *params.loc_nx + jj ].speeds[2] = cells[y_s*params.loc_nx + jj ].speeds[2]; // north
    tmp_cells[ii *params.loc_nx + jj ].speeds[3] = cells[ii *params.loc_nx + x_e].speeds[3]; // west
    tmp_cells[ii *params.loc_nx + jj ].speeds[4] = cells[y_n*params.loc_nx + jj ].speeds[4]; // south
    tmp_cells[ii *params.loc_nx + jj ].speeds[5] = cells[y_s*params.loc_nx + x_w].speeds[5]; // north-east
    tmp_cells[ii *params.loc_nx + jj ].speeds[6] = cells[y_s*params.loc_nx + x_e].speeds[6]; // north-west
    tmp_cells[ii *params.loc_nx + jj ].speeds[7] = cells[y_n*params.loc_nx + x_e].speeds[7]; // south-west
    tmp_cells[ii *params.loc_nx + jj ].speeds[8] = cells[y_n*params.loc_nx + x_w].speeds[8]; // south-east
  }

  /* if the cell contains an obstacle */
  if (obstacles[get_global_y_coord(params, params.rank, ii) * params.nx + jj])
  {
      /* called after propagate, so taking values from scratch space
      ** mirroring, and writing into main grid */
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[1] = tmp_cells[ii*params.loc_nx + jj].speeds[3];
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[2] = tmp_cells[ii*params.loc_nx + jj].speeds[4];
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[3] = tmp_cells[ii*params.loc_nx + jj].speeds[1];
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[4] = tmp_cells[ii*params.loc_nx + jj].speeds[2];
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[5] = tmp_cells[ii*params.loc_nx + jj].speeds[7];
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[6] = tmp_cells[ii*params.loc_nx + jj].speeds[8];
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[7] = tmp_cells[ii*params.loc_nx + jj].speeds[5];
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[8] = tmp_cells[ii*params.loc_nx + jj].speeds[6];
  } else {
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

    /* zero velocity density: weight w0 */
    d_equ[0] = 0.44444444444 * local_density * (1.0 - u_sq * 1.5);
    /* axis speeds: weight w1 */
    d_equ[1] = 0.11111111111 * local_density * (1.0 + u[1] * 3.0
                     + (u[1] * u[1]) * 4.5
                     - u_sq * 1.5);
    d_equ[2] = 0.11111111111 * local_density * (1.0 + u[2] * 3.0
                     + (u[2] * u[2]) * 4.5
                     - u_sq * 1.5);
    d_equ[3] = 0.11111111111 * local_density * (1.0 + u[3] * 3.0
                     + (u[3] * u[3]) * 4.5
                     - u_sq * 1.5);
    d_equ[4] = 0.11111111111 * local_density * (1.0 + u[4] * 3.0
                     + (u[4] * u[4]) * 4.5
                     - u_sq * 1.5);
    /* diagonal speeds: weight w2 */
    d_equ[5] = 0.02777777777 * local_density * (1.0 + u[5] * 3.0
                     + (u[5] * u[5]) * 4.5
                     - u_sq * 1.5);
    d_equ[6] = 0.02777777777 * local_density * (1.0 + u[6] * 3.0
                     + (u[6] * u[6]) * 4.5
                     - u_sq * 1.5);
    d_equ[7] = 0.02777777777 * local_density * (1.0 + u[7] * 3.0
                     + (u[7] * u[7]) * 4.5
                     - u_sq * 1.5);
    d_equ[8] = 0.02777777777 * local_density * (1.0 + u[8] * 3.0
                     + (u[8] * u[8]) * 4.5
                     - u_sq * 1.5);

    // local density total
    local_density = 0.0;
    for (kk = 0; kk < NSPEEDS; kk++)
    {
      // relaxation step
      tmp_tmp_cells[ii*params.loc_nx + jj].speeds[kk] =
          (tmp_cells[ii*params.loc_nx + jj].speeds[kk] + params.omega *
          (d_equ[kk] - tmp_cells[ii*params.loc_nx + jj].speeds[kk]));

      // av_vels step
      local_density += tmp_tmp_cells[ii*params.loc_nx + jj].speeds[kk];
    }

    // x-component of velocity
    u_x = (tmp_tmp_cells[ii*params.loc_nx + jj].speeds[1] +
            tmp_tmp_cells[ii*params.loc_nx + jj].speeds[5] +
            tmp_tmp_cells[ii*params.loc_nx + jj].speeds[8]
        - (tmp_tmp_cells[ii*params.loc_nx + jj].speeds[3] +
            tmp_tmp_cells[ii*params.loc_nx + jj].speeds[6] +
            tmp_tmp_cells[ii*params.loc_nx + jj].speeds[7])) /
        local_density;

    // compute y velocity component
    u_y = (tmp_tmp_cells[ii*params.loc_nx + jj].speeds[2] +
            tmp_tmp_cells[ii*params.loc_nx + jj].speeds[5] +
            tmp_tmp_cells[ii*params.loc_nx + jj].speeds[6]
        - (tmp_tmp_cells[ii*params.loc_nx + jj].speeds[4] +
            tmp_tmp_cells[ii*params.loc_nx + jj].speeds[7] +
            tmp_tmp_cells[ii*params.loc_nx + jj].speeds[8])) /
        local_density;

    // accumulate the norm of x- and y- velocity components
    (*tot_u) += sqrt(u_x*u_x + u_y*u_y);
    // increase counter of inspected cells
    ++(*tot_cells);
  }
}
