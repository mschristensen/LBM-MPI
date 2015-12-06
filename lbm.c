/*
** code to implement a d2q9-bgk lattice boltzmann scheme.
** 'd2' inidates a 2-dimensional grid, and
** 'q9' indicates 9 velocities per grid cell.
** 'bgk' refers to the bhatnagar-gross-krook collision step.
**
** the 'speeds' in each cell are numbered as follows:
**
** 6 2 5
**  \|/
** 3-0-1
**  /|\
** 7 4 8
**
** a 2d grid:
**
**           cols
**       --- --- ---
**      | d | e | f |
** rows  --- --- ---
**      | a | b | c |
**       --- --- ---
**
** 'unwrapped' in row major order to give a 1d array:
**
**  --- --- --- --- --- ---
** | a | b | c | d | e | f |
**  --- --- --- --- --- ---
**
** grid indicies are:
**
**          ny
**          ^       cols(jj)
**          |  ----- ----- -----
**          | | ... | ... | etc |
**          |  ----- ----- -----
** rows(ii) | | 1,0 | 1,1 | 1,2 |
**          |  ----- ----- -----
**          | | 0,0 | 0,1 | 0,2 |
**          |  ----- ----- -----
**          ----------------------> nx
**
** note the names of the input parameter and obstacle files
** are passed on the command line, e.g.:
**
**   d2q9-bgk.exe input.params obstacles.dat
**
** be sure to adjust the grid dimensions in the parameter file
** if you choose a different obstacle file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h> //memcpy

#include "lbm.h"
#include "mpi.h"

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char* argv[])
{
    enum bool {FALSE,TRUE}; /* enumerated type: false = 0, true = 1 */
    char * final_state_file = NULL;
    char * av_vels_file = NULL;
    char * param_file = NULL;

    accel_area_t accel_area;

    param_t  params;              /* struct to hold parameter values */
    speed_t* cells     = NULL;    /* grid containing fluid densities */
    speed_t* tmp_cells = NULL;    /* scratch space */
    int*     obstacles = NULL;    /* grid indicating which cells are blocked */
    float*  av_vels   = NULL;    /* a record of the av. velocity computed for each timestep */

    int    ii;                    /*  generic counter */
    struct timeval timstr;        /* structure to hold elapsed time */
    struct rusage ru;             /* structure to hold CPU time--system and user */
    float tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
    float usrtim;                /* floating point number to record elapsed user CPU time */
    float systim;                /* floating point number to record elapsed system CPU time */

    parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);

    initialise(param_file, &accel_area, &params, &obstacles, &av_vels);

    // Initialize MPI environment.
    MPI_Init(&argc, &argv);
    int flag;
    // Check if initialization was successful.
    MPI_Initialized(&flag);
    if(flag != TRUE) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int strlen;             // length of a character array
    char hostname[MPI_MAX_PROCESSOR_NAME];  // character array to hold hostname running process

    // determine the hostname
    MPI_Get_processor_name(hostname,&strlen);

    /*
    ** determine the SIZE of the group of processes associated with
    ** the 'communicator'.  MPI_COMM_WORLD is the default communicator
    ** consisting of all the processes in the launched MPI 'job'
    */
    MPI_Comm_size( MPI_COMM_WORLD, &(params.size) );

    /* determine the RANK of the current process [0:SIZE-1] */
    MPI_Comm_rank( MPI_COMM_WORLD, &(params.rank) );

    /*
    ** determine process ranks to the left and right of rank
    ** respecting periodic boundary conditions
    */
    params.left = (params.rank == MASTER) ? (params.rank + params.size - 1) : (params.rank - 1);
    params.right = (params.rank + 1) % params.size;

    /*
    ** determine local grid size
    ** each rank gets all the rows, but a subset of the number of columns
    */
    params.loc_ny = params.ny;
    params.loc_nx = calc_ncols_from_rank(params);
    if (params.loc_nx < 1) {
      fprintf(stderr,"Error: too many processes:- local_ncols < 1\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Allocate the local arrays
    allocateLocal(&params, &cells, &tmp_cells);

    printf("Host %s: process %d of %d :: local_cells of size %dx%d plus halos\n", hostname, params.rank, params.size, params.loc_ny, params.loc_nx);

    /* iterate for max_iters timesteps */
    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);


    for (ii = 0; ii < params.max_iters; ii++)
    {
        //TODO: -av_vels reduction
        //      -last calculated av_vels in calc_reynolds
        //      -read back all cell data to big array for write_values
        timestep(params, accel_area, cells, tmp_cells, obstacles);
        float av_vel;
        av_vel = av_velocity(params, cells, obstacles);

        // Reduction
        MPI_Reduce(&av_vel, &(av_vels[ii]), 1, MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);
        av_vels[ii] /= 4.0;
        //if(ii == 5) break;

        #ifdef DEBUG
        printf("==timestep: %d==\n", ii);
        printf("av velocity: %.12E\n", av_vels[ii]);
        printf("tot density: %.12E\n", total_density(params, cells));
        #endif
    }

    // Final read back of all cell data
    //TODO this below malloc is made by all threads... but it needs to be freed
    speed_t* final_cells = (speed_t*)malloc(sizeof(speed_t) * params.nx * params.ny);
    int tag, jj, kk;
    MPI_Status status;     // struct used by MPI_Recv
    if(params.rank == MASTER) {
      // Master first copies his data to the final array
      for (ii = 0; ii < params.loc_ny; ii++)
      {
        for (jj = 0; jj < params.loc_nx; jj++)
        {
          memcpy(final_cells[ii*params.nx + (jj + (params.rank * params.loc_nx))].speeds, cells[ii*params.loc_nx + jj].speeds, sizeof(float)*NSPEEDS);
        }
      }

      // Then he reads data from each of the other ranks and writes that
      for(kk = 1; kk < params.size; kk++) {
        for (ii = 0; ii < params.loc_ny; ii++)
        {
          for (jj = 0; jj < params.loc_nx; jj++)
          {
            MPI_Recv(final_cells[ii*params.nx + (jj + (kk * params.loc_nx))].speeds, NSPEEDS, MPI_FLOAT, kk, tag, MPI_COMM_WORLD, &status);
          }
        }
      }
    } else {
      // Non-master ranks send all their data (not the halos) to master
      for (ii = 0; ii < params.loc_ny; ii++)
      {
        for (jj = 0; jj < params.loc_nx; jj++)
        {
          MPI_Send(cells[ii*params.loc_nx + jj].speeds, NSPEEDS, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
        }
      }
    }

    gettimeofday(&timstr,NULL);
    toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    getrusage(RUSAGE_SELF, &ru);
    timstr=ru.ru_utime;
    usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    timstr=ru.ru_stime;
    systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

    if(params.rank == MASTER)
    {
      printf("==done==\n");
      printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params,final_cells,obstacles));
      printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
      printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
      printf("Elapsed system CPU time:\t%.6f (s)\n", systim);

      printf("Writing results...\n");
      write_values(final_state_file, av_vels_file, params, final_cells, obstacles, av_vels);
    } else {
      printf("Host %s: process %d of %d :: Elapsed time:\t\t\t%.6f (s)\n", hostname, params.rank, params.size, toc-tic);
    }

    // Finalize MPI environment.
    MPI_Finalize();

    MPI_Finalized(&flag);
    if(flag != TRUE) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    finalise(&cells, &tmp_cells, &obstacles, &av_vels);
    free(final_cells);
    return EXIT_SUCCESS;
}

int calc_ncols_from_rank(const param_t params)
{
  int ncols;

  ncols = params.nx / params.size;       /* integer division */
  if ((params.nx % params.size) != 0) {  /* if there is a remainder */
    if (params.rank == params.size - 1)
      ncols += params.nx % params.size;  /* add remainder to last rank */
  }

  printf("Rank %d is has %d cols.\n", params.rank, ncols);

  return ncols;
}

void write_values(const char * final_state_file, const char * av_vels_file,
    const param_t params, speed_t* cells, int* obstacles, float* av_vels)
{
    FILE* fp;                     /* file pointer */
    int ii,jj,kk;                 /* generic counters */
    const float c_sq = 1.0/3.0;  /* sq. of speed of sound */
    float local_density;         /* per grid cell sum of densities */
    float pressure;              /* fluid pressure in grid cell */
    float u_x;                   /* x-component of velocity in grid cell */
    float u_y;                   /* y-component of velocity in grid cell */
    float u;                     /* norm--root of summed squares--of u_x and u_y */

    fp = fopen(final_state_file, "w");

    if (fp == NULL)
    {
        DIE("could not open file output file");
    }

    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* an occupied cell */
            if (obstacles[ii*params.nx + jj])
            {
                u_x = u_y = u = 0.0;
                pressure = params.density * c_sq;
            }
            /* no obstacle */
            else
            {
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += cells[ii*params.nx + jj].speeds[kk];
                }

                /* compute x velocity component */
                u_x = (cells[ii*params.nx + jj].speeds[1] +
                        cells[ii*params.nx + jj].speeds[5] +
                        cells[ii*params.nx + jj].speeds[8]
                    - (cells[ii*params.nx + jj].speeds[3] +
                        cells[ii*params.nx + jj].speeds[6] +
                        cells[ii*params.nx + jj].speeds[7]))
                    / local_density;

                /* compute y velocity component */
                u_y = (cells[ii*params.nx + jj].speeds[2] +
                        cells[ii*params.nx + jj].speeds[5] +
                        cells[ii*params.nx + jj].speeds[6]
                    - (cells[ii*params.nx + jj].speeds[4] +
                        cells[ii*params.nx + jj].speeds[7] +
                        cells[ii*params.nx + jj].speeds[8]))
                    / local_density;

                /* compute norm of velocity */
                u = sqrt((u_x * u_x) + (u_y * u_y));

                /* compute pressure */
                pressure = local_density * c_sq;
            }

            /* write to file */
            fprintf(fp,"%d %d %.12E %.12E %.12E %.12E %d\n",
                jj,ii,u_x,u_y,u,pressure,obstacles[ii*params.nx + jj]);
        }
    }

    fclose(fp);

    fp = fopen(av_vels_file, "w");
    if (fp == NULL)
    {
        DIE("could not open file output file");
    }

    for (ii = 0; ii < params.max_iters; ii++)
    {
        fprintf(fp,"%d:\t%.12E\n", ii, av_vels[ii]);
    }

    fclose(fp);
}

float calc_reynolds(const param_t params, speed_t* cells, int* obstacles)
{
    const float viscosity = 1.0 / 6.0 * (2.0 / params.omega - 1.0);

    return av_velocity(params,cells,obstacles) * params.reynolds_dim / viscosity;
}

float total_density(const param_t params, speed_t* cells)
{
    int ii,jj,kk;        /* generic counters */
    float total = 0.0;  /* accumulator */

    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.ny; jj++)
        {
            for (kk = 0; kk < NSPEEDS; kk++)
            {
                total += cells[ii*params.nx + jj].speeds[kk];
            }
        }
    }

    return total;
}
