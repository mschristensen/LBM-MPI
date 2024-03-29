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

#include <omp.h>
#include "lbm.h"
#include "mpi.h"

void swap(speed_t** one, speed_t** two) {
  speed_t* temp = *one;
  *one = *two;
  *two = temp;
}

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
    speed_t* tmp_tmp_cells = NULL;/* scratch space */
    char*     obstacles = NULL;    /* grid indicating which cells are blocked */
    float*  av_vels   = NULL;     /* a record of the av. velocity computed for each timestep */

    int    ii;                    /*  generic counter */
    struct timeval timstr;        /* structure to hold elapsed time */
    struct rusage ru;             /* structure to hold CPU time--system and user */
    double tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
    float usrtim;                /* floating point number to record elapsed user CPU time */
    float systim;                /* floating point number to record elapsed system CPU time */

    /********************************************/
    /******** Initialize MPI environment ********/
    /********************************************/
    //omp_set_num_threads(2);
    MPI_Init(&argc, &argv);
    /*int provided;
    //omp_set_num_threads(8);
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(provided != MPI_THREAD_FUNNELED) {
      printf("Cannot support funneled threading! Provided: %d :: %d %d %d %d\n", provided, MPI_THREAD_SINGLE, MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }*/

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

    /* determine the SIZE of the group of processes associated with
    ** the 'communicator'.  MPI_COMM_WORLD is the default communicator
    ** consisting of all the processes in the launched MPI 'job' */
    MPI_Comm_size( MPI_COMM_WORLD, &(params.size) );

    /* determine the RANK of the current process [0:SIZE-1] */
    MPI_Comm_rank( MPI_COMM_WORLD, &(params.rank) );

    MPI_Datatype mpi_speed_type;
    const int nblocks = 1;
    int blocklengths[1] = { NSPEEDS };
    MPI_Datatype types[1] = { MPI_FLOAT };
    MPI_Aint offsets[1];
    offsets[0] = offsetof(speed_t, speeds);
    MPI_Type_create_struct(nblocks, blocklengths, offsets, types, &mpi_speed_type);
    MPI_Type_commit(&mpi_speed_type);

    // LBM initialization
    parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);
    initialise(param_file, &accel_area, &params, &obstacles, &av_vels);

    /* determine process ranks to the left and right of rank
    ** respecting periodic boundary conditions */
    params.down = (params.rank == MASTER) ? (params.rank + params.size - 1) : (params.rank - 1);
    params.up = (params.rank + 1) % params.size;

    /* determine local grid size
    ** each rank gets all the rows, but a subset of the number of columns */

    for(ii = 0; ii < params.size; ii++) {
      params.loc_nys[ii] = calc_nrows_from_rank(params, ii);
      params.loc_nxs[ii] = params.nx;
    }
    params.loc_nx = params.loc_nxs[params.rank];
    params.loc_ny = params.loc_nys[params.rank];

    if (params.loc_ny < 1) {
      fprintf(stderr,"Error: too many processes:- local_ncols < 1\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Allocate the local arrays
    allocateLocal(&params, &cells, &tmp_cells, &tmp_tmp_cells);

    printf("Host %s: process %d of %d :: local_cells of size %dx%d\n", hostname, params.rank, params.size, params.loc_ny, params.loc_nx);

    /*
    // CODE TO INSPECT THE MPI + OPENMP SETUP
    #pragma omp parallel
    {
      printf("Host %s:\tProcess: %d/%d\tThread: %d/%d\n", hostname, params.rank, params.size, omp_get_thread_num(), omp_get_max_threads());
    }
    // Finalize MPI environment.
    MPI_Finalize();

    MPI_Finalized(&flag);
    if(flag != TRUE) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    return -1;
    */
    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

    for (ii = 0; ii < params.max_iters; ii++)
    {
      av_vels[ii] = timestep(params, accel_area, cells, tmp_cells, tmp_tmp_cells, obstacles);
      swap(&cells, &tmp_tmp_cells);

      #ifdef DEBUG
      printf("==timestep: %d==\n", ii);
      printf("av velocity: %.12E\n", av_vels[ii]);
      printf("tot density: %.12E\n", total_density(params, cells));
      #endif
    }

    gettimeofday(&timstr,NULL);
    toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    getrusage(RUSAGE_SELF, &ru);
    timstr=ru.ru_utime;
    usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    timstr=ru.ru_stime;
    systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

    // Final read back of all cell data
    speed_t* final_cells;
    if(params.rank == MASTER) final_cells = (speed_t*)malloc(sizeof(speed_t) * params.nx * params.ny);
    int* rcounts = (int*)malloc(sizeof(int) * params.size);
    int* rdispls = (int*)malloc(sizeof(int) * params.size);
    for(ii = 0; ii < params.size; ii++)
    {
      rcounts[ii] = params.loc_nys[ii] * params.loc_nxs[ii];
      if(ii == 0)
      {
        rdispls[ii] = 0;
      } else {
        rdispls[ii] = rdispls[ii-1] + rcounts[ii-1];
      }
    }
    MPI_Gatherv(cells, params.loc_nx * params.loc_ny, mpi_speed_type, final_cells, rcounts, rdispls, mpi_speed_type, MASTER, MPI_COMM_WORLD);


    if(params.rank == MASTER)
    {
      printf("Writing results...\n");
      write_values(final_state_file, av_vels_file, params, final_cells, obstacles, av_vels);

      const float last_av_vel = av_vels[params.max_iters - 1];
      printf("==done==\n");
      printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params,last_av_vel));
      printf("Elapsed time:\t\t\t%.6lf (s)\n", toc-tic);
      printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
      printf("Elapsed system CPU time:\t%.6f (s)\n", systim);
    } else {
      printf("Host %s: process %d of %d :: Elapsed time:\t\t\t%.6f (s)\n", hostname, params.rank, params.size, toc-tic);
    }

    // Finalize MPI environment.
    printf("finalising %d\n", params.rank);
    MPI_Finalize();
    printf("finalised %d\n", params.rank);

    MPI_Finalized(&flag);
    if(flag != TRUE) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    finalise(&cells, &tmp_cells, &tmp_tmp_cells, &obstacles, &av_vels);
    if(params.rank == MASTER) free(final_cells);
    return EXIT_SUCCESS;
}

int calc_nrows_from_rank(const param_t params, const int rank)
{
  int nrows;
  int i;
  nrows = params.ny / params.size;       /* integer division */
  int rem = (params.ny % params.size);
  if (rem != 0) {  /* if there is a remainder */
    for(i = 0; i < rem; i++) {
      if(rank == i) {
        nrows += 1;
      }
      rem -= 1;
    }
  }
  return nrows;
}

void write_values(const char * final_state_file, const char * av_vels_file,
    const param_t params, speed_t* cells, char* obstacles, float* av_vels)
{
    FILE* fp;                     /* file pointer */
    int xx,yy,ii,jj,kk;                 /* generic counters */
    const float c_sq = 1.0/3.0;  /* sq. of speed of sound */
    float local_density;         /* per grid cell sum of densities */
    float pressure;              /* fluid pressure in grid cell */
    float u_x;                   /* x-component of velocity in grid cell */
    float u_y;                   /* y-component of velocity in grid cell */
    float u;                     /* norm--root of summed squares--of u_x and u_y */
    char obs;
    fp = fopen(final_state_file, "w");

    if (fp == NULL)
    {
        DIE("could not open file output file");
    }

    // Loop over original grid size
    for (yy = 0; yy < params.original_grid.y2 + 1; yy++)
    {
        for (xx = 0; xx < params.original_grid.x2 + 1; xx++)
        {
          // if this is a region inside the bounding box then write the cell data
          if( yy >= params.obs_bbox.y1 && yy < params.obs_bbox.y2 && xx >= params.obs_bbox.x1 && xx < params.obs_bbox.x2)
          {
            // translate from bbox to global coords
            ii = yy - params.obs_bbox.y1;
            jj = xx - params.obs_bbox.x1;

            /* an occupied cell */
            if (obstacles[ii*params.nx + jj])
            {
                u_x = u_y = u = 0.0;
                pressure = params.density * c_sq;
                obs = 1;
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
                obs = 0;
            }
          } else {
            // region outside bbox: consider as an obstacle
            u_x = u_y = u = 0.0;
            pressure = params.density * c_sq;
            obs = 1;
          }
          /* write to file */
          fprintf(fp,"%d %d %.12E %.12E %.12E %.12E %d\n",xx,yy,u_x,u_y,u,pressure,obs);
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

float calc_reynolds(const param_t params, const float last_av_vel)
{
    const float viscosity = 1.0 / 6.0 * (2.0 / params.omega - 1.0);

    return last_av_vel * params.reynolds_dim / viscosity;
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
