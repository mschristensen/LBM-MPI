/* utilities for lbm to read files, etc */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>

#include "lbm.h"

void exit_with_error(int line, const char* filename, const char* format, ...)
{
    va_list arglist;

    fprintf(stderr, "Fatal error at line %d in %s: ", line, filename);

    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
    va_end(arglist);

    fprintf(stderr, "\n");

    exit(EXIT_FAILURE);
}

void parse_args (int argc, char* argv[],
    char** final_state_file, char** av_vels_file, char** param_file)
{
    int character;

    *av_vels_file = NULL;
    *final_state_file = NULL;
    *param_file = NULL;

    const char * help_msg =
    "usage: ./lbm [OPTIONS] \n"
    "   -a AV_VELS_FILE\n"
    "       Name of output average velocities file\n"
    "   -f FINAL_STATE_FILE\n"
    "       Name of output final state file\n"
    "   -p PARAM_FILE\n"
    "       Name of input parameter file\n"
    "   -h\n"
    "       Show this message and exit\n";

    /* Used getopt to parse command line arguments for filenames */
    while ((character = getopt(argc, argv, "a:f:p:h")) != -1)
    {
        switch (character)
        {
        case 'a':
            *av_vels_file = optarg;
            break;
        case 'f':
            *final_state_file = optarg;
            break;
        case 'p':
            *param_file = optarg;
            break;
        case 'h':
            fprintf(stderr, "%s", help_msg);
            exit(EXIT_SUCCESS);
            break;
        case '?':
            if (optopt == 'a' || optopt == 'f' || optopt == 'p')
            {
                /* Flag present, but no option specified */
                DIE("No argument specified for '%c'", optopt);
            }
            else if (isprint(optopt))
            {
                DIE("Unknown option %c", optopt);
            }
            break;
        default:
            DIE("Error in getopt");
        }
    }

    /* Make sure they were all present */
    if (NULL == *av_vels_file)
    {
        DIE("No argument specified for av_vels file");
    }
    if (NULL == *final_state_file)
    {
        DIE("No argument specified for final state file");
    }
    if (NULL == *param_file)
    {
        DIE("No argument specified for param file");
    }
}

void initialise(const char* param_file, accel_area_t * accel_area, param_t* params, char** obstacles_ptr, float** av_vels_ptr)
{
  int    ii,jj, kk;          /* generic counters */
  FILE   *fp;            /* file pointer */
  int    retval;         /* to hold return value for checking */

  /* Rectangular obstacles */
  int n_obstacles;
  obstacle_t * obstacles;

  fp = fopen(param_file, "r");

  if (NULL == fp)
  {
      DIE("Unable to open param file %s", param_file);
  }

  /* read in the parameter values */
  retval = fscanf(fp,"%d\n",&(params->nx));
  if (retval != 1) DIE("Could not read param file: nx");
  retval = fscanf(fp,"%d\n",&(params->ny));
  if (retval != 1) DIE("Could not read param file: ny");
  retval = fscanf(fp,"%d\n",&(params->max_iters));
  if (retval != 1) DIE("Could not read param file: max_iters");
  retval = fscanf(fp,"%d\n",&(params->reynolds_dim));
  if (retval != 1) DIE("Could not read param file: reynolds_dim");
  retval = fscanf(fp,"%f\n",&(params->density));
  if (retval != 1) DIE("Could not read param file: density");
  retval = fscanf(fp,"%f\n",&(params->accel));
  if (retval != 1) DIE("Could not read param file: accel");
  retval = fscanf(fp,"%f\n",&(params->omega));
  if (retval != 1) DIE("Could not read param file: omega");

  if (params->nx < 100) DIE("x dimension of grid in input file was too small (must be >100)");
  if (params->ny < 100) DIE("y dimension of grid in input file was too small (must be >100)");

  /* read column/row to accelerate */
  char accel_dir_buf[11];
  int idx;
  retval = fscanf(fp,"%*s %10s %d\n", accel_dir_buf, &idx);
  if (retval != 2) DIE("Could not read param file: could not parse acceleration specification");
  if (idx > 100 || idx < 0) DIE("Acceleration index (%d) out of range (must be bigger than 0 and less than 100)", idx);

  if (!(strcmp(accel_dir_buf, "row")))
  {
      accel_area->col_or_row = ACCEL_ROW;
      accel_area->idx = idx*(params->ny/BOX_Y_SIZE);
  }
  else if (!(strcmp(accel_dir_buf, "column")))
  {
      accel_area->col_or_row = ACCEL_COLUMN;
      accel_area->idx = idx*(params->nx/BOX_X_SIZE);
  }
  else
  {
      DIE("Error reading param file: Unexpected acceleration specification '%s'", accel_dir_buf);
  }


  /* read obstacles */
  retval = fscanf(fp, "%d %*s\n", &n_obstacles);
  if (retval != 1) DIE("Could not read param file: n_obstacles");
  obstacles = (obstacle_t*) malloc(sizeof(obstacle_t)*(n_obstacles));

  for (ii = 0; ii < n_obstacles; ii++)
  {
      retval = fscanf(fp,"%f %f %f %f\n",
          &obstacles[ii].obs_x_min, &obstacles[ii].obs_y_min,
          &obstacles[ii].obs_x_max, &obstacles[ii].obs_y_max);
      if (retval != 4) DIE("Could not read param file: location of obstacle %d", ii + 1);
      if (obstacles[ii].obs_x_min < 0 || obstacles[ii].obs_y_min < 0 ||
          obstacles[ii].obs_x_max > 100 || obstacles[ii].obs_y_max > 100)
      {
          DIE("Obstacle %d out of range (must be bigger than 0 and less than 100)", ii);
      }
      if (obstacles[ii].obs_x_min > obstacles[ii].obs_x_max) DIE("Left x coordinate is bigger than right x coordinate - this will result in no obstacle being made");
      if (obstacles[ii].obs_y_min > obstacles[ii].obs_y_max) DIE("Bottom y coordinate is bigger than top y coordinate - this will result in no obstacle being made");
  }

  /* close file */
  fclose(fp);

  // Allocate arrays
  char* tmp_obstacles = (char*) malloc(sizeof(char)*(params->ny*params->nx));
  if (tmp_obstacles == NULL) DIE("Cannot allocate memory for patches");

  *av_vels_ptr = (float*) malloc(sizeof(float)*(params->max_iters));
  if (*av_vels_ptr == NULL) DIE("Cannot allocate memory for av_vels");

  // Initialise obstacles
  for (ii = 0; ii < params->ny; ii++)
  {
    for (jj = 0; jj < params->nx; jj++)
    {
      (tmp_obstacles)[ii*params->nx + jj] = 0;
    }
  }

  /* Fill in locations of obstacles */
  for (ii = 0; ii < params->ny; ii++)
  {
      for (jj = 0; jj < params->nx; jj++)
      {
          /* coordinates of (jj, ii) scaled to 'real world' terms */
          const float x_pos = jj*(BOX_X_SIZE/params->nx);
          const float y_pos = ii*(BOX_Y_SIZE/params->ny);

          for (kk = 0; kk < n_obstacles; kk++)
          {
              if (x_pos >= obstacles[kk].obs_x_min &&
                  x_pos <  obstacles[kk].obs_x_max &&
                  y_pos >= obstacles[kk].obs_y_min &&
                  y_pos <  obstacles[kk].obs_y_max)
              {
                  (tmp_obstacles)[ii*params->nx + jj] = 1;
              }
          }
      }
  }

  // Redefine params by calculating the bounding box
  bounding_box(tmp_obstacles, params);

  // Allocate actual obstacles array of bounding box size, i.e. ignoring regions outside of bbox
  *obstacles_ptr = (char*) malloc(sizeof(char)*(params->ny*params->nx));
  if (*obstacles_ptr == NULL) DIE("Cannot allocate memory for patches");
  for(ii = 0; ii < params->ny; ii++)
  {
    for(jj = 0; jj < params->nx; jj++)
    {
      (*obstacles_ptr)[ii*params->nx + jj] = (tmp_obstacles)[(ii + params->obs_bbox.y1) * (params->original_grid.x2 + 1) + (jj + params->obs_bbox.x1)];
    }
  }

  // Translate the accelerated index for bbox coords
  if(accel_area->col_or_row == ACCEL_ROW) {
    accel_area->idx = accel_area->idx - (params->obs_bbox.y1);
  } else {
    accel_area->idx = accel_area->idx - (params->obs_bbox.x1);
  }

  free(obstacles);

  // Allocate size arrays
  params->loc_nxs = (int*)malloc(sizeof(int) * params->size);
  params->loc_nys = (int*)malloc(sizeof(int) * params->size);
}

void allocateLocal(param_t* params, speed_t** cells_ptr, speed_t** tmp_cells_ptr, speed_t** tmp_tmp_cells_ptr)
{
    int    ii,jj;          /* generic counters */
    float w0,w1,w2;       /* weighting factors */

    /* Allocate arrays, including halo space */
    *cells_ptr = (speed_t*) malloc(sizeof(speed_t)*(params->loc_ny * params->loc_nx));
    if (*cells_ptr == NULL) DIE("Cannot allocate memory for cells");

    *tmp_cells_ptr = (speed_t*) malloc(sizeof(speed_t)*(params->loc_ny * params->loc_nx));
    if (*tmp_cells_ptr == NULL) DIE("Cannot allocate memory for tmp_cells");

    *tmp_tmp_cells_ptr = (speed_t*) malloc(sizeof(speed_t)*(params->loc_ny * params->loc_nx));
    if (*tmp_tmp_cells_ptr == NULL) DIE("Cannot allocate memory for tmp_tmp_cells");

    w0 = params->density * 4.0/9.0;
    w1 = params->density      /9.0;
    w2 = params->density      /36.0;

    /* Initialise arrays */
    for (ii = 0; ii < params->loc_ny; ii++)
    {
        for (jj = 0; jj < params->loc_nx; jj++)
        {
            /* centre */
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[0] = w0;
            /* axis directions */
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[1] = w1;
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[2] = w1;
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[3] = w1;
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[4] = w1;
            /* diagonals */
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[5] = w2;
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[6] = w2;
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[7] = w2;
            (*cells_ptr)[ii*params->loc_nx + jj].speeds[8] = w2;
        }
    }

    // Allocate send/receive buffers
    params->sendbuf_d = (float*)malloc(sizeof(float) * params->loc_nx * 3);
    params->sendbuf_u = (float*)malloc(sizeof(float) * params->loc_nx * 3);
    params->recvbuf_d = (float*)malloc(sizeof(float) * params->loc_nx * 3);
    params->recvbuf_u = (float*)malloc(sizeof(float) * params->loc_nx * 3);
}

int get_global_y_coord(const param_t params, int rank, int ii) {
  int kk;
  int sum = 0;
  for(kk = 0; kk < rank; kk++) {
    // First read the size of the data to read back
    sum += params.loc_nys[kk];
  }
  return (sum + ii);
}

// Calculate the bounding box of the obstacles so that we dont do any
// unnecessary computation
void bounding_box(char* obstacles, param_t* params)
{
  int ii, jj;
  int min_x = params->nx - 1;
  int min_y = params->nx - 1;
  int max_x = 0;
  int max_y = 0;
  // find the edges of the box
  for (ii = 0; ii < params->ny; ii++)
  {
    for (jj = 0; jj < params->nx; jj++)
    {
      if(obstacles[ii*params->nx + jj] == 0) {
        if(jj < min_x) min_x = jj;
        if(jj > max_x) max_x = jj;
        if(ii < min_y) min_y = ii;
        if(ii > max_y) max_y = ii;
      }
    }
  }

  // then give bbox a border of 1 so that obstacles contain the region if necessary
  bbox_t box;
  min_x = (min_x - 1 < 0) ? 0 : min_x - 1;
  max_x = (max_x + 1 > params->nx - 1) ? params->nx - 1 : max_x + 1;
  min_y = (min_y - 1 < 0) ? 0 : min_y - 1;
  max_y = (max_y + 1 > params->ny - 1) ? params->ny - 1: max_y + 1;
  box.x1 = min_x;
  box.x2 = max_x;
  box.y1 = min_y;
  box.y2 = max_y;
  params->obs_bbox = box;

  // remember the original problem dimensions
  bbox_t orig;
  orig.x1 = 0;
  orig.x2 = params->nx - 1;
  orig.y1 = 0;
  orig.y2 = params->ny - 1;
  params->original_grid = orig;

  // set new size of problem
  // (plus one because of 0 based indexing)
  params->nx = (max_x - min_x) + 1;
  params->ny = (max_y - min_y) + 1;
}

void finalise(speed_t** cells_ptr, speed_t** tmp_cells_ptr, speed_t** tmp_tmp_cells_ptr,
    char** obstacles_ptr, float** av_vels_ptr)
{
    /* Free allocated memory */
    free(*cells_ptr);
    free(*tmp_cells_ptr);
    free(*tmp_tmp_cells_ptr);
    free(*obstacles_ptr);
    free(*av_vels_ptr);
}
