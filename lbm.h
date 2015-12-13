#ifndef LBM_HDR_FILE
#define LBM_HDR_FILE

#define MASTER 0
#define NSPEEDS         9

/* Size of box in imaginary 'units */
#define BOX_X_SIZE (100.0)
#define BOX_Y_SIZE (100.0)

/* struct to hold the 'speed' values */
typedef struct {
    float speeds[NSPEEDS];
} speed_t;

/* struct to hold the parameter values */
typedef struct {
    int nx;            /* no. of cells in x-direction */
    int ny;            /* no. of cells in y-direction */

    int loc_nx;        // no. of cells in local array in x-direction
    int loc_ny;        // no. of cells in local array in y-direction
    int* loc_nxs;      // array of rank sizes in x-dimension
    int* loc_nys;      // array of rank sizes in y-dimension
    int rank;          // 'rank' of process among it's cohort
    int size;          // size of cohort, i.e. num processes started
    int down;          // the rank of the process below
    int up;            // the rank of the process above
    float* sendbuf_d;
    float* sendbuf_u;
    float* recvbuf_d;
    float* recvbuf_u;

    int max_iters;      /* no. of iterations */
    int reynolds_dim;  /* dimension for Reynolds number */
    float density;       /* density per link */
    float accel;         /* density redistribution */
    float omega;         /* relaxation parameter */
} param_t;

/* obstacle positions */
typedef struct {
    float obs_x_min;
    float obs_x_max;
    float obs_y_min;
    float obs_y_max;
} obstacle_t;

typedef enum { ACCEL_ROW, ACCEL_COLUMN } accel_e;
typedef struct {
    accel_e col_or_row;
    int idx;
} accel_area_t;

/* Parse command line arguments to get filenames */
void parse_args (int argc, char* argv[],
    char** final_state_file, char** av_vels_file, char** param_file);

void initialise(const char* param_file, accel_area_t * accel_area, param_t* params, int** obstacles_ptr, float** av_vels_ptr);
void allocateLocal(param_t* params, speed_t** cells_ptr, speed_t** tmp_cells_ptr, speed_t** tmp_tmp_cells_ptr);
int get_global_y_coord(const param_t params, int rank, int ii);

int calc_nrows_from_rank(const param_t params, const int rank);

void write_values(const char * final_state_file, const char * av_vels_file,
    const param_t params, speed_t* cells, int* obstacles, float* av_vels);

void finalise(speed_t** cells_ptr, speed_t** tmp_cells_ptr, speed_t** tmp_tmp_cells_ptr,
    int** obstacles_ptr, float** av_vels_ptr);

float timestep(const param_t params, const accel_area_t accel_area,
    speed_t* cells, speed_t* tmp_cells, speed_t* tmp_tmp_cells, int* obstacles);
void accelerate_flow(const param_t params, const accel_area_t accel_area,
    speed_t* cells, int* obstacles);
float d2q9bgk(const param_t params, speed_t* cells, speed_t* tmp_cells, speed_t* tmp_tmp_cells, int* obstacles);

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const param_t params, speed_t* cells);

/* calculate Reynolds number */
float calc_reynolds(const param_t params, const float last_av_vel);

/* Exit, printing out formatted string */
#define DIE(...) exit_with_error(__LINE__, __FILE__, __VA_ARGS__)
void exit_with_error(int line, const char* filename, const char* format, ...)
__attribute__ ((format (printf, 3, 4)));

#endif
