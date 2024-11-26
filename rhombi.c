#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <errno.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include "mt19937.h"

/* Parameters */
#define pi 3.14159265358979323846
#define hr3 0.86602540378443864676
#define lambda 3.73205080756887729352
#define MAX_HAIR_LENGHT 10
#define LAT_SIZE_MAX 60000
#define NNLIST_SIZE 512 
#define NNDEPTH 6       // should (at least) match biggest n-gon in the tiles
#define MAXERR 100
#define STABLE_EPS_RH 1.0
#define STABLE_GAMMA 4.0
#define ANNEALING_SWEEPS 6000000
// #define PROFILING

#define beta 1.0
// #define D 1

// Simulation state variables
int st_I = 0;
int st_D = 10;

double eps_gamma = 3.0;
double eps_sq = 0.0;
double eps_tr = 0;
double eps_rh = -0.75;

// Simulation runtime variables
int MCsweeps = 24000000; // Total number of sweeps
int output_interval = 15000;        // writing interval
char datadir[256];
char measurementfilename[256];
int top = -1, old_top = -1, backup_top = -1;
int errcount = 0;

void myprintf(const char *format, ...)
{
#ifndef PROFILING
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    va_end(args);
#endif
}

void myfprintf(FILE *f, const char *format, ...)
{
#ifndef PROFILING
    va_list args;
    va_start(args, format);
    vfprintf(f, format, args);
    va_end(args);
#endif
}

void init_datadir(void)
{
    sprintf(datadir, "data_size_%d_%d_gamma_%lf_epsRh_%lf", st_D, st_I, eps_gamma, eps_rh);

    if (mkdir(datadir, S_IRWXU | S_IRWXG | S_IROTH) == -1)
    {
        myprintf("Error: %s\n", strerror(errno));
    }
}

void parse_cli(int argc, char *argv[])
{
    //    return;
    if (argc != 5)
    {
        myprintf("Usage:\n ./rhombi <st_D> <st_I> <gamma> <eps_rh>\n");
        exit(0);
    }
    st_D = atoi(argv[1]);
    st_I = atoi(argv[2]);
    eps_gamma = atof(argv[3]);
    eps_rh = atof(argv[4]);
}

/* Data structures */
typedef struct double2_s
{
    double x, y;
} double2;

typedef struct int2_s
{
    int x, y;
} int2;

typedef struct int4_s
{
    int n0, n1, n2, n3;
} int4;

typedef struct vertex_s
{
    int index;
    double2 par;
    double2 perp;
    int n0, n1, n2, n3;
    int bonds[12];

} vertex;

typedef struct empty_house_s
{
    /*
     *    tip
     *    / \
     *   /   \
     * lm-----rm
     *  |     |
     *  |     |
     *  lb---rb
     */
    int tip, lm, rm, lb, rb;
    int dir_lm_rm;
} empty_house;

typedef struct full_house_s
{
    /*
     *    tip
     *    /|\
     *   / | \
     * lm  o  rm
     *  | / \ |
     *  |/   \|
     *  lb---rb
     */
    int tip, lm, rm, lb, rb, o;
    int dir_o_tip;
} full_house;

/* Constants */
const double2 unit_par[4] = {{1, 0}, {hr3, .5}, {.5, hr3}, {0, 1}};
const double2 unit_perp[4] = {{1, 0}, {-hr3, -.5}, {.5, hr3}, {0, -1}};

const int4 dodecagon_vertices[18] = {
#include "dodecagon_vertices.txt"
};

const int4 unit_vecs[12] = {
#include "unit_vecs.txt"
};

/* Variables */

int N = 0;
int boundary_length = 0;
int new_boundary_length = 0;
int old_boundary_length = 0;
int backup_boundary_length = 0;

int n_sq, n_tr, n_rh;
int old_n_sq, old_n_tr, old_n_rh;
int backup_n_sq, backup_n_tr, backup_n_rh;

int sweep;
int shapes[3] = {0, 0, 0};        // squares, triangles, rhombi
int old_shapes[3] = {0, 0, 0};    // squares, triangles, rhombi
int backup_shapes[3] = {0, 0, 0}; // squares, triangles, rhombi

int shapes_1[3] = {0, 0, 0};
int old_shapes_1[3] = {0, 0, 0};
int backup_shapes_1[3] = {0, 0, 0};

int hist[12];

int countall = 0;

vertex lattice[LAT_SIZE_MAX];
vertex oldlattice[LAT_SIZE_MAX];
vertex backup_lattice[LAT_SIZE_MAX];

empty_house *empty_houses;
int n_empty_houses = 0;

full_house *full_houses;
int n_full_houses = 0;

// Given three collinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
int onSegment(double2 p, double2 q, double2 r)
{
    if (q.x <= fmax(p.x, r.x) && q.x >= fmin(p.x, r.x) &&
        q.y <= fmax(p.y, r.y) && q.y >= fmin(p.y, r.y))
        return 1;

    return 0;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(double2 p, double2 q, double2 r)
{
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0)
        return 0; // collinear

    return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
int doIntersect(double2 p1, double2 q1, double2 p2, double2 q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return 1;

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1))
        return 1;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1))
        return 1;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2))
        return 1;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2))
        return 1;

    return 0; // Doesn't fall in any of the above cases
}

void write_phase_state(void)
{
    char statefile[256];
    sprintf(statefile, "%s/state.csv", datadir);
    FILE *f = fopen(statefile, "w");
    myfprintf(f, "#N, D, I, gamma, eps_sq, eps_tr, eps_rh\n%d,%d,%d,%.10e,%.10e,%.10e,%.10e\n", N, st_D, st_I,
              eps_gamma, eps_sq, eps_tr, eps_rh);

    fclose(f);
}

void verify_lattice(void)
{
    for (int i = 0; i < N; i++)
    {
        for (int n = 0; n < 12; n++)
        {
            int j = lattice[i].bonds[n];
            if (j == -1)
            {
                continue;
            }
            assert(lattice[j].bonds[(n + 6) % 12] == i);
        }
    }
}

void write_list(int id, int list[NNLIST_SIZE], int len)
{
    char fname[256];
    sprintf(fname, "%s/nnlist_%09d.dat", datadir, id);
    FILE *f = fopen(fname, "w");
    for (int i = 0; i < len; i++)
    {
        int n = list[i];
        myfprintf(f, "%d %lf %lf\n", n, lattice[n].par.x, lattice[n].par.y);
    }
    fclose(f);
}

int is_member(int idx, int list[NNLIST_SIZE], int len)
{
    for (int i = 0; i < len; i++)
    {
        if (list[i] == idx)
        {
            return 1;
        }
    }
    return 0;
}

// Prints BFS traversal from a given source s
void Graph_BFS(int s, int max_dept, int visited[NNLIST_SIZE], int *len_visited)
{
    // Create a queue for BFS
    int level = 0;
    int queue[NNLIST_SIZE];
    uint8_t *visited_lut = calloc(N, sizeof(uint8_t));
    int front = 0, rear = 0;
    //    double dist = sqrt(dist2);

    // Mark the current node as visited and enqueue it
    assert(*len_visited < NNLIST_SIZE);
    visited[(*len_visited)++] = s;
    visited_lut[s] = 1;
    queue[rear++] = s;

    while (front != rear)
    {
        int level_size = rear - front;
        for (int _i = 0; _i < level_size; _i++)
        {
            // Dequeue a vertex from queue and print it
            int ss = queue[front++];

            // Get all adjacent vertices of the dequeued
            // vertex s.
            // If an adjacent has not been visited,
            // then mark it visited and enqueue it
            for (int adjacent = 0; adjacent < 12;
                 adjacent++)
            {
                int v = lattice[ss].bonds[adjacent];
                if (v == -1)
                {
                    continue;
                }

                int not_visited = 1 - visited_lut[v];

                if (not_visited)
                {
                    visited[(*len_visited)++] = v;
                    visited_lut[v] = 1;
                    if (rear < NNLIST_SIZE)
                    {
                        queue[rear++] = v;
                    }
                }

                if (*len_visited >= NNLIST_SIZE)
                {
                    free(visited_lut);
                    return;
                }
            }
        }
        level++;
        if (level == max_dept)
        {
            break;
        }
    }
    free(visited_lut);
}

void makeHist(void)
{
    for (int i = 0; i < 12; i++)
        hist[i] = 0;

    for (int n = 0; n < N; n++)
    {
        for (int i = 0; i < 12; i++)
        {
            if (lattice[n].bonds[i] > -1)
                hist[i]++;
        }
    }
}

void init_measurement_file(void)
{
    sprintf(measurementfilename, "%s/measurements.csv", datadir);
    FILE *f = fopen(measurementfilename, "w");

    myfprintf(f, "#MCSweep, n_square, n_triangle, n_rhombus, boundary_length, n_e1, n_e2, n_e3, n_e4, n_e5, n_e6, n_e7, n_e8, n_e9, n_e10, n_e11, n_e12\n");
    fclose(f);
}

void write_measurement(int sweep)
{
    FILE *f = fopen(measurementfilename, "a");

    myfprintf(f, "%d,%d,%d,%d,%d", sweep, n_sq, n_tr, n_rh, boundary_length);
    makeHist();
    for (int i = 0; i < 12; i++)
    {
        myfprintf(f, ",%d", hist[i]);
    }
    myfprintf(f, "\n");
    fclose(f);
}

void save_state(void)
{
    old_boundary_length = boundary_length;
    old_n_sq = n_sq;
    old_n_tr = n_tr;
    old_n_rh = n_rh;
    old_top = top;

    memcpy(old_shapes, shapes, 3 * sizeof(int));
    memcpy(old_shapes_1, shapes_1, 3 * sizeof(int));
    memcpy(oldlattice, lattice, N * sizeof(vertex));
}

void revert_state(void)
{
    boundary_length = old_boundary_length;
    n_sq = old_n_sq;
    n_tr = old_n_tr;
    n_rh = old_n_rh;
    top = old_top;

    memcpy(shapes, old_shapes, 3 * sizeof(int));
    memcpy(shapes_1, old_shapes_1, 3 * sizeof(int));
    memcpy(lattice, oldlattice, N * sizeof(vertex));
}

void revert_state_from_list(int list[NNLIST_SIZE], int len)
{
    boundary_length = old_boundary_length;
    n_sq = old_n_sq;
    n_tr = old_n_tr;
    n_rh = old_n_rh;
    top = old_top;

    memcpy(shapes, old_shapes, 3 * sizeof(int));
    memcpy(shapes_1, old_shapes_1, 3 * sizeof(int));

    for (int i = 0; i < len; i++)
    {
        int n = list[i];
        memcpy(&lattice[n], &oldlattice[n], sizeof(vertex));
    }
}

void save_backup_state(void)
{
    backup_boundary_length = boundary_length;
    backup_n_sq = n_sq;
    backup_n_tr = n_tr;
    backup_n_rh = n_rh;
    backup_top = top;

    memcpy(backup_shapes, shapes, 3 * sizeof(int));
    memcpy(backup_shapes_1, shapes_1, 3 * sizeof(int));
    memcpy(backup_lattice, lattice, N * sizeof(vertex));
}

void revert_backup_state(void)
{
    boundary_length = backup_boundary_length;
    n_sq = backup_n_sq;
    n_tr = backup_n_tr;
    n_rh = backup_n_rh;
    top = backup_top;

    memcpy(shapes, backup_shapes, 3 * sizeof(int));
    memcpy(shapes_1, backup_shapes_1, 3 * sizeof(int));
    memcpy(lattice, backup_lattice, N * sizeof(vertex));
}

void calculateParCoords(vertex *v)
{
    v->par.x = v->n0 * unit_par[0].x + v->n1 * unit_par[1].x + v->n2 * unit_par[2].x + v->n3 * unit_par[3].x;
    v->par.y = v->n0 * unit_par[0].y + v->n1 * unit_par[1].y + v->n2 * unit_par[2].y + v->n3 * unit_par[3].y;
    return;
}

void calculatePerpCoords(vertex *v)
{
    v->perp.x = v->n0 * unit_perp[0].x + v->n1 * unit_perp[1].x + v->n2 * unit_perp[2].x + v->n3 * unit_perp[3].x;
    v->perp.y = v->n0 * unit_perp[0].y + v->n1 * unit_perp[1].y + v->n2 * unit_perp[2].y + v->n3 * unit_perp[3].y;
    return;
}

/* Functions */
void writeLattice(int step)
{

    char datafile1[128];

    // For saving steps seperately
    sprintf(datafile1, "%s/lattice%09d.dat", datadir, step);
    FILE *data1 = fopen(datafile1, "w");

    // // Format:
    // (1)index (2)par-space coords (2)perp-space coords (4)4-space coords (12)neighbour indices
    for (int n = 0; n < N; n++)
    {
        calculatePerpCoords(&lattice[n]);
        myfprintf(data1, "%d %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", n, lattice[n].par.x, lattice[n].par.y, lattice[n].perp.x, lattice[n].perp.y, lattice[n].n0, lattice[n].n1, lattice[n].n2, lattice[n].n3, lattice[n].bonds[0], lattice[n].bonds[1], lattice[n].bonds[2], lattice[n].bonds[3], lattice[n].bonds[4], lattice[n].bonds[5], lattice[n].bonds[6], lattice[n].bonds[7], lattice[n].bonds[8], lattice[n].bonds[9], lattice[n].bonds[10], lattice[n].bonds[11]);
    }

    // fclose(data);
    fclose(data1);
}
void writeoldLattice(int step)
{

    char datafile1[128];

    // For saving steps seperately
    sprintf(datafile1, "%s/lattice%09d.dat", datadir, step);
    FILE *data1 = fopen(datafile1, "w");

    // // Format:
    // (1)index (2)par-space coords (2)perp-space coords (4)4-space coords (12)neighbour indices
    for (int n = 0; n < N; n++)
    {
        calculatePerpCoords(&oldlattice[n]);
        myfprintf(data1, "%d %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", n, oldlattice[n].par.x, oldlattice[n].par.y, oldlattice[n].perp.x, oldlattice[n].perp.y, oldlattice[n].n0, oldlattice[n].n1, oldlattice[n].n2, oldlattice[n].n3, oldlattice[n].bonds[0], oldlattice[n].bonds[1], oldlattice[n].bonds[2], oldlattice[n].bonds[3], oldlattice[n].bonds[4], oldlattice[n].bonds[5], oldlattice[n].bonds[6], oldlattice[n].bonds[7], oldlattice[n].bonds[8], oldlattice[n].bonds[9], oldlattice[n].bonds[10], oldlattice[n].bonds[11]);
    }

    // fclose(data);
    fclose(data1);
}

void genLattice(void)
{
    // Generate a D x D square grid of vertices

    for (int n = 0; n < st_D * st_D; n++)
    {
        int i = n % st_D;
        int j = n / st_D;

        lattice[n].n0 = i;
        lattice[n].n1 = 0;
        lattice[n].n2 = 0;
        lattice[n].n3 = j;

        lattice[n].index = n;
    }

    N += st_D * st_D;
    return;
}

void genChevronLattice(void)
{

    for (int i = 0; i < st_D; i++)
    {
        for (int j = 0; j < st_D; j++)
        {
            int n = i * st_D + j;

            lattice[2 * n].n0 = i;
            lattice[2 * n].n1 = 0;
            lattice[2 * n].n2 = 0;
            lattice[2 * n].n3 = j;

            lattice[2 * n].index = 2 * n;

            lattice[2 * n + 1].n0 = i;
            lattice[2 * n + 1].n1 = 1;
            lattice[2 * n + 1].n2 = 0;
            lattice[2 * n + 1].n3 = j;

            lattice[2 * n + 1].index = 2 * n + 1;
        }
    }

    N = 2 * st_D * st_D;
    return;
}

void genDodecagon1(vertex center)
{
    // This generates a dodecagonal wheel of vertices around the vertex labelled 'center'

    int exists;

    for (int i = 0; i < 18; i++)
    {
        exists = 0;

        // Check if vertex already exists
        // D12vertices gives dodecagon vertex positions relative to center vertex
        for (int j = 0; j < N; j++)
        {
            if (lattice[j].n0 == center.n0 + dodecagon_vertices[i].n0 && lattice[j].n1 == center.n1 + dodecagon_vertices[i].n1)
            {
                if (lattice[j].n2 == center.n2 + dodecagon_vertices[i].n2 && lattice[j].n3 == center.n3 + dodecagon_vertices[i].n3)
                {
                    exists = 1;
                    break;
                }
            }
        }
        // Don't add it if it already exists
        // (Due to neighbouring dodecagon)
        if (exists == 1)
            continue;

        // Add new vertex
        lattice[N].n0 = center.n0 + dodecagon_vertices[i].n0;
        lattice[N].n1 = center.n1 + dodecagon_vertices[i].n1;
        lattice[N].n2 = center.n2 + dodecagon_vertices[i].n2;
        lattice[N].n3 = center.n3 + dodecagon_vertices[i].n3;
        lattice[N].index = N;
        N++;
    }
    return;
}

void genDodecagon2(vertex center)
{
    // This generates a dodecagonal wheel of vertices around the vertex labelled 'center'

    int exists;

    for (int i = 0; i < 18; i++)
    {
        exists = 0;

        // Check if vertex already exists
        // D12vertices gives dodecagon vertex positions relative to center vertex
        for (int j = 0; j < N; j++)
        {
            if (lattice[j].n0 == center.n0 + dodecagon_vertices[i].n3 && lattice[j].n1 == center.n1 + dodecagon_vertices[i].n2)
            {
                if (lattice[j].n2 == center.n2 + dodecagon_vertices[i].n1 && lattice[j].n3 == center.n3 + dodecagon_vertices[i].n0)
                {
                    exists = 1;
                    break;
                }
            }
        }
        // Don't add it if it already exists
        // (Due to neighbouring dodecagon)
        if (exists == 1)
            continue;

        // Add new vertex
        lattice[N].n0 = center.n0 + dodecagon_vertices[i].n3;
        lattice[N].n1 = center.n1 + dodecagon_vertices[i].n2;
        lattice[N].n2 = center.n2 + dodecagon_vertices[i].n1;
        lattice[N].n3 = center.n3 + dodecagon_vertices[i].n0;
        lattice[N].index = N;
        N++;
    }
    return;
}

void inflateLattice(void)
{
    // Performs one step of Stampfli inflation on the existing lattice

    int n0, n1, n2, n3;
    int max = N;

    // Inflate existing vertices
    for (int i = 0; i < max; i++)
    {
        n0 = 2 * lattice[i].n0 + lattice[i].n1 - lattice[i].n3;
        n1 = 2 * lattice[i].n0 + 2 * lattice[i].n1 + lattice[i].n2;
        n2 = lattice[i].n1 + 2 * lattice[i].n2 + 2 * lattice[i].n3;
        n3 = -lattice[i].n0 + lattice[i].n2 + 2 * lattice[i].n3;

        lattice[i].n0 = n0;
        lattice[i].n1 = n1;
        lattice[i].n2 = n2;
        lattice[i].n3 = n3;
    }

    // Replace inflated vertices with dodecagons
    for (int i = 0; i < max; i++)
    {
        if (dsfmt_genrand() < 0.5)
            genDodecagon1(lattice[i]);
        else
            genDodecagon2(lattice[i]);
    }

    return;
}

void removeBondCross(vertex *v)
{
    int n1, n2;
    for (int n = 0; n < 12; n++)
    {
        if ((n1 = v->bonds[n]) != -1)
        {

            if ((n2 = v->bonds[(n + 2) % 12]) != -1 && lattice[n1].bonds[(n + 5) % 12] != -1)
            {
                v->bonds[(n + 2) % 12] = -1;
                lattice[n2].bonds[(n + 8) % 12] = -1;
            }
            if ((n2 = v->bonds[(n + 10) % 12]) != -1 && lattice[n1].bonds[(n + 7) % 12] != -1)
            {
                v->bonds[(n + 10) % 12] = -1;
                lattice[n2].bonds[(n + 4) % 12] = -1;
            }
        }
    }

    return;
}
int initNeighbours(vertex *vertex_i)
{
    /* Finds the neighbours of 'vertex_i'
       and stores their indices in the '.bonds list' */

    int i;
    int n0i, n1i, n2i, n3i;
    int n0j, n1j, n2j, n3j;
    int equal, diff1, diff2;

    vertex *vertex_j;

    // Reset neighbour list
    for (int n = 0; n < 12; n++)
        vertex_i->bonds[n] = -1;

    n0i = vertex_i->n0;
    n1i = vertex_i->n1;
    n2i = vertex_i->n2;
    n3i = vertex_i->n3;
    for (int j = 0; j < N; j++)
    {
        i = vertex_i->index;
        vertex_j = &lattice[j];

        if (i == j)
            continue;

        n0j = vertex_j->n0;
        n1j = vertex_j->n1;
        n2j = vertex_j->n2;
        n3j = vertex_j->n3;

        equal = 0;
        if (n0i == n0j)
            equal += 1;
        if (n1i == n1j)
            equal += 10;
        if (n2i == n2j)
            equal += 100;
        if (n3i == n3j)
            equal += 1000;

        if (equal < 101)
            continue;

        if (equal == 1111)
            return 0;

        if (equal == 1110)
        {
            // n0 is different
            diff1 = n0j - n0i;
            if (diff1 == 1)
            {
                vertex_i->bonds[0] = j;
                vertex_j->bonds[6] = i;
            }
            else if (diff1 == -1)
            {
                vertex_i->bonds[6] = j;
                vertex_j->bonds[0] = i;
            }
        }
        else if (equal == 1101)
        {
            // n1 is different
            diff1 = n1j - n1i;
            if (diff1 == 1)
            {
                vertex_i->bonds[1] = j;
                vertex_j->bonds[7] = i;
            }
            else if (diff1 == -1)
            {
                vertex_i->bonds[7] = j;
                vertex_j->bonds[1] = i;
            }
        }
        else if (equal == 1011)
        {
            // n2 is different
            diff1 = n2j - n2i;
            if (diff1 == 1)
            {
                vertex_i->bonds[2] = j;
                vertex_j->bonds[8] = i;
            }
            else if (diff1 == -1)
            {
                vertex_i->bonds[8] = j;
                vertex_j->bonds[2] = i;
            }
        }
        else if (equal == 111)
        {
            // n3 is different
            diff1 = n3j - n3i;
            if (diff1 == 1)
            {
                vertex_i->bonds[3] = j;
                vertex_j->bonds[9] = i;
            }
            else if (diff1 == -1)
            {
                vertex_i->bonds[9] = j;
                vertex_j->bonds[3] = i;
            }
        }
        else if (equal == 1010)
        {
            // n0 and n2 are different
            diff1 = n0j - n0i;
            diff2 = n2j - n2i;
            if (diff1 == -1 && diff2 == 1)
            {
                vertex_i->bonds[4] = j;
                vertex_j->bonds[10] = i;
            }
            else if (diff1 == 1 && diff2 == -1)
            {
                vertex_i->bonds[10] = j;
                vertex_j->bonds[4] = i;
            }
        }
        else if (equal == 101)
        {
            // n1 and n3 are different
            diff1 = n1j - n1i;
            diff2 = n3j - n3i;
            if (diff1 == -1 && diff2 == 1)
            {
                vertex_i->bonds[5] = j;
                vertex_j->bonds[11] = i;
            }
            else if (diff1 == 1 && diff2 == -1)
            {
                vertex_i->bonds[11] = j;
                vertex_j->bonds[5] = i;
            }
        }
    }

    return 0;
}

int findNeighbours_from_list(vertex *vertex_i, int potential_nbs[NNLIST_SIZE], int len)
{
    /* Finds the neighbours of 'vertex_i'
       and stores their indices in the '.bonds list' */

    int i;
    int n0i, n1i, n2i, n3i;
    int n0j, n1j, n2j, n3j;
    int equal, diff1, diff2;

    vertex *vertex_j;

    // Reset neighbour list
    for (int n = 0; n < 12; n++)
    {
        if (vertex_i->bonds[n] != -1 && is_member(vertex_i->bonds[n], potential_nbs, len))
        {
            lattice[vertex_i->bonds[n]].bonds[(n + 6) % 12] = -1;
            vertex_i->bonds[n] = -1;
        }
    }
    n0i = vertex_i->n0;
    n1i = vertex_i->n1;
    n2i = vertex_i->n2;
    n3i = vertex_i->n3;
    //    for (int j = 0; j < N; j++)
    for (int jj = 0; jj < len; jj++)
    {
        int j = potential_nbs[jj];
        i = vertex_i->index;

        if (i == j)
            continue;

        vertex_j = &lattice[j];

        n0j = vertex_j->n0;
        n1j = vertex_j->n1;
        n2j = vertex_j->n2;
        n3j = vertex_j->n3;

        equal = 0;
        if (n0i == n0j)
            equal += 1;
        if (n1i == n1j)
            equal += 10;
        if (n2i == n2j)
            equal += 100;
        if (n3i == n3j)
            equal += 1000;

        if (equal < 101)
            continue;

        if (equal == 1111)
            return 0;

        if (equal == 1110)
        {
            // n0 is different
            diff1 = n0j - n0i;
            if (diff1 == 1)
            {
                vertex_i->bonds[0] = j;
                vertex_j->bonds[6] = i;
            }
            else if (diff1 == -1)
            {
                vertex_i->bonds[6] = j;
                vertex_j->bonds[0] = i;
            }
        }
        else if (equal == 1101)
        {
            // n1 is different
            diff1 = n1j - n1i;
            if (diff1 == 1)
            {
                vertex_i->bonds[1] = j;
                vertex_j->bonds[7] = i;
            }
            else if (diff1 == -1)
            {
                vertex_i->bonds[7] = j;
                vertex_j->bonds[1] = i;
            }
        }
        else if (equal == 1011)
        {
            // n2 is different
            diff1 = n2j - n2i;
            if (diff1 == 1)
            {
                vertex_i->bonds[2] = j;
                vertex_j->bonds[8] = i;
            }
            else if (diff1 == -1)
            {
                vertex_i->bonds[8] = j;
                vertex_j->bonds[2] = i;
            }
        }
        else if (equal == 111)
        {
            // n3 is different
            diff1 = n3j - n3i;
            if (diff1 == 1)
            {
                vertex_i->bonds[3] = j;
                vertex_j->bonds[9] = i;
            }
            else if (diff1 == -1)
            {
                vertex_i->bonds[9] = j;
                vertex_j->bonds[3] = i;
            }
        }
        else if (equal == 1010)
        {
            // n0 and n2 are different
            diff1 = n0j - n0i;
            diff2 = n2j - n2i;
            if (diff1 == -1 && diff2 == 1)
            {
                vertex_i->bonds[4] = j;
                vertex_j->bonds[10] = i;
            }
            else if (diff1 == 1 && diff2 == -1)
            {
                vertex_i->bonds[10] = j;
                vertex_j->bonds[4] = i;
            }
        }
        else if (equal == 101)
        {
            // n1 and n3 are different
            diff1 = n1j - n1i;
            diff2 = n3j - n3i;
            if (diff1 == -1 && diff2 == 1)
            {
                vertex_i->bonds[5] = j;
                vertex_j->bonds[11] = i;
            }
            else if (diff1 == 1 && diff2 == -1)
            {
                vertex_i->bonds[11] = j;
                vertex_j->bonds[5] = i;
            }
        }
    }

    return 0;
}

void countShapes(vertex v)
{
    // Counts shapes
    int start, length, neighbour, direction = -1, prev_dir;
    ;
    vertex current;
    int first_visited = -1, last_visited = -1;

    shapes_1[0] = 0;
    shapes_1[1] = 0;
    shapes_1[2] = 0;

    for (int tile = 0; tile < 12; tile++)
    {
        length = 0;
        current = v;
        if (current.bonds[tile] == -1)
            continue;
        start = tile;
        prev_dir = start;
        for (int segment = 0; segment < 5; segment++)
        {
            for (direction = start; direction < start + 5; direction++)
            {
                neighbour = current.bonds[direction % 12];
                if (neighbour == -1)
                    continue;

                if (length == 0)
                {
                    first_visited = neighbour;
                    last_visited = neighbour;
                }
                else
                {
                    last_visited = current.index;
                }

                current = lattice[neighbour];
                length++;

                start = (direction + 7) % 12;
                break;
            }
            if (current.index == v.index)
            {
                break;
            }
            prev_dir = direction;
        }
        if (direction == -1)
        {
            continue;
        }
        if (first_visited == -1 || last_visited == -1)
        {
            continue;
        }
        if (first_visited == last_visited)
        {
            continue;
        }
        if (current.index != v.index)
        {
            continue;
        }
        if (length == 4)
        {
            int angle = direction - prev_dir;
            if (angle == -3 || angle == 9)
            {
                shapes_1[0]++;
            }
            if (angle == -5 || angle == 7)
            {
                shapes_1[2]++;
            }
            if (angle == -1 || angle == 11)
            {
                shapes_1[2]++;
            }
        }
        if (length == 3)
        {
            int angle = direction - prev_dir;
            if (angle == -4 || angle == 8)
            {
                shapes_1[1]++;
            }
        }
    }
    return;
}

void countAllShapes(void)
{
    countall++;
    // Counts shapes
    int start, length, neighbour, direction = -1, prev_dir;
    ;
    vertex current, v;

    shapes[0] = 0;
    shapes[1] = 0;
    shapes[2] = 0;

    for (int i = 0; i < N; i++)
    {
        v = lattice[i];
        int first_visited = -1, last_visited = -1;

        for (int tile = 0; tile < 12; tile++)
        {
            length = 0;
            current = v;
            if (current.bonds[tile] == -1)
                continue;
            start = tile;
            prev_dir = start;
            for (int segment = 0; segment < 100; segment++)
            {
                for (direction = start; direction < start + 5; direction++)
                {
                    neighbour = current.bonds[direction % 12];
                    if (neighbour == -1)
                        continue;
                    if (length == 0)
                    {
                        first_visited = neighbour;
                        last_visited = neighbour;
                    }
                    else
                    {
                        last_visited = current.index;
                    }
                    current = lattice[neighbour];
                    length++;

                    start = (direction + 7) % 12;
                    break;
                }
                if (current.index == v.index)
                {
                    break;
                }
                prev_dir = direction;
            }
            if (direction == -1)
            {
                continue;
            }
            if (first_visited == -1 || last_visited == -1)
            {
                continue;
            }
            if (first_visited == last_visited)
            {
                continue;
            }
            if (current.index != v.index)
            {
                continue;
                ;
            }

            if (length == 4)
            {
                int angle = direction - prev_dir;
                if (angle == -3 || angle == 9)
                {
                    shapes[0]++;
                }
                if (angle == -5 || angle == 7)
                {
                    shapes[2]++;
                }
                if (angle == -1 || angle == 11)
                {
                    shapes[2]++;
                }
            }
            if (length == 3)
            {
                int angle = direction - prev_dir;
                if (angle == -4 || angle == 8)
                {
                    shapes[1]++;
                }
            }
        }
    }

    assert(shapes[0] % 4 == 0);
    shapes[0] = shapes[0] / 4;
    assert(shapes[1] % 3 == 0);
    shapes[1] = shapes[1] / 3;
    assert(shapes[2] % 4 == 0);
    shapes[2] = shapes[2] / 4;

    return;
}

int countShapes_from_list(int list[NNLIST_SIZE], int len)
{

    // Counts shapes
    int start, length, neighbour, direction = -1, prev_dir;
    ;
    vertex current, v;

    shapes_1[0] = 0;
    shapes_1[1] = 0;
    shapes_1[2] = 0;

    for (int ii = 0; ii < len; ii++) // loop over every vertex in the list
    {
        int i = list[ii]; // get the index of the vertex
        v = lattice[i];   // v is the vertex we are interested in
        int first_visited = -1, last_visited = -1;

        for (int tile = 0; tile < 12; tile++) // a tile can start in 12 different directions
        {
            length = 0;                    // initialize length of the circumference
            current = v;                   // initialize the current vertex
            if (current.bonds[tile] == -1) // if there is no bond in tile direction then there is no tile starting there
                continue;
            if (is_member(current.bonds[tile], list, len) == 0)
            { // if there  is a bond but it is not a vertex in the list then we are not interested
                continue;
            }
            // if we advenced to here then we have found our first vertex to jump to
            start = tile;
            prev_dir = start;
            for (int segment = 0; segment < 5; segment++) // now we will start hopping
            {
                for (direction = start; direction < start + 5; direction++)
                {
                    neighbour = current.bonds[direction % 12];
                    if (neighbour == -1)
                        continue;
                    if (is_member(neighbour, list, len) == 0)
                    {
                        continue;
                    }

                    if (length == 0)
                    {
                        first_visited = neighbour;
                        last_visited = neighbour;
                    }
                    else
                    {
                        last_visited = current.index;
                    }
                    current = lattice[neighbour];
                    length++;

                    start = (direction + 7) % 12;
                    break;
                }
                if (current.index == v.index)
                {
                    break;
                }
                prev_dir = direction;
            }
            if (direction == -1)
            {
                continue;
            }
            if (first_visited == -1 || last_visited == -1)
            {
                continue;
            }
            if (current.index != v.index)
            {
                continue;
            }
            if (first_visited == last_visited)
            {
                continue;
            }
            if (length == 4)
            {
                int angle = direction - prev_dir;
                if (angle == -3 || angle == 9)
                {
                    shapes_1[0]++;
                }
                if (angle == -5 || angle == 7)
                {
                    shapes_1[2]++;
                }
                if (angle == -1 || angle == 11)
                {
                    shapes_1[2]++;
                }
            }
            if (length == 3)
            {
                int angle = direction - prev_dir;
                if (angle == -4 || angle == 8)
                {
                    shapes_1[1]++;
                }
            }
        }
    }
    if (shapes_1[0] % 4 || shapes_1[1] % 3 || shapes_1[2] % 4)
    {
        return -1;
    }
    assert(shapes_1[0] % 4 == 0);
    shapes_1[0] = shapes_1[0] / 4;
    assert(shapes_1[1] % 3 == 0);
    shapes_1[1] = shapes_1[1] / 3;
    assert(shapes_1[2] % 4 == 0);
    shapes_1[2] = shapes_1[2] / 4;

    return 0;
}

int countNeighbours(vertex *v)
{
    // Counts neighbours
    int count = 0;
    for (int j = 0; j < 12; j++)
    {
        if (v->bonds[j] > -1)
            count++;
    }
    return count;
}

int checkAdjTiles(vertex v, int neighbours[12])
{
    /* This checks whether or not the tiles adjacent to vertex v are
       valid. (squares, triangles, rhombi, perimeter) */

    int start, length, neighbour;
    int first_visited = -1, last_visited = -1;
    int possible_boundary_length = 0;
    vertex current;

    for (int tile = 0; tile < 12; tile++)
    {
        int out_angles = -5;
        length = 0;
        current = v;
        if (current.bonds[tile] == -1)
            continue;
        start = tile;

        for (int segment = 0; segment < N; segment++)
        {
            for (int direction = start; direction < start + 12; direction++)
            {
                neighbour = current.bonds[direction % 12];
                if (neighbour == -1)
                    continue;

                if (length == 0)
                {
                    first_visited = neighbour;
                    last_visited = neighbour;
                }
                else
                {
                    last_visited = current.index;
                }

                current = lattice[neighbour];
                last_visited = neighbour;

                length++;

                out_angles += 6 - (direction - start + 1);

                // Check if we visited an old neighbour
                for (int nb = 0; nb < 12; nb++)
                {
                    if (current.index == lattice[neighbours[nb]].index)
                    {
                        neighbours[nb] = -1;
                    }
                }

                start = (direction + 7) % 12;

                break;
            }
            if (current.index == v.index)
            {
                out_angles += 5 - ((tile - start + 12) % 12);
                break;
            }
        }
        if (first_visited == -1 || last_visited == -1)
        {
            continue;
        }
        if (current.index != v.index)
        {
            continue;
        }
        if (first_visited == last_visited)
        {
            possible_boundary_length += length;
            continue;
        }
        if (length == new_boundary_length)
            continue;
        if (length == 3)
        {
            continue;
        }
        if (length == 4 && out_angles == 12)
        {
            continue;
        }

        return 1;
    }
    if (possible_boundary_length != 0 && possible_boundary_length != new_boundary_length)
    {
        return 1;
    }
    // Reject if not all old neighbours have been visited
    for (int i = 0; i < 12; i++)
    {
        if (neighbours[i] != -1)
        {
            // Check if unvisited vertex shares a neighbour with v
            for (int j = 0; j < 12; j++)
            {
                int k = v.bonds[j];
                if (k != -1)
                {
                    int notvalid = 1;
                    for (int l = 0; l < 12; l++)
                    {
                        if (lattice[k].bonds[l] == neighbours[i])
                        {
                            notvalid = 0;
                            break;
                        }
                    }
                    if (notvalid)
                    {
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

void update_top_vertex(int i)
{
    if (lattice[i].par.y > lattice[top].par.y)
        top = i;
}

void init_top_vertex(void)
{
    top = 0;
    for (int i = 0; i < N; i++)
    {
        update_top_vertex(i);
    }
}

int countBoundaryLength(void)
{
    int length, start, neighbour, start_dir = -1, prev_dir;
    vertex current;

    // Loop around until we find it again
    start = 6;
    prev_dir = start;
    length = 0;
    current = lattice[top];
    for (int step = 0; step < 50000; step++)
    {
        for (int direction = start; direction < start + 12; direction++)
        {
            neighbour = current.bonds[direction % 12];
            if (neighbour == -1)
                continue;
            if (length == 0)
                start_dir = direction;

            current = lattice[neighbour];
            length++;
            if (current.index == top)
            {
                assert(start_dir != -1);
                if (countNeighbours(&lattice[top]) == 1)
                    return length;
                else if (direction != (start_dir + 6) % 12 && direction != (prev_dir + 6) % 12)
                    return length;
            }
            prev_dir = direction;
            start = (direction + 7) % 12;
            break;
        }
    }
    return -1;
}

int checkDistances(vertex v)
{
    /* Reject configurations with distances < 1
       Except for the distance of 0.5 in a small rhombus */
    double dx, dy, dist2;
    for (int i = 0; i < N; i++)
    {
        if (v.index == lattice[i].index)
            continue;
        dx = v.par.x - lattice[i].par.x;
        dy = v.par.y - lattice[i].par.y;
        dist2 = dx * dx + dy * dy;
        if (dist2 < .9 && (dist2 > 0.27 || dist2 < 0.23))
        {
            return 1;
        }
    }
    return 0;
}

int checkDistances_from_list(vertex v, int list[NNLIST_SIZE], int len)
{
    /* Reject configurations with distances < 1
       Except for the distance of 0.5 in a small rhombus */
    double dx, dy, dist2;
    for (int j = 0; j < len; j++)
    {
        int i = list[j];
        if (v.index == lattice[i].index)
            continue;
        dx = v.par.x - lattice[i].par.x;
        dy = v.par.y - lattice[i].par.y;
        dist2 = dx * dx + dy * dy;

        if (dist2 < .999 && (fabs(dist2 - 0.267949) > 0.001))
        {
            return 1;
        }
    }
    return 0;
}

int notValid(vertex *v, int potential_nbs[NNLIST_SIZE], int len /*, int nb, int mv*/)
{
    /* This performs all the checks required to verify
       that a vertex move should be allowed */
    int old_neighbours[12];
    // Store neighbours
    for (int n = 0; n < 12; n++)
        old_neighbours[n] = v->bonds[n];

    // Update neighbours
    for (int n = 0; n < 12; n++)
    {
        int nb = old_neighbours[n];
        if (nb != -1)
        {
            findNeighbours_from_list(&lattice[nb], potential_nbs, len);
        }
    }
    //    findNeighbours(v);
    findNeighbours_from_list(v, potential_nbs, len);
    for (int n = 0; n < 12; n++)
    {
        int nb = v->bonds[n];
        if (nb != -1)
            findNeighbours_from_list(&lattice[nb], potential_nbs, len);
    }

    for (int n = 0; n < 12; n++)
    {
        int nb = old_neighbours[n];
        if (nb != -1)
            removeBondCross(&lattice[nb]);
        nb = v->bonds[n];
        if (nb != -1)
            removeBondCross(&lattice[nb]);
    }
    removeBondCross(v);

    // Check if any vertices have disconnected
    int nbc;
    for (int n = 0; n < 12; n++)
    {
        if (old_neighbours[n] != -1)
        {
            nbc = countNeighbours(&lattice[old_neighbours[n]]);
            if (nbc < 1)
            {
                return 1;
            }
            else if (nbc == 1)
            {
                for (size_t i = 0; i < 12; i++)
                {
                    if (lattice[old_neighbours[n]].bonds[i] != -1)
                    {
                        if (countNeighbours(&lattice[lattice[old_neighbours[n]].bonds[i]]) == 1)
                        {
                            // we have a detached line segment! reject!
                            return 1;
                        }
                        break;
                    }
                }
            }
        }
    }

    nbc = countNeighbours(v);
    if (nbc < 1)
    {
        return 1;
    }
    else if (nbc == 1)
    {
        for (size_t i = 0; i < 12; i++)
        {
            if (v->bonds[i] != -1)
            {
                if (countNeighbours(&lattice[v->bonds[i]]) == 1)
                {
                    // we have a detached line segment! reject!
                    return 1;
                }
                break;
            }
        }
    }

    // Check if any vertices are too close together
    if (checkDistances_from_list(*v, potential_nbs, len))
    {
        return 1;
    }

    // Check if the adjacent tiles are valid
    new_boundary_length = countBoundaryLength();
    if (checkAdjTiles(*v, old_neighbours))
    {
        return 1;
    }

    return 0;
}

int is_boundary(vertex *v)
{
    countShapes(*v);
    int n_edges = countNeighbours(v);
    int n_sq_old = shapes_1[0];
    int n_tr_old = shapes_1[1];
    int n_rh_old = shapes_1[2];

    int n_faces = n_sq_old + n_tr_old + n_rh_old;

    int d_EulerCharacteristic = (-1) - (-n_edges) + (-n_faces + 1);

    return d_EulerCharacteristic != 0;
}

int moveVertex(void)
{
    int random_index, neighbour_direction, move_direction;
    int dn_sq, dn_tr, dn_rh, n_sq_old, n_tr_old, n_rh_old;
    int reject = 0;
    double acc, energy_cost;
    vertex *chosen;
    vertex old, neighbour;

    // Pick random vertex
    random_index = (int)floor(N * dsfmt_genrand());
    chosen = &lattice[random_index];

    // Pick random direction
    neighbour_direction = (int)floor(12 * dsfmt_genrand());

    // Reject if there is no neighour there
    if (chosen->bonds[neighbour_direction] == -1)
    {
        return 0;
    }
    else
        neighbour = lattice[chosen->bonds[neighbour_direction]];

    // Rotate vertex around neighbour (2 possibilities)
    old = *chosen;
    if (dsfmt_genrand() < .5)
        move_direction = (neighbour_direction + 5) % 12;
    else
        move_direction = (neighbour_direction + 7) % 12;

    if (neighbour.bonds[move_direction] != -1)
    {
        // Overlap
        return 0;
    }

    int n_edges = countNeighbours(chosen);

    int potential_nbs[NNLIST_SIZE];
    int len = 0;
    int depth = NNDEPTH;
    Graph_BFS(chosen->index, depth, potential_nbs, &len);
    // Tiles around chosen vertex
    countShapes(*chosen);
    n_sq_old = shapes_1[0];
    n_tr_old = shapes_1[1];
    n_rh_old = shapes_1[2];

    int n_faces = n_sq_old + n_tr_old + n_rh_old;

    int d_EulerCharacteristic = (-1) - (-n_edges) + (-n_faces + 1);
    if (d_EulerCharacteristic != 0) // this is a boundary vertex
    {
        int cerr = countShapes_from_list(potential_nbs, len);
        assert(cerr == 0);

        n_sq_old = shapes_1[0];
        n_tr_old = shapes_1[1];
        n_rh_old = shapes_1[2];
    }

    // Move vertex
    chosen->n0 += unit_vecs[neighbour_direction].n0 + unit_vecs[move_direction].n0;
    chosen->n1 += unit_vecs[neighbour_direction].n1 + unit_vecs[move_direction].n1;
    chosen->n2 += unit_vecs[neighbour_direction].n2 + unit_vecs[move_direction].n2;
    chosen->n3 += unit_vecs[neighbour_direction].n3 + unit_vecs[move_direction].n3;

    calculateParCoords(chosen);
    update_top_vertex(chosen->index);

    // Check if valid
    if (notValid(chosen, potential_nbs, len))
    {
        revert_state_from_list(potential_nbs, len);
        return 0;
    }

    // New tiles around chosen vertex

    if (d_EulerCharacteristic == 0)
    {
        countShapes(*chosen);
    }
    else // this is a boundary vertex
    {
        // count shape also checkes for validity
        if (countShapes_from_list(potential_nbs, len))
        {
            revert_state_from_list(potential_nbs, len);
            return 0;
        }
    }
    dn_sq = shapes_1[0] - n_sq_old;
    dn_tr = shapes_1[1] - n_tr_old;
    dn_rh = shapes_1[2] - n_rh_old;

    int dn_e = (3 * dn_tr + 4 * dn_sq + 4 * dn_rh + new_boundary_length - boundary_length) / 2;
    int dn_f = dn_tr + dn_sq + dn_rh;
    int dec = dn_f - dn_e;

    if (dec != 0)
    {
        revert_state_from_list(potential_nbs, len);
        return 0;
    }

    // Metropolis algorithm
    energy_cost = eps_gamma * (new_boundary_length - boundary_length);
    energy_cost += eps_sq * dn_sq + eps_tr * dn_tr + eps_rh * dn_rh;

    acc = exp(-beta * energy_cost);

    if ((acc < dsfmt_genrand()) || reject)
    {
        revert_state_from_list(potential_nbs, len);
        return 0;
    }

    boundary_length = new_boundary_length;
    n_sq += dn_sq;
    n_tr += dn_tr;
    n_rh += dn_rh;

    return 1;
}

int check_full_house_and_fill(int idx)
/*
 * Checks if idx is the center of a full house and
 * if it is, fills fh (optional) and returns 1.
 * returns 0 otherwise.
 */
{
    vertex *v = &lattice[idx];

    // Count bonds
    int nbond = 0;
    int dir = -1;
    for (int i = 0; i < 12; i++)
    {
        if (v->bonds[i] != -1)
        {
            nbond++;
            dir = i;
        }
    }

    // Bond count has to be 3
    if (nbond != 3)
        return 0;

    assert(v->bonds[dir] != -1);
    // Now we have to check bond angles. 150-60-150 degrees
    // equivalently 5-2-5, 5-5-2, 2-5-5 steps
    int conf525 = v->bonds[(dir + 5) % 12] != -1 && v->bonds[(dir + 7) % 12] != -1;
    int conf552 = v->bonds[(dir + 5) % 12] != -1 && v->bonds[(dir + 10) % 12] != -1;
    int conf255 = v->bonds[(dir + 2) % 12] != -1 && v->bonds[(dir + 7) % 12] != -1;

    int dir_tip, dir_lb, dir_rb;
    if (conf525)
    {
        dir_tip = dir;
        dir_lb = (dir + 5) % 12;
        dir_rb = (dir + 7) % 12;
    }
    else if (conf552)
    {
        dir_rb = dir;
        dir_tip = (dir + 5) % 12;
        dir_lb = (dir + 10) % 12;
    }
    else if (conf255)
    {
        dir_lb = dir;
        dir_rb = (dir + 2) % 12;
        dir_tip = (dir + 7) % 12;
    }
    else
    {
        return 0;
    }

    vertex *tip = &lattice[v->bonds[dir_tip]];
    int invdir_tip = (dir_tip + 6) % 12;
    assert(tip->bonds[invdir_tip] == idx); // Sanity check

    // Check if lm exists
    int lm_idx = tip->bonds[(invdir_tip + 11) % 12];
    if (lm_idx == -1)
        return 0;

    // Check if rm exists
    int rm_idx = tip->bonds[(invdir_tip + 1) % 12];
    if (rm_idx == -1)
        return 0;

    // Check if lm is connected to lb
    vertex *lb = &lattice[v->bonds[dir_lb]];
    int invdir_lb = (dir_lb + 6) % 12;
    assert(lb->bonds[invdir_lb] == idx); // Sanity check

    // assert(lm_idx == lb->bonds[(invdir_lb + 1) % 12]);
    if (lm_idx != lb->bonds[(invdir_lb + 1) % 12])
        return 0;

    // Check if rm is connected to rb
    vertex *rb = &lattice[v->bonds[dir_rb]];
    int invdir_rb = (dir_rb + 6) % 12;
    assert(rb->bonds[invdir_rb] == idx); // Sanity check

    // assert(rm_idx == rb->bonds[(invdir_rb + 11) % 12]);
    if (rm_idx != rb->bonds[(invdir_rb + 11) % 12])
        return 0;

    // Check if lb is connected to rb
    if (v->bonds[dir_rb] != lb->bonds[(invdir_lb + 10) % 12])
        return 0;

    full_house *fh = &full_houses[n_full_houses];

    fh->o = idx;
    fh->tip = v->bonds[dir_tip];
    fh->lm = lm_idx;
    fh->rm = rm_idx;
    fh->lb = v->bonds[dir_lb];
    fh->rb = v->bonds[dir_rb];

    fh->dir_o_tip = dir_tip;

    n_full_houses++;

    assert(fh->dir_o_tip == dir_tip);
    assert((fh->dir_o_tip + 5) % 12 == dir_lb);
    assert((fh->dir_o_tip + 7) % 12 == dir_rb);
    assert(lattice[fh->o].bonds[fh->dir_o_tip] == fh->tip);
    assert(lattice[fh->o].bonds[(fh->dir_o_tip + 5) % 12] == fh->lb);
    assert(lattice[fh->o].bonds[(fh->dir_o_tip + 7) % 12] == fh->rb);

    return 1;
}

int check_empty_house_and_fill(int idx)
/*
 * Checks if idx is the lm of [an] empty house(s) and
 * if it is, fills eh (optional) and returns 1 or 2.
 * returns 0 otherwise.
 */
{
    vertex *v = &lattice[idx];

    // Count bonds
    int nbond = 0;
    int dir = -1;
    for (int i = 0; i < 12; i++)
    {
        if (v->bonds[i] != -1)
        {
            nbond++;
            dir = i;
        }
    }

    // Bond count has to be at least 3
    if (nbond < 3)
        return 0;
    int h_count = 0;
    for (int dir_rm = 0; dir_rm < 12; dir_rm++)
    {
        if (v->bonds[dir_rm] != -1 &&
            v->bonds[(dir_rm + 1) % 12] == -1 &&
            v->bonds[(dir_rm + 2) % 12] != -1 &&
            v->bonds[(dir_rm + 9) % 12] != -1 &&
            v->bonds[(dir_rm + 10) % 12] == -1 &&
            v->bonds[(dir_rm + 11) % 12] == -1)
        {
            int dir_tip = (dir_rm + 2) % 12;
            int dir_lb = (dir_rm + 9) % 12;
            int lb = v->bonds[dir_lb];
            int rb = lattice[lb].bonds[dir_rm];

            if (rb == -1 || rb != lattice[v->bonds[dir_rm]].bonds[dir_lb])
                continue;

            int rm = v->bonds[dir_rm];
            int tip = v->bonds[dir_tip];
            if (lattice[rm].bonds[(dir_rm + 4) % 12] != tip)
                continue;

            // Sanity Checks
            assert(lattice[v->bonds[dir_rm]].bonds[(dir_rm + 6) % 12] == idx);
            assert(lattice[v->bonds[dir_lb]].bonds[(dir_lb + 6) % 12] == idx);
            assert(lattice[v->bonds[dir_tip]].bonds[(dir_tip + 6) % 12] == idx);

            empty_houses[n_empty_houses].dir_lm_rm = dir_rm;
            empty_houses[n_empty_houses].tip = tip;
            empty_houses[n_empty_houses].lm = idx;
            empty_houses[n_empty_houses].lb = v->bonds[dir_lb];
            empty_houses[n_empty_houses].rm = rm;
            empty_houses[n_empty_houses].rb = rb;

            n_empty_houses++;
            h_count++;
        }
    }

    return h_count;
}

void update_full_houses(void)
{
    n_full_houses = 0;
    for (int i = 0; i < N; i++)
        check_full_house_and_fill(i);
}

void update_empty_houses(void)
{
    n_empty_houses = 0;
    int count = 0;
    for (int i = 0; i < N; i++)
        count += check_empty_house_and_fill(i);

    assert(count == n_empty_houses);
}

int not_swapable(const full_house *fh, const empty_house *eh)
{
    if (fh->lb == eh->lm && fh->rb == eh->rm)
    {
        return 1;
    }
    return 0;
}

int swapHouses(void)
{
    save_state();
    update_empty_houses();
    update_full_houses();

    if (n_empty_houses == 0 || n_full_houses == 0)
        return 0;

    int idx_fh = n_full_houses * dsfmt_genrand();
    int idx_eh = n_empty_houses * dsfmt_genrand();
    full_house *fh = &full_houses[idx_fh];
    empty_house *eh = &empty_houses[idx_eh];

    if (not_swapable(fh, eh) ||
        is_boundary(&lattice[eh->lb]) ||
        is_boundary(&lattice[eh->rb]) ||
        is_boundary(&lattice[eh->lm]) ||
        is_boundary(&lattice[eh->rm]) ||
        is_boundary(&lattice[eh->tip]) ||
        is_boundary(&lattice[fh->lb]) ||
        is_boundary(&lattice[fh->rb]) ||
        is_boundary(&lattice[fh->lm]) ||
        is_boundary(&lattice[fh->rm]) ||
        is_boundary(&lattice[fh->tip]) ||
        is_boundary(&lattice[fh->o]))
    {
        return 0;
    }

    // Move is accepted, now we need to perform it
    // Step 1: remove the central vertex from full house
    int o = fh->o;

    assert(lattice[o].bonds[fh->dir_o_tip] != -1);
    assert(lattice[o].bonds[(fh->dir_o_tip + 5) % 12] != -1);
    assert(lattice[o].bonds[(fh->dir_o_tip + 7) % 12] != -1);
    lattice[o].bonds[fh->dir_o_tip] = -1;
    lattice[o].bonds[(fh->dir_o_tip + 5) % 12] = -1;
    lattice[o].bonds[(fh->dir_o_tip + 7) % 12] = -1;

    assert(lattice[fh->tip].bonds[(fh->dir_o_tip + 6) % 12] == o);
    assert(lattice[fh->lb].bonds[(fh->dir_o_tip + 11) % 12] == o);
    assert(lattice[fh->rb].bonds[(fh->dir_o_tip + 1) % 12] == o);
    lattice[fh->tip].bonds[(fh->dir_o_tip + 6) % 12] = -1;
    lattice[fh->lb].bonds[(fh->dir_o_tip + 11) % 12] = -1;
    lattice[fh->rb].bonds[(fh->dir_o_tip + 1) % 12] = -1;

    // Step 1.5: Form the lm-rm bond
    assert(lattice[fh->lm].bonds[(fh->dir_o_tip + 9) % 12] == -1); // check
    assert(lattice[fh->rm].bonds[(fh->dir_o_tip + 3) % 12] == -1); // check
    lattice[fh->lm].bonds[(fh->dir_o_tip + 9) % 12] = fh->rm;
    lattice[fh->rm].bonds[(fh->dir_o_tip + 3) % 12] = fh->lm;

    // Step 2: remove square-triangle edge from the empty house
    assert(lattice[eh->lm].bonds[eh->dir_lm_rm] == eh->rm);
    assert(lattice[eh->rm].bonds[(eh->dir_lm_rm + 6) % 12] == eh->lm);
    lattice[eh->lm].bonds[eh->dir_lm_rm] = -1;
    lattice[eh->rm].bonds[(eh->dir_lm_rm + 6) % 12] = -1;

    // Step 3: add the central vertex
    int dir_o_tip = (eh->dir_lm_rm + 3) % 12;

    lattice[o].bonds[dir_o_tip] = eh->tip;
    lattice[o].bonds[(dir_o_tip + 5) % 12] = eh->lb;
    lattice[o].bonds[(dir_o_tip + 7) % 12] = eh->rb;

    lattice[eh->tip].bonds[(dir_o_tip + 6) % 12] = o;
    lattice[eh->lb].bonds[(dir_o_tip + 11) % 12] = o;
    lattice[eh->rb].bonds[(dir_o_tip + 1) % 12] = o;

    // Step 4: update position
    lattice[o].n0 = lattice[eh->lb].n0 + lattice[eh->tip].n0 - lattice[eh->lm].n0;
    lattice[o].n1 = lattice[eh->lb].n1 + lattice[eh->tip].n1 - lattice[eh->lm].n1;
    lattice[o].n2 = lattice[eh->lb].n2 + lattice[eh->tip].n2 - lattice[eh->lm].n2;
    lattice[o].n3 = lattice[eh->lb].n3 + lattice[eh->tip].n3 - lattice[eh->lm].n3;

    calculateParCoords(&lattice[o]);

    return 1;
}

void handle_sigabrt(int sig)
{
    writeLattice(-1);
    writeoldLattice(-2);
}

int writeBoundary(int step)
{
    char datafile[128];

    sprintf(datafile, "%s/boundary%09d.dat", datadir, step);
    FILE *data = fopen(datafile, "w");

    int length, top, start, neighbour, start_dir = -1, prev_dir;
    vertex current;

    // Find a vertex at the boundary (at the top)
    top = 0;
    for (int i = 0; i < N; i++)
    {
        if (lattice[i].par.y > lattice[top].par.y)
            top = i;
    }

    // Loop around until we find it again
    start = 6;
    prev_dir = start;
    length = 0;
    current = lattice[top];
    calculatePerpCoords(&lattice[top]);
    myfprintf(data, "%d %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", top, lattice[top].par.x, lattice[top].par.y, lattice[top].perp.x, lattice[top].perp.y, lattice[top].n0, lattice[top].n1, lattice[top].n2, lattice[top].n3, lattice[top].bonds[0], lattice[top].bonds[1], lattice[top].bonds[2], lattice[top].bonds[3], lattice[top].bonds[4], lattice[top].bonds[5], lattice[top].bonds[6], lattice[top].bonds[7], lattice[top].bonds[8], lattice[top].bonds[9], lattice[top].bonds[10], lattice[top].bonds[11]);
    for (int step = 0; step < 50000; step++)
    {
        for (int direction = start; direction < start + 12; direction++)
        {
            neighbour = current.bonds[direction % 12];
            if (neighbour == -1)
                continue;
            if (length == 0)
                start_dir = direction;

            current = lattice[neighbour];
            length++;
            calculatePerpCoords(&lattice[neighbour]);
            myfprintf(data, "%d %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", neighbour, lattice[neighbour].par.x, lattice[neighbour].par.y, lattice[neighbour].perp.x, lattice[neighbour].perp.y, lattice[neighbour].n0, lattice[neighbour].n1, lattice[neighbour].n2, lattice[neighbour].n3, lattice[neighbour].bonds[0], lattice[neighbour].bonds[1], lattice[neighbour].bonds[2], lattice[neighbour].bonds[3], lattice[neighbour].bonds[4], lattice[neighbour].bonds[5], lattice[neighbour].bonds[6], lattice[neighbour].bonds[7], lattice[neighbour].bonds[8], lattice[neighbour].bonds[9], lattice[neighbour].bonds[10], lattice[neighbour].bonds[11]);
            if (current.index == top)
            {
                assert(start_dir != -1);
                if (countNeighbours(&lattice[top]) == 1)
                {
                    fclose(data);
                    return length;
                }

                else if (direction != (start_dir + 6) % 12 && direction != (prev_dir + 6) % 12)
                {
                    fclose(data);
                    return length;
                }
            }
            prev_dir = direction;
            start = (direction + 7) % 12;
            break;
        }
    }
    fclose(data);
    return -1;
}

int writeoldBoundary(int step)
{
    char datafile[128];

    sprintf(datafile, "%s/boundary%09d.dat", datadir, step);
    FILE *data = fopen(datafile, "w");

    int length, top, start, neighbour, start_dir = -1, prev_dir, done;
    vertex current;

    // Find a vertex at the boundary (at the top)
    top = 0;
    for (int i = 0; i < N; i++)
    {
        if (oldlattice[i].par.y > oldlattice[top].par.y)
            top = i;
    }

    // Loop around until we find it again
    done = 0;
    start = 6;
    prev_dir = start;
    length = 0;
    current = oldlattice[top];
    calculatePerpCoords(&oldlattice[top]);
    myfprintf(data, "%d %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", top, oldlattice[top].par.x, oldlattice[top].par.y, oldlattice[top].perp.x, oldlattice[top].perp.y, oldlattice[top].n0, oldlattice[top].n1, oldlattice[top].n2, oldlattice[top].n3, oldlattice[top].bonds[0], oldlattice[top].bonds[1], oldlattice[top].bonds[2], oldlattice[top].bonds[3], oldlattice[top].bonds[4], oldlattice[top].bonds[5], oldlattice[top].bonds[6], oldlattice[top].bonds[7], oldlattice[top].bonds[8], oldlattice[top].bonds[9], oldlattice[top].bonds[10], oldlattice[top].bonds[11]);
    for (int step = 0; step < 50000; step++)
    {
        for (int direction = start; direction < start + 12; direction++)
        {
            neighbour = current.bonds[direction % 12];
            if (neighbour == -1)
                continue;
            if (length == 0)
                start_dir = direction;

            current = oldlattice[neighbour];
            length++;
            calculatePerpCoords(&oldlattice[neighbour]);
            myfprintf(data, "%d %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", neighbour, oldlattice[neighbour].par.x, oldlattice[neighbour].par.y, oldlattice[neighbour].perp.x, oldlattice[neighbour].perp.y, oldlattice[neighbour].n0, oldlattice[neighbour].n1, oldlattice[neighbour].n2, oldlattice[neighbour].n3, oldlattice[neighbour].bonds[0], oldlattice[neighbour].bonds[1], oldlattice[neighbour].bonds[2], oldlattice[neighbour].bonds[3], oldlattice[neighbour].bonds[4], oldlattice[neighbour].bonds[5], oldlattice[neighbour].bonds[6], oldlattice[neighbour].bonds[7], oldlattice[neighbour].bonds[8], oldlattice[neighbour].bonds[9], oldlattice[neighbour].bonds[10], oldlattice[neighbour].bonds[11]);
            if (current.index == top)
            {
                assert(start_dir != -1);
                if (countNeighbours(&oldlattice[top]) == 1)
                {
                    fclose(data);
                    return length;
                }

                else if (direction != (start_dir + 6) % 12 && direction != (prev_dir + 6) % 12)
                {
                    fclose(data);
                    return length;
                }
            }
            prev_dir = direction;
            start = (direction + 7) % 12;
            break;
        }
    }
    fclose(data);
    return -1;
}

void go_to_state(double eps_start, double eps_finish, double gamma_start, double gamma_finish, int n_steps)
{
    double deps = (eps_finish - eps_start) / n_steps;
    double dgamma = (gamma_finish - gamma_start) / n_steps;
    int sanity_check_interval = 100;
    for (int step = 1; step < n_steps; step++)
    {
        eps_rh = eps_start + step * deps;
        eps_gamma = gamma_start + step * dgamma;
        for (int nn = 0; nn < N; nn++)
        {
            int acc;
            if (dsfmt_genrand() > 1. / (N + 1))
            {
                acc = moveVertex();
            }
            else
            {
                acc = swapHouses();
            }
            int n_v = N;
            int n_e = (3 * n_tr + 4 * n_sq + 4 * n_rh + boundary_length) / 2;
            int n_f = n_tr + n_sq + n_rh + 1;
            int ec = n_v - n_e + n_f;
            if (ec != 2)
            {
                revert_state();
                acc = 0;
            }
            if (acc)
            {
                save_state();
            }
            // Take measurements every number of steps
            if (step % sanity_check_interval == 0)
            {
                countAllShapes(); // To check if n_sq,n_tr,n_rh are still correct
                if (n_sq != shapes[0] || n_tr != shapes[1] || n_rh != shapes[2])
                {

                    if (errcount == MAXERR)
                    {
                        myfprintf(stderr, "\nShape counting mismatch.\nMaximum number of reversions reached.\n\n");
                        assert(n_sq == shapes[0]);
                        assert(n_tr == shapes[1]);
                        assert(n_rh == shapes[2]);
                    }
                    myfprintf(stderr, "\nShape counting mismatch.\nReverting back to last stable sweep.\n\n");

                    revert_backup_state();
                    step -= sanity_check_interval;
                    errcount++;
                    continue;
                }

                save_backup_state();
            }
        }
    }
    eps_rh = eps_finish;
    eps_gamma = gamma_finish;
}

int main(int argc, char *argv[])
{
    parse_cli(argc, argv);
    signal(SIGABRT, handle_sigabrt);

    int accepted = 0;
    int swap_accepted = 0;
    int proposed = 0;

    // Seed RNG
    uint32_t seed = 3141519;
    FILE *fpp = fopen("/dev/urandom", "r");
    unsigned long tmp = fread(&seed, 1, sizeof(uint32_t), fpp);
    if (tmp != sizeof(uint32_t))
        printf("error with seed\n");
    fclose(fpp);
    myprintf("Seed: %u\n", (uint32_t)seed);
    fflush(stdout);
    dsfmt_seed(seed);

    // Generate D x D lattice
    genLattice();
    //
    //    // Stampfli inflate the lattice st_I times
    for (int inflation_step = 0; inflation_step < st_I; inflation_step++)
          inflateLattice();

    //    genChevronLattice();

    // Find the neighbours for all vertices
    for (int n = 0; n < N; n++)
        calculateParCoords(&lattice[n]);

    for (int i = 0; i < N; i++)
        initNeighbours(&lattice[i]);

    for (int n = 0; n < N; n++)
        removeBondCross(&lattice[n]);

    for (int n = 0; n < N; n++)
        calculatePerpCoords(&lattice[n]);

    empty_houses = calloc(2 * N, sizeof(empty_house));
    full_houses = calloc(2 * N, sizeof(full_house));

    init_datadir();
    makeHist();
    writeLattice(0);

    init_top_vertex();
    boundary_length = countBoundaryLength();
    countAllShapes();
    n_sq = shapes[0];
    n_tr = shapes[1];
    n_rh = shapes[2];

    write_phase_state();
    init_measurement_file();

    myprintf("Initial boundary length = %d\n", boundary_length);
    myprintf("Number of vertices = %d\n", N);
    update_full_houses();
    update_empty_houses();
    myprintf("# of fHouse \t\t= %d\n", n_full_houses);
    myprintf("# of eHouse \t\t= %d\n", n_empty_houses);

    save_state();
    save_backup_state();
    go_to_state(STABLE_EPS_RH, eps_rh, STABLE_GAMMA, eps_gamma, ANNEALING_SWEEPS);
    write_measurement(0);
    writeLattice(0);
    int ec = 2, oec = 2;
    // Move vertices
    clock_t t = clock();
    for (sweep = 1; sweep < MCsweeps + 1; sweep++)
    {
        for (int nn = 0; nn < N; nn++)
        {
            proposed++;
            oec = ec;
            int acc;
            if (dsfmt_genrand() > 1. / (N + 1))
            {
                acc = moveVertex();
            }
            else
            {
                acc = swapHouses();
                swap_accepted += acc;
            }

            int n_v = N;
            int n_e = (3 * n_tr + 4 * n_sq + 4 * n_rh + boundary_length) / 2;
            int n_f = n_tr + n_sq + n_rh + 1;
            ec = n_v - n_e + n_f;
            if (ec != 2)
            {
                revert_state();
                acc = 0;
            }
            accepted += acc;
            if (acc)
            {
                save_state();
            }
        }

        // Take measurements every number of steps
        if (sweep % output_interval == 0)
        {
            countAllShapes(); // To check if n_sq,n_tr,n_rh are still correct
            if (n_sq != shapes[0] || n_tr != shapes[1] || n_rh != shapes[2])
            {

                if (errcount == MAXERR)
                {
                    myfprintf(stderr, "\nShape counting mismatch.\nMaximum number of reversions reached.\n\n");
                    assert(n_sq == shapes[0]);
                    assert(n_tr == shapes[1]);
                    assert(n_rh == shapes[2]);
                }
                myfprintf(stderr, "\nShape counting mismatch.\nReverting back to last output sweep.\n\n");

                revert_backup_state();
                accepted = 0;
                proposed = 0;
                swap_accepted = 0;
                countall = 0;
                sweep -= output_interval;
                errcount++;
                continue;
            }

            makeHist();
            writeLattice(sweep);
            write_measurement(sweep);
            update_full_houses();
            update_empty_houses();

            myprintf("Gamma = %.1lf\n", eps_gamma);
            myprintf("Move %d:\n", sweep);

            myprintf("\tAcceptance ratio \t= %lf\n", (double)accepted / proposed);
            myprintf("\tSwap acceptance \t= %lf\n", (double)swap_accepted / proposed);

            myprintf("\tBoundary length \t= %d\n", boundary_length);
            myprintf("\t# of squares \t\t= %d (%d)\n", n_sq, shapes[0]);
            myprintf("\t# of triangles \t\t= %d (%d)\n", n_tr, shapes[1]);
            myprintf("\t# of rhombi \t\t= %d (%d)\n", n_rh, shapes[2]);
            myprintf("\t# of fHouse \t\t= %d\n", n_full_houses);
            myprintf("\t# of eHouse \t\t= %d\n", n_empty_houses);
            myprintf("\tV - E + F \t\t\t= %d\n", ec);

            assert(n_sq == shapes[0]);
            assert(n_tr == shapes[1]);
            assert(n_rh == shapes[2]);

            myprintf("\tCumulative runtime \t= %.2lf s\n", (clock() - t) / (double)(CLOCKS_PER_SEC));
            myprintf("\n");

            accepted = 0;
            proposed = 0;
            swap_accepted = 0;
            countall = 0;
            save_backup_state();
        }
    }
    myprintf("Total runtime: %.2lf s\n", (clock() - t) / (double)(CLOCKS_PER_SEC));

    FILE *ff = fopen("finished.txt", "a");
    myfprintf(ff, "%s\n", datadir);
    fclose(ff);
    return 0;
}
