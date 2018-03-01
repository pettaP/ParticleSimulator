#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct
{
    double x;
    double y;
    double vx;
    double vy;
    double ax;
    double ay;
    bool inMiddle;
} particle_t;

//
//particle node structure
//
struct particle_node_t
{
    particle_t * p;
    particle_node_t * next;
};

//
//square data structure
//
typedef struct
{
    bool occupied;
    bool trueNeighbours;
    particle_node_t *particles;
} square_t;

//
//
//

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );

void initSquares();
void clearEnvironment();
void freeNodes(particle_node_t* destroyNode);
void putInSquare(particle_t *particle);
void applyForces(particle_t *particle);

void apply_force( particle_t &particle, particle_t &neighbor );
void move( particle_t &p );

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
