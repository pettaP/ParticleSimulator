#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;
double intervall;
const int sizesteps = 100;

square_t squares[sizesteps][sizesteps];
square_t *previousSquares[1000];
int squareCounter = 0;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
    intervall = size / sizesteps;
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );

    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;

    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;

    for( int i = 0; i < n; i++ )
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];

        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

void initSquares(){
    for(int i = 0; i < sizesteps; i++){
        for(int j = 0; j < sizesteps; j++){
            squares[i][j].trueNeighbours = false;
            squares[i][j].occupied = false;
            squares[i][j].particles = nullptr;
        }
    }
}

void clearEnvironment(){
    for(int i = 0; i < squareCounter; i++){
        previousSquares[i]->occupied = false;
        previousSquares[i]->trueNeighbours = false;
        freeNodes(previousSquares[i]->particles);
    }

    squareCounter = 0;
}

void freeNodes(particle_node_t* destroyNode){
    if(destroyNode->next == nullptr) free(destroyNode);
    else{
        freeNodes(destroyNode->next);
        free(destroyNode);
    }
}

void putInSquare(particle_t *particle){
    auto x = static_cast<int>(std::floor(particle->x / intervall));
    auto y = static_cast<int>(std::floor(particle->y / intervall));
    particle->inMiddle = true;
    //mutex
    particle_node_t * rest;
    rest = squares[x][y].particles->next;
    particle_node_t * ny;
    ny = (particle_node_t*) malloc(sizeof(particle_node_t));
    ny->p = particle;
    ny->next = rest;

    if(!squares[x][y].occupied) {
        squares[x][y].occupied = true;
        previousSquares[squareCounter++] = &squares[x][y];
    }

    squares[x][y].particles = ny;
    //unmutex

    if((x * intervall <= (particle->x)) && ((particle->x) <= (x * intervall) + cutoff*cutoff)){
        squares[x][y].trueNeighbours = true;
        particle->inMiddle = false;
    } else if ((((x + 1) * intervall - cutoff*cutoff) <= (particle->x)) && ((particle->x) <= ((x + 1) * intervall))){
        squares[x][y].trueNeighbours = true;
        particle->inMiddle = false;
    }
    if((y * intervall <= (particle->y)) && ((particle->y) <= (y * intervall) + cutoff*cutoff)){
        squares[x][y].trueNeighbours = true;
        particle->inMiddle = false;
    } else if ((((y + 1) * intervall - cutoff*cutoff) <= (particle->y)) && ((particle->y) <= ((y + 1) * intervall))){
        squares[x][y].trueNeighbours = true;
        particle->inMiddle = false;
    }

}

void applyForces(particle_t *particle){
    auto x = static_cast<int>(std::floor(particle->x / intervall));
    auto y = static_cast<int>(std::floor(particle->y / intervall));
    particle->ax = particle-> ay = 0;
    particle_node_t *temp;

    if(particle->inMiddle) {
        temp = squares[x][y].particles;
        while (temp != nullptr) {
            apply_force(*particle, *temp->p);
            temp = temp->next;
        }
    }else{
        int tempX;
        int tempY;
        int maxX,maxY;
        if(x > 0) tempX = x - 1; else tempX = x;
        if(y > 0) tempY = y - 1; else tempY = y;
        if(x < sizesteps - 1) maxX = x + 1; else maxX = x;
        if(y < sizesteps - 1) maxY = y + 1; else maxY = y;

        for (int i = tempX; i < maxX; i++) {
            for (int j = tempY; j < maxY; j++) {
                if (squares[i][j].trueNeighbours) {
                    temp = squares[i][j].particles;
                    while (temp != nullptr) {
                        apply_force(*particle, *temp->p);
                        temp = temp->next;
                    }
                }
            }
        }
    }
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor )
{
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff ) //0.0001
        return;
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
