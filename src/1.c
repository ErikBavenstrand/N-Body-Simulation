#ifndef _REENTRANT
#define _REENTRANT
#endif

/* includes */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/times.h>
#include <sys/stat.h>

/* defines */
#define MAX_BODIES 1000000
#define MAX_STEPS 500000
#define MAX_SIZE 1000000
#define MIN_DISTANCE_CAL 0.0001
#define G 6.67E-8
#define DT 1

/* a single point */
typedef struct {
    double x;
    double y;
} Point;

/* a single body */
typedef struct {
    Point position;
    Point velocity;
    Point force;
    double mass;
} Body;

/* global variables */
Body *bodies;

int num_bodies;
int num_steps;
int total_size;
int seed;

static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;

void calculateForces();
void moveBodies();
void initializeBodies();
void initializeSolarsystem();
void start_clock();
void end_clock();

int main(int argc, char *argv[]) {
    int i, j;/* for iterating */

    /* create folder for output */
    #ifdef VISUALIZE
    struct stat st = {0};
    if (stat("./outfiles", &st) == -1) {
        mkdir("./outfiles", 0700);
    }
    FILE *fp = fopen("./outfiles/1.out","w");
    #endif

    /* check if arg count is correct */
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_bodies num_steps total_size seed\n", argv[0]);
        fflush(stderr);
        exit(-1);
    }

    /* assign values from input */
    num_bodies = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    total_size = atoi(argv[3]);
    seed = atoi(argv[4]);

    /* sanity check all important input values */
    if (num_bodies > MAX_BODIES || num_bodies <= 0) num_bodies = MAX_BODIES;
    if (num_steps > MAX_STEPS || num_steps <= 0) num_steps = MAX_STEPS;
    if (total_size > MAX_SIZE || total_size <= 0) total_size = MAX_SIZE;

    /* initialize bodies */
    #ifdef SOLAR_SYSTEM
    initializeSolarsystem();
    #else
    initializeBodies();
    #endif

    /* print initial line to file */
    #ifdef VISUALIZE
    fprintf(fp, "%d %d\n", num_bodies, 1);
    #endif
    
    start_clock();/* start time */
    for (i = 0; i < num_steps; i++) {
        calculateForces();/* calculate forces on all the bodies */
        moveBodies();/* move the bodies according to the forces */

        /* print locations to file */
        #ifdef VISUALIZE
        for (j = 0; j < num_bodies; j++) {
            fprintf(fp, "%lf %lf\n", bodies[j].position.x / 1000, bodies[j].position.y / 1000);
            fflush(fp);
        }
        #endif
    }
    end_clock();/* end time */

    /* close file */
    #ifdef VISUALIZE
    fclose(fp);
    #endif

    free(bodies);/* free the heap */
    return 0;
}

/* calculate forces on all bodies */
void calculateForces() {
    int i, j;/* for iteration */
    double distance;
    double magnitude;
    Point direction;
    /* loop all bodies for each body */
    for (i = 0; i < num_bodies - 1; i++) {
        for (j = i + 1; j < num_bodies; j++) {
            distance = sqrt(pow(bodies[i].position.x - bodies[j].position.x, 2) + 
                            pow(bodies[i].position.y - bodies[j].position.y, 2));
            
            /* if distance is too small */
            if (distance < MIN_DISTANCE_CAL) {
                distance = MIN_DISTANCE_CAL;
            }
            
            magnitude = (G * bodies[i].mass * bodies[j].mass) / pow(distance, 2);
            direction.x = bodies[j].position.x - bodies[i].position.x;
            direction.y = bodies[j].position.y - bodies[i].position.y;

            bodies[i].force.x += magnitude * direction.x / distance;
            bodies[i].force.y += magnitude * direction.y / distance;

            bodies[j].force.x -= magnitude * direction.x / distance;
            bodies[j].force.y -= magnitude * direction.y / distance;
        }
    }
}

/* move all bodies according to their forces */
void moveBodies() {
    int i;/* for iteration */
    Point delta_v;
    Point delta_p;
    for (i = 0; i < num_bodies; i++) {
        delta_v.x = (bodies[i].force.x / bodies[i].mass) * DT;
        delta_v.y = (bodies[i].force.y / bodies[i].mass) * DT;

        delta_p.x = (bodies[i].velocity.x + delta_v.x / 2) * DT;
        delta_p.y = (bodies[i].velocity.y + delta_v.y / 2) * DT;

        bodies[i].velocity.x += delta_v.x;
        bodies[i].velocity.y += delta_v.y;

        bodies[i].position.x += delta_p.x;
        bodies[i].position.y += delta_p.y;

        /* keep bodies inside the space */
        #ifdef BOUNDS
        if (bodies[i].position.x > total_size) {
            bodies[i].velocity.x *= -1;
            bodies[i].position.x = total_size;
        }
        if (bodies[i].position.x < 0) {
            bodies[i].velocity.x *= -1;
            bodies[i].position.x = 0;
        }
        if (bodies[i].position.y > total_size) {
            bodies[i].velocity.y *= -1;
            bodies[i].position.y = total_size;
        } 
        if (bodies[i].position.y < 0) {
            bodies[i].velocity.y *= -1;
            bodies[i].position.y = 0;
        }
        #endif 

        /* important to resent the forces */
        bodies[i].force.x = 0;
        bodies[i].force.y = 0;
    }
}

/* initialize the bodies */
void initializeBodies() {
    int i;/* for iteration */
    bodies = (Body *)malloc(sizeof(Body) * num_bodies);

    /* seed the random number generator */
    #ifdef RAND
    time_t t;
    srand((unsigned) time(&t));
    #else
    srand(seed);
    #endif

    /* gererate */
    for (i = 0; i < num_bodies; i++) {
        bodies[i].position.x = rand()%total_size;
        bodies[i].position.y = rand()%total_size;
        bodies[i].mass = rand()%10000000;
        bodies[i].force.x = 0;
        bodies[i].force.y = 0;
        bodies[i].velocity.x = 0;
        bodies[i].velocity.y = 0;
    }
}

/* initialize a solar system */
void initializeSolarsystem() {
    num_bodies = 2;
    bodies = (Body *)malloc(sizeof(Body) * num_bodies);
    bodies[0].position.x = 500;
    bodies[0].position.y = 500;
    bodies[0].mass = 5000000;
    bodies[0].force.x = 0;
    bodies[0].force.y = 0;
    bodies[0].velocity.x = 0;
    bodies[0].velocity.y = 0;
    bodies[1].position.x = 450;
    bodies[1].position.y = 500;
    bodies[1].mass = 50000;
    bodies[1].force.x = 0;
    bodies[1].force.y = 0;
    bodies[1].velocity.x = 0;
    bodies[1].velocity.y = 0.1;
}

/* start the clock */
void start_clock() {
    st_time = times(&st_cpu);
}

/* end the clock */
void end_clock() {
    en_time = times(&en_cpu);
    printf("%g\n", (double)(en_time - st_time)/sysconf(_SC_CLK_TCK));
}