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
#include <pthread.h>

/* defines */
#define MAX_BODIES 1000000
#define MAX_STEPS 500000
#define MAX_WORKERS 16
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
Point **forces;

#ifdef VISUALIZE
FILE *fp;
#endif

int num_bodies;
int num_steps;
int num_workers;
int num_arrived;
int total_size;
int seed;

pthread_mutex_t barrier_lock;
pthread_cond_t go;

static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;

void *calculationWrapper(void *arg);
void calculateForces(long id);
void moveBodies(long id);
void initializeBodies();
void initializeSolarsystem();
void barrier();
void start_clock();
void end_clock();

int main(int argc, char *argv[]) {
    long i;/* for iteration */

    /* create folder for output */
    #ifdef VISUALIZE
    struct stat st = {0};
    if (stat("./outfiles", &st) == -1) {
        mkdir("./outfiles", 0700);
    }
    fp = fopen("./outfiles/2.out","w");
    #endif

    /* check if arg count is correct */
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_bodies num_steps num_workers total_size seed\n", argv[0]);
        fflush(stderr);
        exit(-1);
    }

    /* set values from input */
    num_bodies = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    num_workers = atoi(argv[3]);
    total_size = atoi(argv[4]);
    seed = atoi(argv[5]);

    /* sanity check all input values */
    if (num_bodies > MAX_BODIES || num_bodies <= 0) num_bodies = MAX_BODIES;
    if (num_steps > MAX_STEPS || num_steps <= 0) num_steps = MAX_STEPS;
    if (num_workers > MAX_WORKERS || num_workers <= 0) num_workers = MAX_WORKERS;
    if (total_size > MAX_SIZE || total_size <= 0) total_size = MAX_SIZE;

    /* initialize workers etc. */
    pthread_t workers[num_workers];
    pthread_mutex_init(&barrier_lock, NULL);
    pthread_cond_init(&go, NULL);

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

    /* create workers */
    for (i = 0; i < num_workers; i++) {
        pthread_create(&workers[i], NULL, calculationWrapper, (void *)i);
    }

    /* wait for completiton */
    for (i = 0; i < num_workers; i++) {
        pthread_join(workers[i], NULL);
    }
    end_clock();/* end time */

    #ifdef VISUALIZE
    fclose(fp);
    #endif

    /* free allocated memory from heap */
    for (i = 0; i < num_workers; i++) {
        free(forces[i]);
    }
    free(forces);
    free(bodies);
    return 0;
}

/* wrapper for the calculations */
void *calculationWrapper(void *arg) {
    long id = (long)arg;/* id of process */
    int i, j;/* for iteration */
    for (i = 0; i < num_steps; i++) {
        calculateForces(id);/* calculate forces on bodies */
        barrier();/* wait for all threads to complete */
        moveBodies(id);/* move the bodies */
        barrier();/* wait for all threads to complete */

        /* print positions to file */
        #ifdef VISUALIZE
        if (id == 0) {
            for (j = 0; j < num_bodies; j++) {
                fprintf(fp, "%lf %lf\n", bodies[j].position.x / 1000, bodies[j].position.y / 1000);
                fflush(fp);
            }
        }
        #endif
    }
    return 0;
}

/* calculate forces on bodies */
void calculateForces(long id) {
    long i, j;/* for iteration */
    double distance;
    double magnitude;
    Point direction;
    
    /* each thread computes every num_workers'th body and stores results in own section of forces array */
    for (i = id; i < num_bodies; i += num_workers) {
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

            forces[id][i].x += magnitude * direction.x / distance;
            forces[id][i].y += magnitude * direction.y / distance;

            forces[id][j].x -= magnitude * direction.x / distance;
            forces[id][j].y -= magnitude * direction.y / distance;
        }
    }
}

/* move bodies according to forces */
void moveBodies(long id) {
    long i, j;/* for iteration */
    Point delta_v;
    Point delta_p;
    for (i = id; i < num_bodies; i += num_workers) {
        
        /* summarize forces from force array and reset */
        for (j = 0; j < num_workers; j++) {
            bodies[i].force.x += forces[j][i].x;
            bodies[i].force.y += forces[j][i].y;
            forces[j][i].x = 0;
            forces[j][i].y = 0;
        }
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

/* initialize bodies */
void initializeBodies() {
    int i, j;/* for iteration */

    /* allocate bodies and forces array */
    bodies = (Body *)malloc(sizeof(Body) * num_bodies);
    forces = malloc(sizeof(*forces) * num_workers);
    if (forces) {
        for (i = 0; i < num_workers; i++) {
            forces[i] = malloc(sizeof(*forces[i]) * num_bodies);
            for (j = 0; j < num_bodies; j++) {
                forces[i][j].x = 0;
                forces[i][j].y = 0;
            }
        }
    }

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
        bodies[i].mass = rand()%100000;
        bodies[i].force.x = 0;
        bodies[i].force.y = 0;
        bodies[i].velocity.x = 0;
        bodies[i].velocity.y = 0;
    }
}


/* initialize a solar system */
void initializeSolarsystem() {
    int i, j;
    num_bodies = 2;
    bodies = (Body *)malloc(sizeof(Body) * num_bodies);
    forces = malloc(sizeof(*forces) * num_workers);
    if (forces) {
        for (i = 0; i < num_workers; i++) {
            forces[i] = malloc(sizeof(*forces[i]) * num_bodies);
            for (j = 0; j < num_bodies; j++) {
                forces[i][j].x = 0;
                forces[i][j].y = 0;
            }
        }
    }
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

/* barrier for waiting for threads */
void barrier() {
  pthread_mutex_lock(&barrier_lock);
  num_arrived++;
  if (num_arrived == num_workers) {
    num_arrived = 0;
    pthread_cond_broadcast(&go);
  } else
    pthread_cond_wait(&go, &barrier_lock);
  pthread_mutex_unlock(&barrier_lock);
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