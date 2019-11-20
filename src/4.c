#ifndef _REENTRANT
#define _REENTRANT
#endif

/* includes */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <pthread.h>

/* defines */
#define MAX_BODIES 1000000
#define MAX_STEPS 1000000
#define MAX_WORKERS 16
#define MAX_FAR 1000000
#define MAX_SIZE 1000000
#define MIN_DISTANCE_CAL 0.0001
#define MIN_DISTANCE_APP 1
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

/* a node in quad tree */
typedef struct Node{
    struct Node *ne;
    struct Node *nw;
    struct Node *sw;
    struct Node *se;
    Body body;
    Point center_of_mass;
    double mass;
    double size;
    bool isLeaf;
    bool hasBody;
    bool bodyIsHere;
    Point position;
    int id;
} Node;

/* global variables */
Body *bodies;
Node *root;
Node *childNodes[MAX_WORKERS];

int num_bodies;
int num_steps;
int far;
int total_size;
int num_workers;
int num_arrived;
int seed;

#ifdef VISUALIZE
FILE *fp;
#endif

pthread_mutex_t barrier_lock;
pthread_cond_t go;

static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;

void *calculationWrapper(void *arg);
void calculateForces(Body *body, Node *node);
void moveBodies(Body *body);
void calculateCenterOfMasses(Node *node);
void constructQuadTree(Body *body, Node *node);
Node *getBodyQuadrantInNode(Point position, Node *node);
bool isCorrectSubtree(long id_worker, int id_node);
bool checkFarAway(Point position);
void summarizeWeights();
void initializeBodies();
void initializeSolarsystem();
void initializeTree();
void initializeChildNodes(Node *node);
void freeNodes(Node* node);
void hasBodyRemove(Node* node);
void freeWasteNodes(Node* node);
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
    fp = fopen("./outfiles/4.out","w");
    #endif

    /* check if arg count is correct */
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_bodies num_steps num_workers far size seed\n", argv[0]);
        fflush(stderr);
        exit(-1);
    }

    /* set values from input */
    num_bodies = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    num_workers = atoi(argv[3]);
    far = atoi(argv[4]);
    total_size = atoi(argv[5]);
    seed = atoi(argv[6]);

    /* sanity check all input values */
    if (num_bodies > MAX_BODIES || num_bodies <= 0) num_bodies = MAX_BODIES;
    if (num_steps > MAX_STEPS || num_steps <= 0) num_steps = MAX_STEPS;
    if (num_workers > MAX_WORKERS || num_workers <= 0) num_workers = MAX_WORKERS;
    if (far > MAX_FAR || far < 0) far = MAX_FAR;
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

    start_clock();/* start timer */
    /* create workers */
    for (i = 0; i < num_workers; i++) {
        pthread_create(&workers[i], NULL, calculationWrapper, (void *)i);
    }
    pthread_exit(NULL);
    
}

void *calculationWrapper(void *arg) {
    long id = (long)arg;/* worker id */
    long i, j;/* for iteration */
    if (id == 0) {
        initializeTree();/* initialize the tree */
    }
    for (i = 0; i < num_steps; i++) {
        barrier();

        /* parallel construction of the quad tree */
        for (j = 0; j < num_bodies; j++) {
            Node *n = getBodyQuadrantInNode(bodies[j].position, getBodyQuadrantInNode(bodies[j].position, root)); /* find the correct second level quadrant */
            
            /* check if woker is assigned to this quadrant and is inside given area */
            if (isCorrectSubtree(id, n->id) && checkFarAway(bodies[j].position)) {
                constructQuadTree(&bodies[j], n);
            }
        }

        barrier();

        for (j = 0; j < MAX_WORKERS; j++) {
            if (isCorrectSubtree(id, childNodes[j]->id)) {
                freeWasteNodes(childNodes[j]);
            }
        }

        barrier();

        if (id == 0) {
            summarizeWeights();
            calculateCenterOfMasses(root);
        }

        barrier();

        /* calculate forces for bodies in parallel */
        for (j = id; j < num_bodies; j += num_workers) {
            calculateForces(&bodies[j], root);
        }

        barrier();

        /* move bodies according to forces in parallel */
        for (j = id; j < num_bodies; j += num_workers) {
            moveBodies(&bodies[j]);
        }

        barrier();

        for (j = 0; j < MAX_WORKERS; j++) {
            if (isCorrectSubtree(id, childNodes[j]->id)) {
                hasBodyRemove(childNodes[j]);
            }
        }

        barrier();
        if (id == 0) {
            hasBodyRemove(root);
            #ifdef VISUALIZE
            for (j = 0; j < num_bodies; j++) {
                fprintf(fp, "%lf %lf\n", bodies[j].position.x / 1000, bodies[j].position.y / 1000);
                fflush(fp);
            }
            #endif
        }
    }
    barrier();
    if (id == 0) {
        end_clock();/* end timer */

        /* close file */
        #ifdef VISUALIZE
        fclose(fp);
        #endif
        freeNodes(root);
        free(bodies);
        return 0;
    }
    pthread_exit(NULL);
}

/* calculate forces on a body */
void calculateForces(Body *body, Node *node) {
    double distance;
    double magnitude;
    Point direction;

    /* if node is a leaf and contains a body */
    if (node->isLeaf && node->hasBody) {
        distance = sqrt(pow(body->position.x - node->body.position.x, 2) +
                        pow(body->position.y - node->body.position.y, 2));
        
        /* if distance is too small */
        if (distance < MIN_DISTANCE_CAL) {
            distance = MIN_DISTANCE_CAL;
        }

        magnitude = (G * body->mass * node->body.mass) / pow(distance, 2);
        direction.x = node->body.position.x - body->position.x;
        direction.y = node->body.position.y - body->position.y;

        body->force.x += magnitude * direction.x / distance;
        body->force.y += magnitude * direction.y / distance;

    }
    /* else if not is not a leaf but contains a body somewhere */
    else if (!node->isLeaf && node->hasBody) {
        distance = sqrt(pow(body->position.x - node->center_of_mass.x, 2) +
                        pow(body->position.y - node->center_of_mass.y, 2));

        /* if distance is too small */
        if (distance < MIN_DISTANCE_APP) {
            distance = MIN_DISTANCE_APP;
        }

        /* if distance is far away */
        if (distance > far) {
            magnitude = (G * body->mass * node->mass) / pow(distance, 2);
            direction.x = node->center_of_mass.x - body->position.x;
            direction.y = node->center_of_mass.y - body->position.y;

            body->force.x += magnitude * direction.x / distance;
            body->force.y += magnitude * direction.y / distance;
        }
        /* if distance is close, recursively look in the children */
        else {
            calculateForces(body, node->ne);
            calculateForces(body, node->nw);
            calculateForces(body, node->sw);
            calculateForces(body, node->se);
        }
    }
}

/* move the bodies according to the forces */
void moveBodies(Body *body) {
    Point delta_v;
    Point delta_p;

    delta_v.x = (body->force.x / body->mass) * DT;
    delta_v.y = (body->force.y / body->mass) * DT;

    delta_p.x = (body->velocity.x + delta_v.x / 2) * DT;
    delta_p.y = (body->velocity.y + delta_v.y / 2) * DT;

    body->velocity.x += delta_v.x;
    body->velocity.y += delta_v.y;

    body->position.x += delta_p.x;
    body->position.y += delta_p.y;

    /* keep bodies inside the space */
    #ifdef BOUNDS
    if (body->position.x > total_size) {
        body->velocity.x *= -1;
        body->position.x = total_size;
    }
    if (body->position.x < 0) {
        body->velocity.x *= -1;
        body->position.x = 0;
    }
    if (body->position.y > total_size) {
        body->velocity.y *= -1;
        body->position.y = total_size;
    } 
    if (body->position.y < 0) {
        body->velocity.y *= -1;
        body->position.y = 0;
    }
    #endif

    /* important to resent the forces */
    body->force.x = 0;
    body->force.y = 0;
}

/* calculate the center of mass */
void calculateCenterOfMasses(Node *node) {
    Point pom = {.x = 0, .y = 0};
    /* if node is a leaf and contains a body */
    if (node->isLeaf && node->hasBody) {
        node->center_of_mass.x = node->body.position.x;
        node->center_of_mass.y = node->body.position.y;
    }
    /* if node is not a leaf but contains a body */
    else if (!node->isLeaf && node->hasBody) {
        /* calculate center of mass for all child nodes */
        calculateCenterOfMasses(node->ne);
        calculateCenterOfMasses(node->nw);
        calculateCenterOfMasses(node->sw);
        calculateCenterOfMasses(node->se);

        /* summarize child nodes */
        pom.x = (node->ne->center_of_mass.x * node->ne->mass) + 
                (node->nw->center_of_mass.x * node->nw->mass) + 
                (node->sw->center_of_mass.x * node->sw->mass) + 
                (node->se->center_of_mass.x * node->se->mass);
        pom.y = (node->ne->center_of_mass.y * node->ne->mass) + 
                (node->nw->center_of_mass.y * node->nw->mass) + 
                (node->sw->center_of_mass.y * node->sw->mass) + 
                (node->se->center_of_mass.y * node->se->mass);

        /* set center of mass for parent node */
        node->center_of_mass.x = pom.x / node->mass;
        node->center_of_mass.y = pom.y / node->mass;
    }
}

/* construct a quad tree */
void constructQuadTree(Body *body, Node *node) {
    node->mass += body->mass;/* increase the mass of the tree */

    if (node->bodyIsHere && node->hasBody) {
        /* check if the two bodies occupy the same point */
        if ((node->body.position.x == body->position.x) && (node->body.position.y == body->position.y)) {
            fprintf(stderr, "Error: particles x:%lf y:%lf occupy the same point\n", body->position.x, body->position.y);
            exit(-1);
        }

        node->isLeaf = false;
        node->bodyIsHere = false;
        initializeChildNodes(node);

        Node *child_node_for_existing_body = getBodyQuadrantInNode(node->body.position, node);
        constructQuadTree(&(node->body), child_node_for_existing_body);

        Node *child_node_for_new_body = getBodyQuadrantInNode(body->position, node);
        constructQuadTree(body, child_node_for_new_body);
    }
    else if (!node->bodyIsHere && !node->isLeaf) {
        node->hasBody = true;
        Node *child_node = getBodyQuadrantInNode(body->position, node);
        constructQuadTree(body, child_node);
    }
    else if (!node->bodyIsHere && !node->hasBody && node->isLeaf) {
        node->hasBody = true;
        node->bodyIsHere = true;
        node->body = *body;
    }
}

/* returns the correct child quadrant that contains a point */
Node *getBodyQuadrantInNode(Point position, Node *node) {
    /* if x values are on left side of middle*/
    if (position.x <= node->position.x + (node->size / 2)) {

        /* if y values are on lower side of middle */
        if (position.y <= node->position.y + (node->size / 2)) {
            return node->sw;
        }
        /* if y values are on upper side of middle */
        else {
            return node->nw;
        }
    }
    /* if x values are on right side of middle */
    else {

        /* if y values are on lower side of middle */
        if (position.y <= node->position.y + (node->size / 2)) {
            return node->se;
        }
        /* if y values are on upper side of middle */
        else {
            return node->ne;
        }
    }
}

/* check if worker is assigned to this node */
bool isCorrectSubtree(long id_worker, int id_node) {
    int i;/* for iteration */

    /* check if id_worker+(num_workers * i)==id_node : 0<i<MAX_WORKERS */
    for (i = 0; i < MAX_WORKERS; i++) {
        if (id_worker + (num_workers * i) > 15) {
            return false;
        }
        if (id_worker + (num_workers * i) == id_node) {
            return true;
        }
    }
    return false;
}

/* check if body is far away and can be ignored */
bool checkFarAway(Point position) {
    if (position.x < 0 || position.y < 0 || position.x > total_size || position.y > total_size) {
        return false;
    }
    return true;
}

/* needed for calculating the top two layers of the tree */
void summarizeWeights() {
    /* check if ne child contains a body */
    if (root->ne->ne->hasBody || root->ne->nw->hasBody || root->ne->sw->hasBody || root->ne->se->hasBody) {
        root->ne->hasBody = true;
    }

    /* check if nw child contains a body */
    if (root->nw->ne->hasBody || root->nw->nw->hasBody || root->nw->sw->hasBody || root->nw->se->hasBody) {
        root->nw->hasBody = true;
    }

    /* check if sw child contains a body */
    if (root->sw->ne->hasBody || root->sw->nw->hasBody || root->sw->sw->hasBody || root->sw->se->hasBody) {
        root->sw->hasBody = true;
    }

    /* check if se child contains a body */
    if (root->se->ne->hasBody || root->se->nw->hasBody || root->se->sw->hasBody || root->se->se->hasBody) {
        root->se->hasBody = true;
    }

    /* check if root contains a body */
    if (root->ne->hasBody || root->nw->hasBody || root->sw->hasBody || root->se->hasBody) {
        root->hasBody = true;
    }

    /* summarize ne weight */
    root->ne->mass += root->ne->ne->mass;
    root->ne->mass += root->ne->nw->mass;
    root->ne->mass += root->ne->sw->mass;
    root->ne->mass += root->ne->se->mass;

    /* summarize nw weight */
    root->nw->mass += root->nw->ne->mass;
    root->nw->mass += root->nw->nw->mass;
    root->nw->mass += root->nw->sw->mass;
    root->nw->mass += root->nw->se->mass;

    /* summarize sw weight */
    root->sw->mass += root->sw->ne->mass;
    root->sw->mass += root->sw->nw->mass;
    root->sw->mass += root->sw->sw->mass;
    root->sw->mass += root->sw->se->mass;

    /* summarize se weight */
    root->se->mass += root->se->ne->mass;
    root->se->mass += root->se->nw->mass;
    root->se->mass += root->se->sw->mass;
    root->se->mass += root->se->se->mass;

    /* summarize root weight */
    root->mass += root->ne->mass;
    root->mass += root->nw->mass;
    root->mass += root->sw->mass;
    root->mass += root->se->mass;
}

/* initialize bodies */
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

    /* generate */
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

/* initialize the tree */
void initializeTree() {
    /* allocate space on the heap */
    root = (Node *)malloc(sizeof(Node));
    root->size = total_size;
    root->mass = 0;
    root->center_of_mass.x = 0;
    root->center_of_mass.y = 0;
    root->position.x = 0;
    root->position.y = 0;
    root->hasBody = false;
    root->bodyIsHere = false;
    root->isLeaf = false;/* leaf to false since we initialize children now */

    /* initialize 1 order children */
    initializeChildNodes(root);
    root->ne->isLeaf = false;
    root->nw->isLeaf = false;
    root->sw->isLeaf = false;
    root->se->isLeaf = false;

    /* initialize 2 order children and set their id */
    initializeChildNodes(root->ne);
    initializeChildNodes(root->nw);
    initializeChildNodes(root->sw);
    initializeChildNodes(root->se);
    root->ne->ne->id = 0;
    root->ne->nw->id = 1;
    root->ne->sw->id = 2;
    root->ne->se->id = 3;

    root->nw->ne->id = 4;
    root->nw->nw->id = 5;
    root->nw->sw->id = 6;
    root->nw->se->id = 7;

    root->sw->ne->id = 8;
    root->sw->nw->id = 9;
    root->sw->sw->id = 10;
    root->sw->se->id = 11;

    root->se->ne->id = 12;
    root->se->nw->id = 13;
    root->se->sw->id = 14;
    root->se->se->id = 15;

    childNodes[0] = root->ne->ne;
    childNodes[1] = root->ne->nw;
    childNodes[2] = root->ne->sw;
    childNodes[3] = root->ne->se;

    childNodes[4] = root->nw->ne;
    childNodes[5] = root->nw->nw;
    childNodes[6] = root->nw->sw;
    childNodes[7] = root->nw->se;

    childNodes[8] = root->sw->ne;
    childNodes[9] = root->sw->nw;
    childNodes[10] = root->sw->sw;
    childNodes[11] = root->sw->se;

    childNodes[12] = root->se->ne;
    childNodes[13] = root->se->nw;
    childNodes[14] = root->se->sw;
    childNodes[15] = root->se->se;
}

/* initialize the child nodes of a parent node */
void initializeChildNodes(Node *node) {
    /* allocate space on the heap */
    node->ne = (Node *)malloc(sizeof(Node));
    node->nw = (Node *)malloc(sizeof(Node));
    node->sw = (Node *)malloc(sizeof(Node));
    node->se = (Node *)malloc(sizeof(Node));

    /* lenght of children is half of parent */
    double child_size = node->size / 2;
    node->ne->size = child_size;
    node->nw->size = child_size;
    node->sw->size = child_size;
    node->se->size = child_size;

    /* set values for ne child */
    node->ne->mass = 0;
    node->ne->center_of_mass.x = 0;
    node->ne->center_of_mass.y = 0;
    node->ne->position.x = node->position.x + child_size;
    node->ne->position.y = node->position.y + child_size;
    node->ne->hasBody = false;
    node->ne->bodyIsHere = false;
    node->ne->isLeaf = true;
    node->ne->id = -1;

    /* set values for nw child */
    node->nw->mass = 0;
    node->nw->center_of_mass.x = 0;
    node->nw->center_of_mass.y = 0;
    node->nw->position.x = node->position.x;
    node->nw->position.y = node->position.y + child_size;
    node->nw->hasBody = false;
    node->nw->bodyIsHere = false;
    node->nw->isLeaf = true;
    node->nw->id = -1;

    /* set values for sw child */
    node->sw->mass = 0;
    node->sw->center_of_mass.x = 0;
    node->sw->center_of_mass.y = 0;
    node->sw->position.x = node->position.x;
    node->sw->position.y = node->position.y;
    node->sw->hasBody = false;
    node->sw->bodyIsHere = false;
    node->sw->isLeaf = true;
    node->sw->id = -1;

    /* set values for se child */
    node->se->mass = 0;
    node->se->center_of_mass.x = 0;
    node->se->center_of_mass.y = 0;
    node->se->position.x = node->position.x + child_size;
    node->se->position.y = node->position.y;
    node->se->hasBody = false;
    node->se->bodyIsHere = false;
    node->se->isLeaf = true;
    node->se->id = -1;
}

/* free nodes recursively */
void freeNodes(Node* node) {
    if (!node->isLeaf) {
        freeNodes(node->ne);
        freeNodes(node->nw);
        freeNodes(node->sw);
        freeNodes(node->se);
    }
    free(node);
}

/* set all nodes to not contain a body */
void hasBodyRemove(Node* node) {
    if (!node->isLeaf) {
        if (node->ne->hasBody) {
            hasBodyRemove(node->ne);
        }
        if (node->nw->hasBody) {
            hasBodyRemove(node->nw);
        }
        if (node->sw->hasBody) {
            hasBodyRemove(node->sw);
        }
        if (node->se->hasBody) {
            hasBodyRemove(node->se);
        }
    }
    node->hasBody = false;
    node->bodyIsHere = false;
    node->mass = 0;
}

/* free all unused nodes */
void freeWasteNodes(Node* node) {
    if (!node->hasBody && !node->isLeaf) {
            freeNodes(node->ne);
            freeNodes(node->nw);
            freeNodes(node->sw);
            freeNodes(node->se);
            node->isLeaf = true;
    }
    else if (node->hasBody && !node->isLeaf) {
        freeWasteNodes(node->ne);
        freeWasteNodes(node->nw);
        freeWasteNodes(node->sw);
        freeWasteNodes(node->se);
    }
}

/* barrier for waiting for workers */
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

/* start the timer */
void start_clock() {
    st_time = times(&st_cpu);
}

/* end the timer */
void end_clock() {
    en_time = times(&en_cpu);
    printf("%g\n", (double)(en_time - st_time)/sysconf(_SC_CLK_TCK));
}

/* start the timer */
clock_t start_clock_lap() {
    return times(&st_cpu);
}

/* end the timer */
clock_t end_clock_lap(clock_t start_time) {
    return times(&en_cpu) - start_time;
}