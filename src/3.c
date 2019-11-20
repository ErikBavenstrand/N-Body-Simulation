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

/* defines */
#define MAX_BODIES 1000000
#define MAX_STEPS 1000000
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
} Node;

/* global variables */
Body *bodies;
Node *root;

int num_bodies;
int num_steps;
int far;
int total_size;
int seed;

static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;

void calculateForces(Body *body, Node *node);
void moveBodies(Body *body);
void initializeTree();
void initializeChildNodes(Node *node);
void constructQuadTree(Body *body, Node *node);
Node *getBodyQuadrantInNode(Point position, Node *node);
void calculateCenterOfMasses(Node *node);
bool checkFarAway(Point position);
void freeNodes(Node* node);
void initializeBodies();
void initializeSolarsystem();
void start_clock();
void end_clock();
void freeWasteNodes(Node* node);
void hasBodyRemove(Node* node);


int main(int argc, char *argv[]) {
    int i, j;/* for iteration */

    /* create folder for output */
    #ifdef VISUALIZE
    struct stat st = {0};
    if (stat("./outfiles", &st) == -1) {
        mkdir("./outfiles", 0700);
    }
    FILE *fp = fopen("./outfiles/3.out","w");
    #endif

    /* check if arg count is correct */
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_bodies num_steps far size seed\n", argv[0]);
        fflush(stderr);
        exit(-1);
    }

    /* set values from input */
    num_bodies = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    far = atoi(argv[3]);
    total_size = atoi(argv[4]);
    seed = atoi(argv[5]);

    /* sanity check all input values */
    if (num_bodies > MAX_BODIES || num_bodies <= 0) num_bodies = MAX_BODIES;
    if (num_steps > MAX_STEPS || num_steps <= 0) num_steps = MAX_STEPS;
    if (far > MAX_FAR || far < 0) far = MAX_FAR;
    if (total_size > MAX_SIZE || total_size <= 0) total_size = MAX_SIZE;

    /* initialize bodies */
    #ifdef SOLAR_SYSTEM
    initializeSolarsystem();
    #else
    initializeBodies();
    #endif

    /* print initial lines to file */
    #ifdef VISUALIZE
    fprintf(fp, "%d %d\n", num_bodies, 1);
    #endif

    start_clock();/* start time */
    initializeTree();/* initialize the tree */
    for (i = 0; i < num_steps; i++) {
        /* construct the tree */
        for (j = 0; j < num_bodies; j++) {
            /* check if body is inside given area */
            if (checkFarAway(bodies[j].position)) {
                constructQuadTree(&bodies[j], root);
            }
        }

        freeWasteNodes(root);

        calculateCenterOfMasses(root);/* calculate the center of mass for all nodes */

        /* calculate the forces for all bodies */
        for (j = 0; j < num_bodies; j++) {
            calculateForces(&bodies[j], root);
        }

        /* move all bodies according to forces */
        for (j = 0; j < num_bodies; j++) {
            moveBodies(&bodies[j]);
        }

        /* print positions to file */
        #ifdef VISUALIZE
        for (j = 0; j < num_bodies; j++) {
            fprintf(fp, "%lf %lf\n", bodies[j].position.x / 1000, bodies[j].position.y / 1000);
            fflush(fp);
        }
        #endif

        /* free allocated memory from the heap */
        hasBodyRemove(root);
    }
    end_clock();/* end time */

    #ifdef VISUALIZE
    fclose(fp);
    #endif

    free(bodies);

    return 0;
}

/* calculate forces for all bodies */
void calculateForces(Body *body, Node *node) {
    double distance;
    double magnitude;
    Point direction;

    /* if the node is a leaf and contains a body */
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
    /* if body is not a leaf but has a body somewhere inside it*/
    else if (!node->isLeaf && node->hasBody) {
        distance = sqrt(pow(body->position.x - node->center_of_mass.x, 2) +
                        pow(body->position.y - node->center_of_mass.y, 2));

        /* if distance is too small */
        if (distance < MIN_DISTANCE_APP) {
            distance = MIN_DISTANCE_APP;
        }

        /* if distance is far away, use center of mass for that node */
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

/* initialize the tree */
void initializeTree() {
    root = (Node *)malloc(sizeof(Node));
    root->size = total_size;
    root->mass = 0;
    root->position.x = 0;
    root->position.y = 0;
    root->hasBody = false;
    root->bodyIsHere = false;
    root->isLeaf = true;
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

    /* set values for nw child */
    node->nw->mass = 0;
    node->nw->center_of_mass.x = 0;
    node->nw->center_of_mass.y = 0;
    node->nw->position.x = node->position.x;
    node->nw->position.y = node->position.y + child_size;
    node->nw->hasBody = false;
    node->nw->bodyIsHere = false;
    node->nw->isLeaf = true;

    /* set values for sw child */
    node->sw->mass = 0;
    node->sw->center_of_mass.x = 0;
    node->sw->center_of_mass.y = 0;
    node->sw->position.x = node->position.x;
    node->sw->position.y = node->position.y;
    node->sw->hasBody = false;
    node->sw->bodyIsHere = false;
    node->sw->isLeaf = true;

    /* set values for se child */
    node->se->mass = 0;
    node->se->center_of_mass.x = 0;
    node->se->center_of_mass.y = 0;
    node->se->position.x = node->position.x + child_size;
    node->se->position.y = node->position.y;
    node->se->hasBody = false;
    node->se->bodyIsHere = false;
    node->se->isLeaf = true;
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

/* return the correct child from a position and a parent */
Node *getBodyQuadrantInNode(Point position, Node *node) {
    if (position.x <= node->position.x + (node->size / 2)) {
        if (position.y <= node->position.y + (node->size / 2)) {
            return node->sw;
        }
        else {
            return node->nw;
        }
    }
    else {
        if (position.y <= node->position.y + (node->size / 2)) {
            return node->se;
        }
        else {
            return node->ne;
        }
    }
}

/* calculate the center of mass */
void calculateCenterOfMasses(Node *node) {
    Point pom = {.x = 0, .y = 0};
    /* if node is a leaf and has a body */
    if (node->isLeaf && node->hasBody) {
        node->center_of_mass.x = node->body.position.x;
        node->center_of_mass.y = node->body.position.y;
    }
    /* if body is not a leaf and has a body somewhere */
    else if (!node->isLeaf && node->hasBody) {
        
        /* calculate the child nodes */
        calculateCenterOfMasses(node->ne);
        calculateCenterOfMasses(node->nw);
        calculateCenterOfMasses(node->sw);
        calculateCenterOfMasses(node->se);

        /* summarize the child nodes */
        pom.x = (node->ne->center_of_mass.x * node->ne->mass) + 
                (node->nw->center_of_mass.x * node->nw->mass) + 
                (node->sw->center_of_mass.x * node->sw->mass) + 
                (node->se->center_of_mass.x * node->se->mass);
        pom.y = (node->ne->center_of_mass.y * node->ne->mass) + 
                (node->nw->center_of_mass.y * node->nw->mass) + 
                (node->sw->center_of_mass.y * node->sw->mass) + 
                (node->se->center_of_mass.y * node->se->mass);
        
        /* set the center of mass */
        node->center_of_mass.x = pom.x / node->mass;
        node->center_of_mass.y = pom.y / node->mass;
    }
}

/* check if body is far away and can be ignored */
bool checkFarAway(Point position) {
    if (position.x < 0 || position.y < 0 || position.x > total_size || position.y > total_size) {
        return false;
    }
    return true;
}

/* free all allocated nodes recursivly */
void freeNodes(Node* node) {
    if (!node->isLeaf) {
        freeNodes(node->ne);
        freeNodes(node->nw);
        freeNodes(node->sw);
        freeNodes(node->se);
    }
    free(node);
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

/* start the timer */
void start_clock() {
    st_time = times(&st_cpu);
}

/* end the timer */
void end_clock() {
    en_time = times(&en_cpu);
    printf("%g\n", (double)(en_time - st_time)/sysconf(_SC_CLK_TCK));
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