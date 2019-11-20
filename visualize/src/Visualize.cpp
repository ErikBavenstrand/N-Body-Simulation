/* 
    Erik Båvenstrand, 2019
*/

/* included in files */
#include <GLFW/glfw3.h>

/* standard imports */
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <sys/time.h>
#include <vector>

/* defines */
#define eps 0.1 /*  */
#define SCALE 800 /* zoom of window */
#define FPS 1000 /* speed of animation */
#define MIN_SIZE 100 /* minimum window size */

double read_timer();
static void error_callback(int error, const char* description);
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void fixRenderOnMac(GLFWwindow* window);

/* struct for positions of particles */
struct particle_t { 
    float x, y; 
};

int main(int argc, char *argv[]) {
    /* input file */
    const char *filename = argv[1];
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        printf("failed to find %s\n", filename);
        return 1;
    }
    
    int n;/* number of particles */
    float size;/* size of area they exist within */
    fscanf(f, "%d%g", &n, &size);/* get values from first line in file */
	
    /* read entire file and store in vector */
    particle_t p;
    std::vector<particle_t> particles;
    while(fscanf(f, "%g%g", &p.x, &p.y) == 2) {
        particles.push_back(p);
    }
    fclose(f);
	
    /* total number of positions for each particle */
    int nframes = particles.size()/n;
    if(nframes == 0) {
        return 2;
    }
    
    /* calculate window size */
    int window_size = (int)((size+2*eps)*SCALE);
    window_size = window_size > MIN_SIZE ? window_size : MIN_SIZE;
	
    /* initialize and create window */
    GLFWwindow* window;
    glfwSetErrorCallback(error_callback);
    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }
    window = glfwCreateWindow(window_size, window_size, "Visualizer", NULL, NULL);

    /* check if window creation was a success */
    if (!window) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    /* bind glfw to the currently opened OpenGL window */
    glfwMakeContextCurrent(window);

    /* set function as key input handler */
    glfwSetKeyCallback(window, key_callback);

    glEnable(GL_POINT_SMOOTH);/* draw points with proper filtering */
    glEnable(GL_BLEND);/* blend the computed fragment color values with the values in the color buffers */
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);/* specifies the index of the draw buffer for which to set the blend function */
	
    glPointSize(8);/* size of the points */
    glClearColor(1, 1, 1, 1);/* color of the screen when clear function is called. rgba */
   
    glViewport(0, 0, window_size, window_size);/* set viewport size */
    glMatrixMode(GL_PROJECTION);/* applies subsequent matrix operations to the projection matrix stack */
    glLoadIdentity();/* replace the current matrix with the identity matrix */

    /* while window close flag is not up */
    while (!glfwWindowShouldClose(window)) {
        /* get current size of window and adjust the viewport */
        float ratio;
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        ratio = width / (float) height;
        glViewport(0, 0, width, height);
        

        glClear(GL_COLOR_BUFFER_BIT);/* clear the window to previously set color */
        glLoadIdentity();
        glOrtho(-eps, size+eps, -eps, size+eps, 1.f, -1.f);

        /* draw the edges in almost black */
        glColor3f( 0.75, 0.75, 0.75 );
        glBegin(GL_LINE_LOOP);
        glBegin( GL_LINE_LOOP );
        glVertex2d( 0, 0 );
        glVertex2d( size, 0 );
        glVertex2d( size, size );
        glVertex2d( 0, size );
        glEnd( );

        /* calculate current frame */
        int iframe = (int)(read_timer()*FPS) % nframes;
        particle_t *p = &particles[iframe*n];
		
        /* draw each of the positions for this frame */
        glColor3f( 0, 0, 0 );
        glBegin( GL_POINTS );
        for( int i = 0; i < n; i++ )
            glVertex2fv( &p[i].x );
        glEnd( );

        /* swap the currently displayed buffer to the new one  */
        glfwSwapBuffers(window);
        glfwPollEvents();

        fixRenderOnMac(window);/* mac has an issue it seems */
    }
    /* exit */
    glfwDestroyWindow(window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
}

double read_timer() {
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if(!initialized) {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday(&end, NULL);
	
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

static void error_callback(int error, const char* description) {
    fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

void fixRenderOnMac(GLFWwindow* window) {
#ifdef __APPLE__
    static bool macMoved = false;
    if(!macMoved) {
        int x, y;
        glfwGetWindowPos(window, &x, &y);
        glfwSetWindowPos(window, ++x, y);
        macMoved = true;
    }
#endif
}