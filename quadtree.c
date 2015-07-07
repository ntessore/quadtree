#include <stdlib.h>
#include <stdio.h>
#include <math.h>


//--------------------------------------
// basic definitions
//--------------------------------------

#define WIDTH   20
#define HEIGHT  20
#define N       10
#define THRESH  1.0

#define LENS_X  11.23
#define LENS_Y  9.87
#define LENS_B  6.34
#define LENS_PA 34.56
#define LENS_Q  0.78


//--------------------------------------
// auxiliary stuff
//--------------------------------------

// a simple two-dimensional point
typedef struct
{
    double x;
    double y;
} double2;

// a simple SIE lens
void lens(double2* p)
{
    double x, y, r;
    double ax, ay;
    
    // position
    double x0 = LENS_X;
    double y0 = LENS_Y;
    
    // scale radius
    double b = LENS_B;
    
    // axis ratio
    double q = LENS_Q;
    
    // position angle
    double c = cos(LENS_PA*M_PI/180);
    double s = sin(LENS_PA*M_PI/180);
    
    // in central & rotated coordinate system
    x = (p->x - x0)*c - (p->y - y0)*s;
    y = (p->x - x0)*s + (p->y - y0)*c;
    
    // elliptical radius
    r = sqrt(q*q*x*x + y*y);
    
    // deflection angle
    ax = b*sqrt(q)/sqrt(1-q*q)*atan(x*sqrt(1-q*q)/r);
    ay = b*sqrt(q)/sqrt(1-q*q)*atanh(y*sqrt(1-q*q)/r);
    
    // apply deflection
    p->x -= ax*c + ay*s;
    p->y -= ay*c - ax*s;
}


//--------------------------------------
// the quadtree
// -------------------------------------

// forward declaration
typedef struct quadtree quadtree;

// the quadtree struct with position,
// dimensions, child nodes, and points
// to be sorted into the tree
struct quadtree
{
    double      x, y;
    double      w, h;
    
    quadtree*   child;
    
    size_t      size;
    double2*    points;
};

// add a point to the quadtree, takes
// care of enlarging the points array
void quadtree_add(quadtree* q, double2* p)
{
    // allocate space for N more points if full
    if(q->size%N == 0)
    {
        void* tmp = realloc(q->points, (q->size + N)*sizeof(double2));
        if(tmp == NULL)
        {
            perror("could not allocate points");
            exit(EXIT_FAILURE);
        }
        q->points = tmp;
    }
    
    // add point
    q->points[q->size].x = p->x;
    q->points[q->size].y = p->y;
    q->size += 1;
}

// recursively refine the tree node
void quadtree_refine(quadtree* q)
{
    // indices for child node
    int i, j;
    
    // refine if number of points exceeds threshold
    if(q->size > THRESH*N*N)
    {
        // create child nodes
        q->child = malloc(4*sizeof(quadtree));
        if(q->child == NULL)
        {
            perror("could not allocate children");
            exit(EXIT_FAILURE);
        }
        
        // initialise child nodes
        for(j = 0; j < 2; ++j)
        for(i = 0; i < 2; ++i)
        {
            q->child[j*2 + i].x      = q->x + (2*i - 1)*0.25*q->w;
            q->child[j*2 + i].y      = q->y + (2*j - 1)*0.25*q->h;
            q->child[j*2 + i].w      = 0.5*q->w;
            q->child[j*2 + i].h      = 0.5*q->h;
            q->child[j*2 + i].child  = NULL;
            q->child[j*2 + i].size   = 0;
            q->child[j*2 + i].points = NULL;
        }
        
        // sort points into child nodes
        for(int k = 0; k < q->size; ++k)
        {
            // get quadrant of point
            i = q->points[k].x > q->x;
            j = q->points[k].y > q->y;
            
            // add point to quadrant
            quadtree_add(&q->child[j*2 + i], &q->points[k]);
        }
        
        // free the points, they have moved to children
        free(q->points);
        q->points = NULL;
        q->size = 0;
        
        // refine child nodes
        for(int k = 0; k < 4; ++k)
            quadtree_refine(q->child + k);
    }
}

// apply a function to leaf nodes of the tree
void quadtree_apply_leaves(quadtree* q, void (*f)(quadtree*))
{
    // if node has children, recurse
    if(q->child)
    {
        for(int k = 0; k < 4; ++k)
            quadtree_apply_leaves(q->child + k, f);
    }
    else
    {
        // apply to this leaf node
        f(q);
    }
}

// recursively free the tree node and its children
void quadtree_free(quadtree* q)
{
    if(q->child)
    {
        for(int i = 0; i < 4; ++i)
            quadtree_free(q->child + i);
        free(q->child);
    }
}


//--------------------------------------
// driver program
//--------------------------------------

// function to print out a source grid cell
void grid_print(quadtree* q)
{
    static int counter = 0;
    printf("%10d%10g%10g%10g%10g%10zu\n",
           ++counter, q->x, q->y, q->w, q->h, q->size);
}

// create and output quadtree
int main(int argc, char* argv[])
{
    // a point for lensing
    double2 point;
    
    // indices
    int i, j;
    
    // the tree
    quadtree* q = malloc(WIDTH*HEIGHT*sizeof(quadtree));
    if(q == NULL)
    {
        perror("could not allocate quadtree");
        exit(EXIT_FAILURE);
    }
    
    // create tree roots at integer points 1...WIDTH, 1...HEIGHT
    for(int n = 0; n < WIDTH*HEIGHT; ++n)
    {
        q[n].x      = n%WIDTH + 1;
        q[n].y      = n/WIDTH + 1;
        q[n].w      = 1;
        q[n].h      = 1;
        q[n].child  = NULL;
        q[n].size   = 0;
        q[n].points = NULL;
    }
    
    // sample points in each pixel
    for(int n = 0; n < WIDTH*HEIGHT; ++n)
    {
        // N*N grid of points
        for(int k = 0; k < N*N; ++k)
        {
            // grid point
            point.x = (n%WIDTH) + 0.5 + ((k%N) + 0.5)/N;
            point.y = (n/WIDTH) + 0.5 + ((k/N) + 0.5)/N;
            
            // lens the point
            lens(&point);
            
            // skip out-of-range points
            if(point.x < 0.5 || point.x >= WIDTH + 0.5 || point.y < 0.5 || point.y >= HEIGHT + 0.5)
                continue;
            
            // indices for lensed point
            i = point.x - 0.5;
            j = point.y - 0.5;
            
            // sort point into tree root
            quadtree_add(&q[j*WIDTH + i], &point);
        }
    }
    
    // at this point, the WIDTH*HEIGHT tree nodes are filled with the lensed
    // grid points, and now need to be refined into a quadtree
    
    // recursively refine the grid
    for(int n = 0; n < WIDTH*HEIGHT; ++n)
        quadtree_refine(q + n);
    
    // the quadtree was built and the source plane grid is done
    
    // recursively output each grid cell (i.e. leaf nodes of quadtree)
    for(int n = 0; n < WIDTH*HEIGHT; ++n)
        quadtree_apply_leaves(q + n, grid_print);
    
    // free the tree
    for(int n = 0; n < WIDTH*HEIGHT; ++n)
        quadtree_free(q + n);
    free(q);
    
    return EXIT_SUCCESS;
}
