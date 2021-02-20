/*
 SU Project -- Taniya -- Cuda
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <cuda_runtime.h>

/* Bounds of the Mandelbrot set */
#define X_MIN -1.78
#define X_MAX 0.78
#define Y_MIN -0.961
#define Y_MAX 0.961

__global__ void pixel_calculation(double dx, double dy, char * pixels, int nb_iter, double x_min, double y_max, int columns);

typedef struct {

  int nb_rows, nb_columns; /* Dimensions */
  char * pixels; /* Linearized matrix of pixels */

} Image;

static void error_options () {

  fprintf (stderr, "Use : ./mandel [options]\n\n");
  fprintf (stderr, "Options \t Meaning \t\t Default val.\n\n");
  fprintf (stderr, "-n \t\t Nb iter. \t\t 100\n");
  fprintf (stderr, "-b \t\t Bounds \t\t -1.78 0.78 -0.961 0.961\n");
  fprintf (stderr, "-d \t\t Dimensions \t\t 1024 768\n");
  fprintf (stderr, "-f \t\t File \t\t /tmp/mandel.ppm\n");
  exit (1);
}

static void analyzis (int argc, char * * argv, int * nb_iter, double * x_min, double * x_max, double * y_min, double * y_max, int * width, int * height, char * * path) {

  const char * opt = "b:d:n:f:" ;
  int c ;

  /* Default values */
  * nb_iter = 100;
  * x_min = X_MIN;
  * x_max = X_MAX;
  * y_min = Y_MIN;
  * y_max = Y_MAX;
  * width = 1024;
  * height = 768;
  * path = "mandel.ppm";

  /* Analysis of arguments */
  while ((c = getopt (argc, argv, opt)) != EOF) {
    
    switch (c) {
      
    case 'b':
      sscanf (optarg, "%lf", x_min);
      sscanf (argv [optind ++], "%lf", x_max);
      sscanf (argv [optind ++], "%lf", y_min);
      sscanf (argv [optind ++], "%lf", y_max);
      break ;
    case 'd': /* width */
      sscanf (optarg, "%d", width);
      sscanf (argv [optind ++], "%d", height);
      break;
    case 'n': /* Number of iterations */
      * nb_iter = atoi (optarg);
      break;
    case 'f': /* Output file */
      * path = optarg;
      break;
    default :
      error_options ();
    };
  }  
}

static void initialization (Image * im, int nb_columns, int nb_rows) {
  im -> nb_rows = nb_rows;
  im -> nb_columns = nb_columns;
  im -> pixels = (char *) malloc (sizeof (char) * nb_rows * nb_columns); /* Space memory allocation */
} 

static void save (const Image * im, const char * path) {
  /* Image saving using the ASCII format'.PPM' */
  unsigned i;
  FILE * f = fopen (path, "w");  
  fprintf (f, "P6\n%d %d\n255\n", im -> nb_columns, im -> nb_rows); 
  for (i = 0; i < im -> nb_columns * im -> nb_rows; i ++) {
    char c = im -> pixels [i];
    fprintf (f, "%c%c%c", c, c, c); /* Monochrome weight */
  }
  fclose (f);
}

static void Compute (Image * im, int nb_iter, double x_min, double x_max, double y_min, double y_max) {
  
  double dx = (x_max - x_min) / im -> nb_columns, dy = (y_max - y_min) / im -> nb_rows; /* Discretization */
    
     int rownum = im->nb_rows, colnum = im-> nb_columns;

    dim3 blocksize(16,16,1); // 16 blocks of 16 threads each
    dim3 nblocks(rownum/16, colnum/16, 1);
    
    char * im_pixels_d;
    cudaMalloc(&im_pixels_d, sizeof(char)*rownum*colnum);
    cudaMemcpy(im_pixels_d, im->pixels , sizeof(char) * rownum * colnum,cudaMemcpyHostToDevice);
    pixel_calculation<<< nblocks, blocksize >>> (dx, dy, im_pixels_d, nb_iter, x_min, y_max, colnum);
    
    cudaMemcpy(im -> pixels, im_pixels_d, sizeof(char)*rownum*colnum,cudaMemcpyDeviceToHost);
    cudaFree(im_pixels_d);
}

__global__ void pixel_calculation(double dx, double dy, char * pixels, int nb_iter, double x_min, double y_max, int colnum)
{
    
    int id_x = blockIdx.x *blockDim.x + threadIdx.x;
    int id_y = blockIdx.y *blockDim.y + threadIdx.y;
    
    double a = x_min + id_y * dx, b = y_max - id_x * dy, x = 0, y = 0;
     int i=0;
      while (i < nb_iter) {
	double tmp = x;
	x = x * x - y * y + a;
	y = 2 * tmp * y + b;
	if (x * x + y * y > 4) /* Divergence ! */
	  break; 
	else
	  i++;
      }
      
      pixels [id_x*colnum+id_y] = (double) i / nb_iter * 255;
    }

int main (int argc, char * * argv) {
  
  int nb_iter, width, height; /* Degree of precision, dimensions of the image */  
  double x_min, x_max, y_min, y_max; /* Bounds of representation */
  char * path; /* File destination */
  Image im;

  analyzis(argc, argv, & nb_iter, & x_min, & x_max, & y_min, & y_max, & width, & height, & path);
  initialization (& im, width, height);
  Compute (& im, nb_iter, x_min, x_max, y_min, y_max);
  save (& im, path);
  return 0 ;
}
