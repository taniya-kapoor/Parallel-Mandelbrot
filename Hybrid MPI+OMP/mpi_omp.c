/*
 SU Project -- Taniya -- Hybrid MPI + OMP
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

/* Bounds of the Mandelbrot set */
#define X_MIN -1.78
#define X_MAX 0.78
#define Y_MIN -0.961
#define Y_MAX 0.961

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
  * path = "q3-1.ppm";

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
    fprintf (f, "%c%c%c", 2*c, 3*c, 4*c); /* Monochrome weight */
  }
  fclose (f);
}

static void Compute (Image * im, int nb_iter, double x_min, double x_max, double y_min, double y_max) {
    
    //int pos = 0;
    
    int l, c, i = 0;
    
    double dx = (x_max - x_min) / im -> nb_columns, dy = (y_max - y_min) / im -> nb_rows; /* Discretization */
    
    double x, y, a, b, tmp;
#pragma omp parallel for private(l,c,a,b,x,y,i,tmp) schedule(static)
    for (l = 0; l < im -> nb_rows; l ++) {
        
        for (c = 0; c < im -> nb_columns; c ++) {
            
            /* Computation at each point of the image */
            
            a = x_min + c * dx;
            b = y_max - l * dy;
            x = 0; y = 0;
            i=0;
            while (i < nb_iter) {
                tmp = x;
                x = x * x - y * y + a;
                y = 2 * tmp * y + b;
                if (x * x + y * y > 4) /* Divergence ! */
                    break;
                else
                    i++;
            }
            
            im -> pixels [im -> nb_columns*l+c] = (double) i / nb_iter * 255;
        }
    }
}

int main (int argc, char * * argv) {
    
    int nbProc, myRank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nbProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Status status;
    
    int nb_iter, width, height;
    double x_min, x_max, y_min, y_max;
    char * path;
    Image imloc;
    Image im;
    
    analyzis(argc, argv, & nb_iter, & x_min, & x_max, & y_min, & y_max, & width, & height, & path);
    
    int locheight = height/nbProc;
    
    initialization (& imloc, width, locheight);
    initialization (& im, width, height);
    
    double arr[nbProc+1];
    int i;
    double d = (y_max-y_min)/nbProc;
    
    double min = y_min;
    
    for (i = 0; i<=nbProc; i++){
        arr[i] = min;
        min = min+d;
    }
    
    double starttime, endtime, time_parallel, time_serial, speed_up, final_speed_up;
    starttime = MPI_Wtime();
    
    Compute (& imloc, nb_iter, x_min, x_max, arr[nbProc-myRank-1], arr[nbProc-myRank]);
    
    MPI_Gather(& imloc.pixels[0],width*locheight,MPI_CHAR,& im.pixels[0],width*locheight,MPI_CHAR,0,MPI_COMM_WORLD);
    
    if (myRank == 0) {
        save (& im, path);
    }
    
    endtime = MPI_Wtime();
    printf("That took %f seconds in processor %d\n",endtime-starttime, myRank);
    
    time_serial = 0.242479876449473;//from g5k serial run
    time_parallel = endtime - starttime;
    speed_up = time_serial/time_parallel;
    printf("time taken in parallel is: %f ",time_parallel);
    printf("The speed up obtained is: %f ",speed_up);
    
    MPI_Finalize();
    return 0 ;
}

