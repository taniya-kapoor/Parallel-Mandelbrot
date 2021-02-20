/*
 SU project Taniya -- MPI_Pack/Unpack Routine
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>

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
  fprintf (stderr, "-f \t\t File \t\t mandel.ppm\n");
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
  
  int pos = 0;

  int l, c, i = 0;
  
  double dx = (x_max - x_min) / im -> nb_columns, dy = (y_max - y_min) / im -> nb_rows; /* Discretization */

  for (l = 0; l < im -> nb_rows; l ++) {
    
    for (c = 0; c < im -> nb_columns; c ++) {

      /* Computation at each point of the image */

      double a = x_min + c * dx, b = y_max - l * dy, x = 0, y = 0;      
      i=0;
      while (i < nb_iter) {
	double tmp = x;
	x = x * x - y * y + a;
	y = 2 * tmp * y + b;
	if (x * x + y * y > 4) /* Divergence ! */
	  break; 
	else
	  i++;
      }
      
      im -> pixels [pos ++] = (double) i / nb_iter * 255;    
    }
  }
}

int main (int argc, char * * argv) {
  
  int nb_iter, width, height; /* Degree of precision, dimensions of the image */  
  double x_min, x_max, y_min, y_max; /* Bounds of representation */
  char * path; /* File destination */
  Image im, imloc;
    
  MPI_Status status;
  int id = 0;
  char buffer[1<<20]; // power raised for large memory call
  int *tag
    
  int size, rank, local_size, locheight, lsize;
  double dy, loc_ymin, loc_ymax;
  int l;
    
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  analyzis(argc, argv, & nb_iter, & x_min, & x_max, & y_min, & y_max, & width, & height, & path);
  initialization (& im, width, height);

  locheight = 1;
  dy = (y_max-y_min)/height;
  local_size = (int)(width*locheight);
  initialization (& imloc, width, locheight);

  tag = malloc(sizeof(int)*size);
  for(int i=0; i<size; i++)
    tag[i] = (i+1)*100;
  
  lsize = height/size;
  
    
    double starttime, endtime, time_parallel, time_serial, speed_up, final_speed_up;
    
    starttime = MPI_Wtime();

  for(int jj=0; jj<lsize; jj++)
    {
      if(rank!=0)
	{
	  loc_ymax = y_max-dy*(rank+size*jj);
	  loc_ymin = loc_ymax-dy;

	  Compute (& imloc, nb_iter, x_min, x_max, loc_ymin, loc_ymax);
      
	  MPI_Pack(&(imloc.pixels[0]), local_size+1, MPI_CHAR, buffer, 1<<20, &id, MPI_COMM_WORLD); // each non zero processor packs the result in buffer.
	}
      else if(rank==0)
	{
	  loc_ymax = y_max-dy*(size*jj);
	  loc_ymin = loc_ymax-dy;

	  Compute (& imloc, nb_iter, x_min, x_max, loc_ymin, loc_ymax); //proc 0 computes its local alloted image
	  l = (jj*size)*local_size;
	  for(int i=0; i<width; i++)
	    im.pixels[l+i] = imloc.pixels[i];
	}
    }

  if(rank!=0)
    MPI_Send(buffer, 1<<20, MPI_PACKED, 0, tag[rank], MPI_COMM_WORLD); // slave processors send the packed buffer to master processor
  else if(rank==0)
    {
      for(int r=1; r<size; r++)
	{
	  MPI_Recv(buffer, 1<<20, MPI_PACKED, r, tag[r], MPI_COMM_WORLD, &status); // buffer received at proc 0

	  id = 0;
	  for(int jj=0; jj<lsize; jj++)
	    {
	      l = (jj*size+r)*local_size;
	      MPI_Unpack(buffer, 1<<20, &id, &(im.pixels[l]), local_size+1, MPI_CHAR, MPI_COMM_WORLD); // proc 0 unpacks the local images and stores in the image location.

	    }
	} 
    }
    
    endtime = MPI_Wtime();
    
    if(rank==0)
        save (& im, path);
    
    printf("That took %f seconds in processor %d\n",endtime-starttime, rank);
    
    time_serial = 0.242479876449473;//from g5k serial run
    time_parallel = endtime - starttime;
    speed_up = time_serial/time_parallel;
    printf("The speed up obtained is: %f ",speed_up);
    
  MPI_Finalize();
}


