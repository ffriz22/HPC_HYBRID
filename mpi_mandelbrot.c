//
//  mandelbrot.c
//  
//
//  The Mandelbrot calculation is to iterate the equation
//  z = z*z + c, where z and c are complex numbers, z is initially
//  zero, and c is the coordinate of the point being tested. If
//  the magnitude of z remains less than 2 for ever, then the point
//  c is in the Mandelbrot set. In this code We write out the number of iterations
//  before the magnitude of z exceeds 2, or UCHAR_MAX, whichever is
//  smaller.//
//
//

/*
    New on this version (v9):
	- Print num_thread in subtask results
*/


#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#include <omp.h>
#include <mpi.h>

int omp_get_thread_num();

const debug = 0;		// Debug messages must be printed

/* Structure definition for complex numbers  */

typedef struct {
    double real, imag;
} complex;

typedef unsigned char colors[3];

/*
void color(int red, int green, int blue)
{
    fputc((char)red, stdout);
    fputc((char)green, stdout);
    fputc((char)blue, stdout);
}
*/

//Editado version antigua es: void color_print(int x, int y, int red, int green, int blue)
void color_print(int x, int y, char red)
{
    char green, blue;

    if (red != 0) {
	green = red;
	blue = (char)255;
    } else {
	green = red;
	blue = red;
    }
    fputc((char)red, stdout);
    fputc((char)green, stdout);
    fputc((char)blue, stdout);

    if (debug) fprintf(stderr, "(%d,%d)=%d ", x, y, (int)red, omp_get_thread_num());
}


int main(int argc, char **argv)
{
	//each iteration, it calculates: newz = oldz*oldz + p, where p is the current pixel, and oldz stars at the origin
    
    int x, y;			// Pixel position   
    complex p;			// pixel p	
    complex newZ, oldZ;		// new and old z
    double zoom = 1, moveX = -0.5, moveY = 0; // zoom and position
    int ompThreads = 4;		// Number of OMP threads in every MPI task 
    double begin, end;
    double time_spent, time_calc;
    int yini, ycount;		// Initial y position and number of rows
    char colour;
    int brightness;

    /* MPI variables	*/
    int numTasks, namelen, rank, dest = 1, tag = 111, source = 0, size, rc, i;
    int sendcount, reccount;
//    int *sendbuf;
//    int *recbuf;
    char hostname[256];
    MPI_Request request1, request2;     
    MPI_Status status;
    
    /* Process command-line arguments  */
    int maxIterations; 		// Number of iterations in every pixel
    int w, h;			// Witdh and Height of the image
    int numThreads;		// Total threads (OMP + MPI)
    
    if (argc < 4) {
    	printf( "usage:  %s maxiter width height numthreads\n");
	return 0;
    }
    maxIterations = atoi(argv[1]);
    w = atoi(argv[2]);
    h = atoi(argv[3]);
    numThreads = atoi(argv[4]);
    //numTasks = (int) (numThreads / 4);

    int *sposbuf= (int *)malloc(sizeof(int)*2);	// Buffer to send the position to calculate
    int *rposbuf = (int *)malloc(sizeof(int)*2); // Buffer to receive the position to calculate
	   
    /* MPI configuration	*/
    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);   // get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // get current process id
    MPI_Get_processor_name(hostname, &namelen); // get CPU name

    if (rank == 0) {
    	fprintf(stderr,"Iteration : %d\n", maxIterations);
    	fprintf(stderr,"Width     : %d\n", w);    
    	fprintf(stderr,"Height    : %d\n", h);
    	fprintf(stderr,"NumThreads: %d\n", numThreads);
    	fprintf(stderr,"NumTasks  : %d\n", numTasks);
 
    	printf("P6\n# CREATOR: Franco Friz & David Sanchez \n");
    	printf("%d %d\n255\n",w,h);
    }

    fprintf(stderr, "Hostname: %s, numTasks: %d, rank: %d\n", hostname, numTasks, rank);

    //omp_set_num_threads(ompThreads);
    
    begin = MPI_Wtime();

    // Send y positions
    //
    sendcount = (int) h/numTasks;
    if (rank == 0)	// Master task 
    {
	fprintf(stderr, "Task 0 sends to task %d: y= %d and count= %d \n", rank, 0, sendcount);

    	for (i = 1; i < numTasks; i++) { // Send message to all the task with
	    y = i * sendcount;	    // Initial y position
	    sposbuf[0] = y;
	    if (y + sendcount < h)  // Number of rows to calculate
	        sposbuf[1] = sendcount;
	    else
		sposbuf[1] = h - y;
	    fprintf(stderr, "Task 0 sends to task %d: y= %d and count= %d \n", i, sposbuf[0], sposbuf[1]);
            MPI_Isend(sposbuf, 2, MPI_INT, i, tag, MPI_COMM_WORLD, &request1);
	    MPI_Wait(&request1, &status);	    
 	}
	
	fprintf(stderr, "Task %d has received:: y= %d and count= %d \n", rank, 0, sendcount);
	yini = 0;
	ycount = sendcount;


    	end = MPI_Wtime();
    	time_spent = (double)(end - begin);
	fprintf(stderr, "\nSending time: %.2lf seconds.\n", time_spent);


    } 
    else		// Slave tasks
    {
	MPI_Irecv(rposbuf, 2, MPI_INT, 0, tag, MPI_COMM_WORLD, &request2);
	MPI_Wait(&request2, &status);

	yini = rposbuf[0];
	ycount = rposbuf[1];
        fprintf(stderr, "Task %d has received: y= %d and count= %d \n", rank, yini, ycount);	
    }
    free(sposbuf);
    free(rposbuf);

//	Calculate results

//	NEW
    char *sresulbuf = (char *)malloc(sizeof(char) * w * ycount);  // Buffer to return the value of points in a row
    char *rresulbuf = (char *)malloc(sizeof(char)* w * ycount * numTasks);  // Buffer to return the value of points in a row
 
    if (debug) fprintf(stderr, "Results: \n"); 
    
    omp_set_num_threads(numThreads);
    for (y = yini; y < yini + ycount; y++) {	//  y positions calculated by this task 
//        y = yini + i;

	#pragma omp parallel for schedule(dynamic) private (x, p, newZ, oldZ)
        for(x = 0; x < w; x++)		// x positions
        {
	    int numthread = omp_get_thread_num();

            //calculate the initial real and imaginary part of z, based on the pixel location and zoom and position values
            p.real = 1.5 * (x - w / 2) / (0.5 * zoom * w) + moveX;	
            p.imag = (y - h / 2) / (0.5 * zoom * h) + moveY;	
            
            newZ.real = newZ.imag = oldZ.real = oldZ.imag = 0; //these should start at 0,0
            
            //"i" will represent the number of iterations
            int i;
            //start the iteration process
            for(i = 0; i < maxIterations; i++)
            {
                //remember value of previous iteration
                oldZ = newZ;
                
                //the actual iteration, the real and imaginary part are calculated
                newZ.real = oldZ.real * oldZ.real - oldZ.imag * oldZ.imag + p.real;
                newZ.imag = 2 * oldZ.real * oldZ.imag + p.imag;
                
                if((newZ.real * newZ.real + newZ.imag * newZ.imag) > 4) break;
            }
            
            if(i == maxIterations) 
		    {
				// 	black color
		     brightness = 0;
		     colour= (char) brightness;	//	NEW
		     sresulbuf[(y - yini)*w + x] = colour;  // NEW	guardar resultado brightness (en negro)	     
//		     sresulbuf[y - yini][x] = brightness;
	    	}
	    else
	        {
	            double z = sqrt(newZ.real * newZ.real + newZ.imag * newZ.imag);
	            brightness = 256 * log2(1.75 + i - log2(log2(z))) / log2((double)maxIterations);
		
	            colour= (char) brightness;	// NEW
	     	    sresulbuf[(y - yini)*w + x] = colour; // NEW	guardar resultado brightness (real)	     
	        }
	        if (debug) fprintf(stderr, "(%d,%d):%d=%d ", x, y, numthread, brightness);
        }
	if (debug) fprintf(stderr, "\n");		
    }
    if (debug) fprintf(stderr, "\n");

    if (rank == 0) {
    	end = MPI_Wtime();
    	time_spent = (double)(end - begin);
	fprintf(stderr, "\nCalculation time: %.2lf seconds.\n", time_spent);
    }

//	Return results

    MPI_Gather(sresulbuf, w * ycount, MPI_CHAR, rresulbuf, w * ycount, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
    	end = MPI_Wtime();
    	time_spent = (double)(end - begin);
	fprintf(stderr, "\nReturn time: %.2lf seconds.\n", time_spent);
    }

// Guardar en un fichero el mpi: NEW

    colors *result= (colors *)malloc(sizeof(colors)*w*h);

    if (rank == 0) {
	if (debug) fprintf(stderr, "\nFINAL RESULT: \n", i);

	for (i = 0; i < w * h; i++) {
	   y = (int) i / w;
	   x = i - y * w;
	   //color_print(x, y, rresulbuf[i]);
        if (rresulbuf[i] != 0) {
            result[i][0] = rresulbuf[i];
            result[i][1] = rresulbuf[i];
            result[i][2] = (char)255;
        } else {
            result[i][0] = rresulbuf[i];
            result[i][1] = rresulbuf[i];
            result[i][2] = rresulbuf[i];
        }
	}

    fwrite(result,sizeof(colors),w*h,stdout);	

/*
	    for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
                // Posicion x,y
                // Contenido rresulbuf[y*w+x] -> red
                
                    color_print(x, y, rresulbuf[y*w+x]);
		// rresulbuf[y][x];
		}
	        if (debug) fprintf(stderr, "\n");
	    }
	    fprintf(stderr, "\n");
*/
    }
   
    free(sresulbuf);
    free(rresulbuf);

    if (rank == 0) {
    	end = MPI_Wtime();
    	time_spent = (double)(end - begin);
	fprintf(stderr, "\nExecution time: %.2lf seconds.\n", time_spent);
    }
    fprintf(stderr," \n");    
    MPI_Finalize ();

    return 0;
}
