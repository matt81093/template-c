#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include "mpi.h"

double Timer();
float iteration(float **, float **, int, int, int);

float iterationDEBUG(float **, float **, int, int, int);

int main(int argc, char **argv)
{
	int i, rc, numranks, myrank;
	char * buff;
	int * len;
	
	int rows = atoi(argv[1]);
	int cols = atoi(argv[2]);
	float top = (float)atof(argv[3]);
	float left = (float)atof(argv[4]);
	float right = (float)atof(argv[5]);
	float bot = (float)atof(argv[6]);
	float e = (float)atof(argv[7]);
	
	float avg = top * (cols-2) + left * (rows-1) + right *(rows-1) + bot * cols;
	float val = rows*2 + cols*2 - 4;
	avg = avg/val;
	
    rc = MPI_Init(&argc,&argv);  // NULL,NULL ; rc == MPI_SUCCESS
    rc = MPI_Comm_size(MPI_COMM_WORLD,&numranks);
    rc = MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Status status;
	if(myrank==0)
	{
		double start = Timer();
		int counter = 0;
		float difference = 99.99;
		
		while(difference > e)
		{
			for(int i = 1; i < numranks; i++)
			{
				float temp;
				MPI_Recv(&temp, 1, MPI_FLOAT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status);
				if(temp < difference)
				{
					difference = temp;
				}
			}
				printf("%6d %.6f\n", counter, difference);
			MPI_Barrier(MPI_COMM_WORLD);
			if (counter && (ceil(log2(counter)) == log2(counter)))
				printf("%6d %.6f\n", counter, difference);
			
			if (difference <= e)
			{
				printf("%6d %.6f\n", counter, difference);
			}
			counter++;
				printf("not here0\n");
		}
		double end = Timer();
		end -= start;
		printf("TIME: %f\n", end);
		
		MPI_Abort(MPI_COMM_WORLD,0);//task failed successfully
	}
	else
	{
		int start,end;
		start = (rows-1)/(numranks-1);
		end = start*(myrank-1);
		start = start*(myrank-2)+1;
		if(myrank!=numranks)
		{
			end+=1;
			start-=1;
		}
		if((rows-1)%(numranks-1) && myrank == (numranks-1))
		{
			end+=(rows-1)%(numranks-1);
		}
		end = end-start;
		start = 0;
		
		float ** A = malloc(end*sizeof(float*));
		float ** B = malloc(end*sizeof(float*));
		
		for(int i = 0; i < end; i++)
		{
			A[i] = malloc(cols*sizeof(float));
			B[i] = malloc(cols*sizeof(float));
		}
		
		float * myDiff = malloc(4);
		if(myrank == 1)
		{
			for(int i = 0; i < cols; i++)
			{
				A[0][i] = top;
			}
			for(int i = 0; i < end; i++)
			{
				A[i][0] = left;
				A[i][cols-1] = right;
			}
			for(int i = 1; i < end; i++)
			{
				for(int j = 1; j < cols-1; j++)
				{
					A[i][j] = avg;
				}
			}
		}
		else if(myrank == numranks-1)
		{
			for(int i = 0; i < cols; i++)
			{
				A[end-1][i] = bot;
			}
			for(int i = 0; i < end; i++)
			{
				A[i][0] = left;
				A[i][cols-1] = right;
			}
			for(int i = 0; i < end-1; i++)
			{
				for(int j = 1; j < cols-1; j++)
				{
					A[i][j] = avg;
				}
			}
		}
		else
		{
			for(int i = 0; i < end; i++)
			{
				A[i][0] = left;
				A[i][cols-1] = right;
			}
			for(int i = 0; i < end; i++)
			{
				for(int j = 1; j < cols-1; j++)
				{
					A[i][j] = avg;
				}
			}
		}
		/*if(myrank == 1)
		{
			for(int i = 0; i < end; i++)
			{
				for(int j = 0; j < cols; j++)
				{
					printf("%.1f ",A[i][j]);
				}printf("\n");
			}
		}
			*/
		for (int i = 0; i < 9001; i++)
		{
			if(i%2 == 0)
			{
				*myDiff = iteration(A,B,start+1,end-1,cols);
			}
			else
			{
				*myDiff = iteration(B,A,start+1,end-1,cols);
			}
			
			printf("%d, %.6f\n",myrank,*myDiff);
			
			MPI_Send(myDiff,1,MPI_FLOAT,0,myrank,MPI_COMM_WORLD);
			if(numranks>2)
			{
				if(i%2 == 0)
				{
					if(myrank == 1)
					{
						MPI_Recv(&(B[end-1]), cols, MPI_FLOAT, myrank+1, 1, MPI_COMM_WORLD, &status);//receive from top
						MPI_Send(B[end-2], cols, MPI_FLOAT, myrank+1, 1, MPI_COMM_WORLD);//send to top
					}
					else if(myrank == numranks-1)
					{
						MPI_Send(B[0], cols, MPI_FLOAT, myrank-1, 1, MPI_COMM_WORLD);//send to bottom 
						MPI_Recv(&(B[1]), cols, MPI_FLOAT, myrank-1, 1, MPI_COMM_WORLD, &status);//receive from bottom
					}
					else
					{
						MPI_Recv(&(B[end-1]), cols, MPI_FLOAT, myrank+1, 1, MPI_COMM_WORLD, &status);//receive from top
						printf("%d,%.1f\n",myrank,(float)B[end-1][1]);
						MPI_Send(B[1], cols, MPI_FLOAT, myrank-1, 1, MPI_COMM_WORLD);//send to bottom
						MPI_Send(B[end-2], cols, MPI_FLOAT, myrank+1, 1, MPI_COMM_WORLD);//send to top
						MPI_Recv(&(B[0]), cols, MPI_FLOAT, myrank-1, 1, MPI_COMM_WORLD, &status);//receive from bottom
					}
				}
				else
				{
					if(myrank == 1)
					{
						MPI_Recv(&(A[end-1]), cols, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD, &status);//receive from top
						MPI_Send(A[end-2], cols, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);//send to top
					}
					else if(myrank == numranks-1)
					{
						MPI_Send(A[1], cols, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);//send to bottom 
						MPI_Recv(&(A[0]), cols, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD, &status);//receive from bottom
					}
					else
					{
						MPI_Recv(&(A[end-1]), cols, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD, &status);//receive from top
						MPI_Send(A[1], cols, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);//send to bottom
						MPI_Send(A[end-2], cols, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);//send to top
						MPI_Recv(&(A[0]), cols, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD, &status);//receive from bottom
					}
				}
			}
			
			
			MPI_Barrier(MPI_COMM_WORLD);
			printf("am i getting here%d\n",myrank);
		}
	}
	
	
	MPI_Finalize();
}


float iteration(float ** A, float ** B, int rowS, int rowE, int c)
{
	float differential = 0;
	for(int i = rowS; i < rowE; i++)
	{
		for(int j = 1; j < c-1; j++)
		{
			float sum = A[i-1][j] + A[i+1][j] + A[i][j-1] + A[i][j+1];
			float nVal = sum/4.0;
			float curDiff = fabs(nVal-A[i][j]);
			if (curDiff > differential)
			{
				differential = curDiff;
			}
			B[i][j] = nVal;
		}
	}
	return differential;
}
float iterationDEBUG(float ** A, float ** B, int rowS, int rowE, int c)
{
	float differential = 0;
	for(int i = rowS; i < rowE; i++)
	{
		for(int j = 1; j < c-1; j++)
		{
			float sum = A[i-1][j] + A[i+1][j] + A[i][j-1] + A[i][j+1];
			float nVal = sum/4.0;
			float curDiff = fabs(nVal-A[i][j]);
			if (curDiff > differential)
			{
				differential = curDiff;
			}
			B[i][j] = nVal;
		}
	}
	return differential;
}

double Timer()
{
    struct timeval tv;

    gettimeofday( &tv, ( struct timezone * ) 0 );
    return ( (double) (tv.tv_sec + (tv.tv_usec / 1000000.0)) );
}
