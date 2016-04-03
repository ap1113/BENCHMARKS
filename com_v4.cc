#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include "carbon_user.h"
#include <time.h>
#include <sys/timeb.h>

#define BILLION 1E9

typedef struct		//define a structure (class) with attributes
{
  int*      Q;
  // float**     W;
  // int**     W_index;
  int       tid;
  int       P;
  int       N;
  int       DEG;
  int 		iterations;
  //int**     sigma;
	  float*      mod_gain;
	  float*      total_mod_gain;
	  long*        comm;
	  float*      C;
  pthread_barrier_t* barrier;
} thread_arg_t;

// HERERERERERE //////
// int start = 64;
int next_source = -1;
thread_arg_t thread_arg[1024];
pthread_t thread_handle[1024];
 // int* avg;
 // double* delta;

pthread_barrier_t  barrier;       //barrier to synchronize threads
pthread_mutex_t id;               //single lock for 1 variable
pthread_mutex_t* locks;   			//locks for each vertex in the graph
int tid_g = -1;                     //to increment global thread IDs later

int* D;              
int* Q;
//int* avg;              
//double* delta;
float W[131072][32] = {0};  ;                          //graph weights 
long W_index[131072][32] = {0};                   //graph connections (adjacency list)
//int**     sigma; 
//int Total_tid[1024] = {0};        //triangles for each thread
int out_table[1024] = {0};
int in_table[1024] = {0};
int iterations = 0;
//int up = -1;
float* mod_gain;
float* total_mod_gain;
long* comm;
float* C;


long N=131072;                          //global, otherwise softbound complains
int DEG=32;
int P=0;
//long long Total = 0;              //Total triangles (a large number so long long)

volatile int passed = 0; // Number of barriers, passed by all threads.
int bar=0;
int PMAX=0;
void* barrier_wait();		//Initialize wait function
void* threadMain(void* arg);                //function initialization

void* barrier_wait()
{
    int passed_old = passed; // Should be evaluated before incrementing *bar*!

    if(__sync_fetch_and_add(&bar,1) == (PMAX - 1)) // Not happening ...
    {   
        // The last thread, faced barrier.
        bar = 0;
        // *bar* should be reseted strictly before updating of barriers counter.
        __sync_synchronize(); 
        passed++; // Mark barrier as passed.
    }   
    else
    {   
        // Not the last thread. Wait others.
        while(passed == passed_old) {}; 
        // Need to synchronize cache with other threads, passed barrier.
        __sync_synchronize();
    } 
	return 0;
}



void* threadMain(void* arg)                 //main parallel function
{
  //each thread gets it IDs
  // pthread_mutex_lock(&id);
	// tid_g++;
  // int tid = tid_g;
   int tid=0;
  tid = __sync_fetch_and_add(&tid_g,1);
  tid++;
  printf("\n TID: %d",tid);
  // pthread_mutex_unlock(&id);

	// pthread_barrier_wait(&barrier);
	// barrier_wait();

   // //divide graph among threads
	long start = tid    *  (N) / (P);
	long stop = (tid+1) *  (N) / (P);
	// // long node = 0;
	
	// //OTHER STUFF
	  float modularity = 0;
	  long index = 0;
	  int index_id = 0;
	  float sum_tot = 0;
	  float sum_in = 0;
	  int local_count = 0;

	  // //float total_edges = N*DEG;
	  float mod_gain_temp = 0;
	  float mod_gain_temp_temp = 0;
	  // //float inv_total_edges = 2/total_edges;
   
   // // pthread_barrier_wait(&barrier);
   printf("\n %d b0",tid);
	// barrier_wait(); //this 
	 

   // //add edges to D array
	for(long uu=start;uu<stop;uu++)
	  {
		  comm[uu] = uu;
	  }
	// //up = 0;
	// // pthread_barrier_wait(&barrier);
	barrier_wait();
   

	while(local_count<iterations)
	{
		for(long uu=start;uu<stop;uu++)
	   {
			for(int i = 0; i < DEG; i++)
			{
			    float total_edges = (N-1)*DEG;
				float inv_total_edges = 2/total_edges;
				int tempo = (inv_total_edges)*(inv_total_edges);
				tempo = tempo*2;
				float subtr = sum_tot*W[uu][i];
				subtr = subtr*tempo;
				subtr = W[uu][i] - subtr;
				mod_gain_temp_temp = (inv_total_edges)*(subtr);
				//total_mod_gain[uu] = total_mod_gain[uu] + mod_gain[uu];

				if(mod_gain_temp_temp>mod_gain_temp) 
				{
				  mod_gain_temp = mod_gain_temp_temp;
				  index = W_index[uu][i];
				  index_id = i;
				}
			}
		   mod_gain[uu] = mod_gain_temp;
		  comm[uu] = index;   //cvk


		  // //update individual  sums
		  // //pthread_mutex_lock(&lock);
		  sum_tot = sum_tot + W[uu][DEG-1];
		  sum_in = sum_in + W[uu][DEG-1];	   
	   
	    }
		
		// // pthread_barrier_wait(&barrier);
		barrier_wait();
		for(long uu=start;uu<stop;uu++)
		{
		  for(int i=0;i<DEG;i++)
		  {
			long neighbor = W_index[uu][i];  //this 
			// //pthread_mutex_lock(&locks[neighbor]);
			  // //W_index[uu][i] = comm[neighbor];
			  W[uu][i] = comm[uu] - comm[neighbor];  //this
			// //pthread_mutex_unlock(&locks[neighbor]);
		  }
		}
	  
	    // // pthread_barrier_wait(&barrier);
	barrier_wait(); //this
	    local_count++;
			
	}
	
	for(long i=stop;i<N;i++)
      {
        // //for(int j=1;j<P;j++)
        // //{
          comm[i] = comm[((1) *  (N) / (P))-1];
        // //}
      }
  
	// // pthread_barrier_wait(&barrier);
    printf("\n %d bf",tid);
	barrier_wait();
    return 0;
}

int main(int argc, char *argv[])
{

   int N1 = 16384; //atoi(argv[3]);  //1024;         //input arguments
   int DEG1 = 32; //atoi(argv[4]); //16;
   int ITER1 = 2; //atoi(argv[2]);
   int P1 = atoi(argv[1]); //8;

     // N=N1;
	 // DEG=DEG1;
	 P=P1;
	 PMAX = P;		//for the new barrier wait issue
	 iterations = ITER1;
	 if (DEG > N)
	{
		fprintf(stderr, "DEG cannot be greater than N (0 P N DEG)\n");
		exit(EXIT_FAILURE);
	}

   // D = (int*) malloc((N)*sizeof(int));            //memory allocations
   int D1[N];
   if(!D) {D = D1;}
   // Q = (int*) malloc((N)*sizeof(int));
   int Q1[N];
   if(!Q) {Q = Q1;}
   // comm = (int*) malloc((N)*sizeof(int));            
   long comm1[N];
   if(!comm) {comm = comm1;}
   // C = (float*) malloc((N)*sizeof(float));
   float C1[N];
   if(!C) {C = C1;}
   // mod_gain = (float*) malloc((N)*sizeof(float));
   float mod_gain1[N];
   if(!mod_gain) {mod_gain = mod_gain1;}
   // total_mod_gain = (float*) malloc((N)*sizeof(float));
   float total_mod_gain1[N];
   if(!total_mod_gain) {total_mod_gain = total_mod_gain1;}
   // W = (float**)malloc((N)*sizeof(float*));
   // W_index = (int**)malloc((N)*sizeof(int*));
   //sigma = (int**) malloc(P*sizeof(int*));
   locks = new pthread_mutex_t[N];

   // for (int i = 0; i < N; i++)
   // {
	  // W[i] = (float*) malloc((DEG)*sizeof(float));
	  // W_index[i] = (int*) malloc((DEG)*sizeof(int));
	  // //sigma[i] = (int*) malloc((DEG)*sizeof(int));
   // }

pthread_mutex_init(&id, NULL);                         //barrier and lock inits
pthread_barrier_init(&barrier, NULL, P);
    for(long i=0; i<N; i++)
		pthread_mutex_init(&locks[i], NULL);

//fill D and Q
for(long i = 0; i < N; i++)
{
    D[i] = 0;
    Q[i] = 1;
	// delta[i]=0;
	// avg[i]=0;
}

//Now fill the graph with synthetic stuff (same as init_weights in orig)

          int range = DEG + (N - DEG)/16;
  
          // Initialize to -1
          for(long i = 0; i < N; i++)
                  for(int j = 0; j < DEG; j++)
                          W_index[i][j]= -1;

          // Populate Index Array
          for(long i = 0; i < N; i++)
          {
                  long last = 0;
                  int min = 0;
                  int max = DEG;
                  for(int j = 0; j < DEG; j++)
                  {
                          if(W_index[i][j] == -1)
                          {        
								long neighbor = i+j;																	
                                  if(neighbor > last)
                                  {
                                          W_index[i][j] = neighbor;
                                          last = W_index[i][j];
                                  }
                                  else
                                  {
                                          if(last < (N-1))
                                          {
                                                  W_index[i][j] = (last + 1);
                                                  last = W_index[i][j];
                                          }
                                  }
                          }
                          else
                          {
                                  last = W_index[i][j];
                          }
                          if(W_index[i][j]>=N)
                          {
                                  W_index[i][j] = N-1;
                          }
                  }
          }
  
          // Populate Cost Array
          for(long i = 0; i < N; i++)
          {
                  for(int j = 0; j < DEG; j++)
                  {
                          double v = drand48();
  
                          if(W_index[i][j] == i)
                                  W[i][j] = 0;
  
                          else
                                  W[i][j] = (int) (v*100) + 1;	
                  }
          }
		  
//Graph filling finished


   // Enable performance and energy models
   CarbonEnableModels();
   //printf("turned on\n");
   struct timespec requestStart, requestEnd;
   clock_gettime(CLOCK_REALTIME, &requestStart);  //to measure time
	
   pthread_t thread_handles[P];              //spawning threads
   //printf("threads array generated\n");
   for (int i = 1; i < P; i++)
   {
       // printf("creating tid: %d\n", i);
     pthread_create(thread_handle+i, NULL, threadMain, (void*)&thread_arg[i]); 
	 // printf("created tid: %d\n", i);
   }
	 threadMain((void*) &thread_arg[0]);

   /*for (int i = 1; i < P; i++)
   {
      pthread_join(thread_handles[i], NULL);
   }*/

   
   
   	clock_gettime(CLOCK_REALTIME, &requestEnd);
	  double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;
	  printf( "\nTime:%lf\n", accum );   //prints time
// Disable performance and energy models
   CarbonDisableModels();

   printf("Community: Done.\n");
   
   // Free memory

   //CarbonStopSim();
   return 0;
}
