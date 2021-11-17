#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#define TIME 100
#define BINS 20
#define ROOT 0
#define R 15
void prop(int*, double*);

/* ********** Stochastic Malaria Simulation* ********* */
int main(int argc, char *argv[]) {

    if (argc != 2) {
    printf("Error, expected input: mpirun -np p ./malaria N \n");
    printf("[p = processes, N = experiments] \n");
    return -1; }

    /* ---- Initialize MPI ---- */
    int size, rank;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* ---- Each process executes SSA n times, N/p = n ---- */
    int N = atoi(argv[1]);
    if ((N % size) != 0) {
        printf("Error, N not divisible by p. \n");
        return -1; }
    int n = N/size; // Runs per process

    /* ---- Time ---- */
    double t = 0;
    double t_final = TIME;
    double start, end;

    /* ---- Parameters ---- */
    double a0, u1, u2, tau;
    double max, min, interval;
    int local_max, local_min, r;

    /* ---- Vectors, matrices ---- */
    double w[15];      // Propensities
    int x[7];          // State vector
    int result[n][7];  // Result matrix
    int* susceptible_humans = malloc(n*sizeof(int));
    int bins[BINS+1];
    int sparse_P[15][7] = {
        { 1, 0, 0, 0, 0, 0, 0},
        {-1, 0, 0, 0, 0, 0, 0},
        {-1, 0, 1, 0, 0, 0, 0},
        {0,  1, 0, 0, 0, 0, 0},
        {0, -1, 0, 0, 0, 0, 0},
        {0, -1, 0, 1, 0, 0, 0},
        {0, 0, -1, 0, 0, 0, 0},
        {0, 0, -1, 0, 1, 0, 0},
        {0, 0, 0, -1, 0, 0, 0},
        {0, 0, 0, -1, 0, 1, 0},
        {0, 0, 0, 0, -1, 0, 0},
        {0, 0, 0, 0, -1, 0, 1},
        {0, 0, 0, 0, 0, -1, 0},
        {1, 0, 0, 0, 0, 0, -1},
        {0, 0, 0, 0, 0, 0, -1}
    };

    srand(time(NULL) + rank); // <-- Important: random seed for each process!
    start = MPI_Wtime();

    /* --- Each process does n number SSA runs ----*/
    for (int i = 0; i < n; i++) {
        t = 0;
        /* Inital State */
        x[0] = 900;
        x[1] = 900;
        x[2] = 30;
        x[3] = 330;
        x[4] = 50;
        x[5] = 270;
        x[6] = 20;
        /* SSA loop */
    	while(t < t_final) {
    	    /* Propensities */
    	    prop(x, w);

    	    /* Compute a0 */
    	    a0 = 0;
    	    for(int j = 0; j < R; j++) {
                a0 += w[j];
    	    }

    	    /* Generate two uniform random numbers u1, u2 */
    	    u1 = rand() / (float) RAND_MAX;
    	    u2 = rand() / (float) RAND_MAX;

    	    /* Set tau */
    	    tau = -log(u1)/a0;

            /* Find r */
            double sum = 0;
            for(int i = 0; i < R; i++) {
                sum += w[i];
                if (sum >= (a0*u2)) { r = i; break; }
            }
            /* Update state vector x */
            x[0] += sparse_P[r][0];
            x[1] += sparse_P[r][1];
            x[2] += sparse_P[r][2];
            x[3] += sparse_P[r][3];
            x[4] += sparse_P[r][4];
            x[5] += sparse_P[r][5];
            x[6] += sparse_P[r][6];
    	    t = t + tau; // Update current time
    	}
        /* Add x vector to result matrix */
        result[i][0] = x[0];
        result[i][1] = x[1];
        result[i][2] = x[2];
        result[i][3] = x[3];
        result[i][4] = x[4];
        result[i][5] = x[5];
        result[i][6] = x[6];

    }
    /* ---- Store first part of x (susceptible humans)  ---- */
    for (int i = 0; i < n; i++) {
        susceptible_humans[i] = result[i][0];
    }

    /* ---- Find max and min values locally  ---- */
    int c;
    int max_location = 0;
    int min_location = 0;
    for (c = 1; c < n; c++) {
        if (susceptible_humans[c] > susceptible_humans[max_location]) { // MAX
            max_location = c; }
        if (susceptible_humans[c] < susceptible_humans[min_location]) { // MIN
            min_location = c; }
    }
    local_max = susceptible_humans[max_location];
    local_min = susceptible_humans[min_location];


    /* ---- Collect global min and max values and distribute to all procs  ---- */
    MPI_Allreduce(&local_max, &(bins[BINS+1]), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min, &(bins[0]), 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    /* ---- Compute bin sizes  ---- */
    max = (int) bins[BINS+1];
    min = (int) bins[0];
    interval = (max - min)/BINS;

    /* ---- Add remaining bins to bin array ---- */
    for (int i = 0; i < BINS+1; i++) {
	       bins[i] = min + (i * interval);
    }

    /* ---- Count number of elements in each bin group locally ---- */
    int local_freq[BINS] = {0};
    int count = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < BINS; j++) {
            if (susceptible_humans[i] < bins[j+1] && susceptible_humans[i] >= bins[j]) {
                local_freq[j]++; count++; break; }
        }
    } local_freq[BINS-1] = local_freq[BINS-1] + (n-count);

    /* ---- Sum the elements in the bin groups globally ---- */
    int* frequencies;
    if (rank == ROOT) {
	       frequencies = malloc(sizeof(int)*BINS);
    }
    MPI_Reduce(local_freq, frequencies, BINS, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);

    end = MPI_Wtime();

    /* ---- Find mean and standard deviation  ---- */
    /*      Uncomment the code to get mean and SD  */
    // int local_sum = 0;
    // int global_sum;
    // for (int i = 0; i < n; i++) {
    //     local_sum = local_sum + susceptible_humans[i]; }
    // MPI_Allreduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // float mean = global_sum / N;
    // float local_sq_diff = 0;
    // for (int i = 0; i < n; i++) {
    //     local_sq_diff += (susceptible_humans[i] - mean) * (susceptible_humans[i] - mean);
    // }
    // float global_sq_diff;
    // MPI_Reduce(&local_sq_diff, &global_sq_diff, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);


    /* ----  Print results ---- */
    if(rank == 0){
        printf("Time: %0.4f\n", end-start); // Time
        // float stddev = sqrt(global_sq_diff/N);
        // printf("Mean: %0.3f, Standard deviation = %f \n", mean, stddev); // Mean, std

        /* Print output */
        printf("Frequency: \n");
        int s = 0;
        for (int i = 0; i < BINS; i++) {
            printf("%d ", frequencies[i]);
        }
        printf("\nIntervals: \n");
        for (int i = 0; i < BINS+1; i++) {
            printf("%d ", bins[i]); }
        printf("\n");
	free(frequencies);
    }
    /* ---- Finalize MPI ---- */
    free(susceptible_humans);
    MPI_Finalize();
    return 0;
}
