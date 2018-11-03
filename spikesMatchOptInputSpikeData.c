#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
typedef int bool;
#define true 1
#define false 0



// membrane resistance
double R = 100000000;

// for generating random numbers 
int scaler = 1000;

int numIter = 200;


// specify input current distribution

double inputMean =   .00000000025;

double inputVar  =   .0000000001;

double refScaling = .23;


//network size; ie. number of neurons

int size = 500;

// voltage threshold in Volts

double voltageThreshold = .03;

// Reciprocal of membrane time constant in Hz

double Tau = 30.0;


//simulation step size in seconds

double step = 0.003;

//number of time steps

//int numSteps = 10000;

int numSteps = 10000;

// for generating random numbers


//double factor = .000005;

double factor =     .000001;

//double factor = 0;


int maxLoopCount = 100;


// for reading in a weight matrix from a 'weights.csv' file in the directory if provided
double weightScale = 3;


// maximum size of an individual weight modification that can be made during optimization


//double epsilon =.0000008;

double epsilon =  .00000001;

double weightMean = .0004;
//double weightVar = 0;
double weightVar  = .0004;

//double weightMean = 0;
//double weightVar = 0;

//double rescale = .000000001;

//double rescale =   .0000000001;
//double rescale   =   .0000000025;

//double rescale = .0000000001;

double rescale =   .000000001;
int maxSepDistDecrementer = 15;

// simulate the network with given input stimulating currents. Outputs are spike trains.

void simulate2(double ** outputs, double ** network, double ** inputs, double ** vstatemat, double * synCurrentRatio, double * synPercent);
void simulate1(double ** outputs, double ** network, double ** inputs, double ** vstatemat);


// perform basic matrix-vector multiplication
void matrixMult(double ** A, double * B);


// optTrain: a spike train from the optimized network; reftrain: the corresponding spike train in the reference network
// maxSpikesAndPairs is a double pointer: the first location contains a simple 1 if the opt spike train has more spikes and a 0 if the reference network spiked more 
// using the information from maxSpikesandPairs, we can compute the error and which weights to reduce/increase
double pairSpikes(double * optTrain, double * refTrain, double * maxSpikesAndPairs, int * a, int * maxSepDist, int ** d, int **M);


// to load network in 
void loadNetwork(double ** network, int size);

// to load a set of spikes in
void loadSpikes(double ** inputSpikes);



int main()
{
	
	
        time_t t;

        srand(time(NULL));

        int i, j, k;
        double **inputs = (double **) malloc(size*sizeof(double *));
        for(i=0; i<size; i++)
        {
                inputs[i] = (double *) malloc(numSteps*sizeof(double));

        }
        double **outputs = (double **) malloc(size*sizeof(double *));
        for(i=0; i<size; i++)
        {
                outputs[i] = (double *) malloc(numSteps*sizeof(double));
        }

        double **optOutputs = (double **) malloc(size*sizeof(double *));
        for(i=0; i<size; i++)
        {
                optOutputs[i] = (double *) malloc(numSteps*sizeof(double));
        }


        double **outputsCopy = (double **) malloc(size*sizeof(double *));
        for(i=0; i<size; i++)
        {
                outputsCopy[i] = (double *) malloc(numSteps*sizeof(double));
        }

        double **optOutputsCopy = (double **) malloc(size*sizeof(double *));
        for(i=0; i<size; i++)
        {
                optOutputsCopy[i] = (double *) malloc(numSteps*sizeof(double));
        }

	
	double **refNet = (double **) malloc(size*sizeof(double *));
        for(i=0; i<size; i++)
        {
               refNet[i] = (double *) malloc(size*sizeof(double));
        }


        double ** optNet = (double **) malloc(size*sizeof(double *));
        for(i=0; i<size; i++)
        {
               optNet[i] = (double *) malloc(size*sizeof(double));
        }
	
        
	double ** optNetControl = (double **) malloc(size*sizeof(double *));
	for(i=0; i<size; i++)
	{
		optNetControl[i] = (double *) malloc(size*sizeof(double));
	}

        double ** spikeDomAndPairs = (double **) malloc(size*sizeof(double *));

        for(i=0; i<size; i++)
        {
                spikeDomAndPairs[i] = (double *) malloc(2*(numSteps)*sizeof(double));
                
        }

        double ** optNetCopy = (double **) malloc(size*sizeof(double*));
        for(i=0; i<size; i++)
        {
                optNetCopy[i] = (double *) malloc(size*sizeof(double));
                for(j=0; j<size;j++)
                {
                        optNetCopy[i][j] = optNet[i][j];
                }
        }

        int ** d = (int **) malloc(numSteps*sizeof(int *));
        for(i=0; i<numSteps; i++)
        {
                d[i] = (int *) malloc(numSteps*sizeof(int));
        }

        int ** M = (int **) malloc(numSteps*sizeof(int *));
        for(i=0; i<numSteps; i++)
        {
                M[i] = (int *) malloc(numSteps*sizeof(int));
        }


	double ** refstatemat = (double **) malloc(size*sizeof(double *));

        for(i=0; i<size; i++)
        {
                refstatemat[i] = (double *) malloc(numSteps*sizeof(double));
        }


	
	// initialize reference and naive networks. Here the optimized network is drawn from a distribution with upper bound 1.2 times that of ref's. In addition, 40% of reference network's initial weights set to 0.


	// in this file: same inhibitory neurons for opt and ref
	
	
//	int numInhibitory = (int) (size * .2);

//	int numInhibitory = (int) (.2*size);

	int numInhibitory = 0;

	int * inhibIndices = (int *) malloc(numInhibitory * sizeof(int));
	//set all inhib indices to 0


	for(i=0; i<numInhibitory; i++)
	{
		int index = (rand() % size);
		
		int inList = 1;
		while(inList == 1)
		{
			inList = 0;
			
			for(j=0; j<numInhibitory; j++)
			{
				if(inhibIndices[j] == index)
				{
					inList = 1;
					index = rand() % size;
				}
			}
		}
		inhibIndices[i] = index;

		printf("inhibIndices[%d] = %d\n", i, inhibIndices[i]); 
	}


	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			optNet[i][j] = (rand() % scaler) * factor;
			refNet[i][j] = (rand() % scaler) * factor;
			
			
			

			for(k = 0; k<numInhibitory; k++)
			{
				if(j == inhibIndices[k])
				{
					optNet[i][j] = -optNet[i][j];
					refNet[i][j] = -refNet[i][j];
				}
			}	
			

			optNetControl[i][j] = optNet[i][j];
	//		printf("optNet[%d][%d] = %lf\n", j, i, optNet[j][i]);
		}
	}



        for(i=0; i<size; i++)
        {
                for(j=0; j<size; j++)
                {
                        double U1 = 1.0*(rand() % 1000) / 1000.0;
                        double U2 = 1.0 * (rand() % 1000)/1000.0;

                        while(U1 == 0)
                        {
                                U1 = 1.0*(rand() % 1000) / 1000.0;

                        }
                        double Z1 = sqrt(-2*log(U1))*cos(2*3.14159 * U2);

                        double val = weightVar * Z1 + weightMean;
                        optNet[i][j] = val;
			optNetControl[i][j] = optNet[i][j];
                }
        }

	









        FILE * naiveweights;
        naiveweights = fopen("naiveweights.txt","w+");
        for(i=0;i<size;i++)
        {
                for(j=0;j<size;j++)
                {
                        fprintf(naiveweights,"%lf,",optNet[i][j]);
                }
                fprintf(naiveweights,"\n");
        }
        fclose(naiveweights);
	







	// generate Gaussian input currents
	
	for(i=0; i<size; i++)
	{
		for(j=0; j<numSteps; j++)
		{
			double U1 = 1.0*(rand() % 1000) / 1000.0;
			double U2 = 1.0 * (rand() % 1000)/1000.0;
			
			while(U1 == 0)
			{
				U1 = 1.0*(rand() % 1000) / 1000.0;

			}
			double Z1 = sqrt(-2*log(U1))*cos(2*3.14159 * U2);
			
			double val = inputVar * Z1 + inputMean;
			inputs[i][j] = val;
		}
	}


	double * refSynCurrentRatio = (double *) malloc(size*sizeof(double));
	double * refSynPercent = (double *) malloc(size*sizeof(double));
	double * synCurrentRatio = (double *) malloc(size*sizeof(double));
	double * synPercent = (double *) malloc(size*sizeof(double));

	
	for(i=0; i<size; i++)
	{
		refSynCurrentRatio[i] = 0;
		refSynPercent[i] = 0;
		synCurrentRatio[i] = 0;
		synPercent[i] = 0;
	}
	

	// generate outputs by simulating ref
//	simulate2(outputs, refNet, inputs, refstatemat, refSynCurrentRatio, refSynPercent);

	

	loadSpikes(outputs);


	int numSpikesOutputs = 0;

	for(i=0; i<size; i++)
	{
		for(j=0; j<numSteps; j++)
		{
			numSpikesOutputs += outputs[i][j];
		}
	}
	double spikeRate = numSpikesOutputs / (size*numSteps*step);	
	printf("average spike rate in ref: %lf\n Hz", spikeRate);


	FILE * fer = fopen("refSynCurrentRatio.txt","w+");
	
	FILE * fer2 = fopen("refSynPercent.txt","w+");



	double ** rememberOutputs = (double **) malloc(size * sizeof(double *));
        double ** rememberOptOutputs = (double **) malloc(size * sizeof(double *));

	for(i=0; i<size; i++)
	{
		rememberOutputs[i] = (double *) malloc(numSteps * sizeof(double));
                rememberOptOutputs[i] = (double *) malloc(numSteps * sizeof(double));

		for(j=0; j<numSteps; j++)
		{
	                rememberOutputs[i][j] = outputs[i][j];
	                outputsCopy[i][j] = outputs[i][j];

		}
	}
	
	printf("simulated refNet\n");

double error = 500;
double threshold = 100;
int loopcount = 0;
int r,startError,endError;


// z is the number of timesteps prior to a spike event that spikes in the presynaptic pool will be considered as causally responsible for that spike event
int * z = (int *)malloc(sizeof(int)); 
*z = 15;


// maxSepDist is the greatest number of timesteps between two spikes that will allow them to be paired.
int * maxSepDist = (int *) malloc(sizeof(int));
*maxSepDist = 15;

int optUnpairedSpikeNum = 0;
int refUnpairedSpikeNum = 0;


// simulate optNet

simulate2(optOutputs, optNet, inputs, refstatemat, synCurrentRatio, synPercent);

int countspikes = 0;
for(i=0; i<size; i++)
{
	for(j=0; j<numSteps; j++)
	{
		if(optOutputs[i][j]==1)
		{
			countspikes++;
		}
	}
}

printf("countspikes: %d\n", countspikes);
printf("spikerate: %lf\n", countspikes/(size*numSteps*step));

FILE * initialoutputs = fopen("initialoutputs.txt","w+");

for(i=0; i<size; i++)
{
	for(j=0; j<numSteps; j++)
	{
		fprintf(initialoutputs, "%lf,",optOutputs[i][j]);

	}
	fprintf(initialoutputs, "\n");
}

int originalSpikeNum = 0;
for(i=0; i<size; i++)
{
	for(j=0; j<numSteps; j++)
	{
		originalSpikeNum += optOutputs[i][j];
		rememberOptOutputs[i][j] = optOutputs[i][j];
	}
}



double initialError = 0;



int * a = (int *) malloc(sizeof(int));
*a = *maxSepDist + 1;


initialError += pairSpikes(optOutputs[0], outputs[0], spikeDomAndPairs[0], a, maxSepDist, d, M);







FILE * unmatchedspikesf = fopen("unmatchedspikes.txt","w+");
FILE * spikechangesf = fopen("spikechanges.txt","w+");
FILE * spikeDomAndPairsf = fopen("spikeDomAndPairs.txt","w+");
FILE * errorProgressionf = fopen("errorprogression.txt", "w+");



int ** previousSpikeDomAndPairs = (int **) malloc(size*sizeof(int *));
for(i=0; i<size; i++)
{
	previousSpikeDomAndPairs[i] = (int *) malloc(2*numSteps*sizeof(int));
}

fprintf(spikeDomAndPairsf, "before optimization \n");


// begin optimization procedure; process will be repeated maxLoopCount times. 

while(loopcount < maxLoopCount)
{


	optUnpairedSpikeNum = 0;
	refUnpairedSpikeNum = 0;



	if(loopcount % maxSepDistDecrementer == 0 && *maxSepDist > 0)
	{
		*maxSepDist=(*maxSepDist)-1;
	}
	printf("loopcount: %d\n",loopcount);
	loopcount++;

	// simulate opt	

	simulate2(optOutputs, optNet, inputs, refstatemat, synCurrentRatio, synPercent);

	printf("simulated\n");
	// make a copy of outputs before calculating the total error and determining which weights need to be reduced/increased
       

	

	
	double toteError = 0;
	int * a = (int *) malloc(sizeof(int));

	// make a > maxSepDist
	*a = *maxSepDist + 1;
	for(i=0; i<size; i++)
	{
		toteError += pairSpikes(optOutputs[i], outputsCopy[i], spikeDomAndPairs[i], a, maxSepDist, d, M);
	}
	
	error = toteError;
	fprintf(errorProgressionf, "%lf\n",error);



	if(loopcount == maxLoopCount)
	{
		endError = error;	
	}

	printf("toteError at loop %d = %lf\n", loopcount,error);
	
	
	int optCounter = 0;

	int numWeightReductions = 0;
	int numWeightIncreases = 0;

	printf("starting alterations\n");
	for(i=0; i<size; i++)
	{

		optCounter = 0;
	//	refCounter = 0;		
		for(j=0; j<numSteps; j++)
		{
			while(optOutputs[i][optCounter] != 1 && optCounter < numSteps)
			{
				optCounter++;
			}
			if(optCounter >= numSteps)
			{
				break;
			}
			
				//unmatched spike in opt
			if(spikeDomAndPairs[i][j] == -1)
			{
				optUnpairedSpikeNum++;
				numWeightReductions++;
					// reduce presynaptic weights to neuron i, reduce a random synapse in the control network by same amount
				int x, y;
			
				for(x=0; x<size; x++)
				{
					if(optCounter > *z)
					{
						for(y = optCounter - *z; y < optCounter; y++)
						{
							if(optOutputs[x][y] == 1)
							{
															
								int rand1 = rand() % size;
								double q = epsilon *(rand() % 10);

                                                                optNet[i][x] -= q;
								
								int b;
								int xinhib = 0;

								for(b=0; b<numInhibitory; b++)
								{
									if(inhibIndices[b] == x)
									{
										xinhib = 1;
									}
								}								
		
								if(xinhib == 0 && optNet[i][x] < 0)
								{
									optNet[i][x] = 0;
								}
								if(xinhib == 1 && optNet[i][x] > 0)
								{
									optNet[i][x] = 0;
								}
																
									
								int rand2 = rand() % size;

								
								optNetControl[rand1][rand2] -= q;
								

								int rand2inhib = 0;
								int u;
								for(u=0; u<numInhibitory; u++)
								{
									if(inhibIndices[u] == rand2)
									{
										rand2inhib = 1;
									}
								}
								if(rand2inhib == 0 && optNetControl[rand1][rand2] < 0)
								{
									optNetControl[rand1][rand2] = 0;
								}
								if(rand2inhib == 1 && optNetControl[rand1][rand2] > 0)
								{
									optNetControl[rand1][rand2] = 0;
								}




								break;

							}
						}
					}

				}
					
					
			}
			

			optCounter++;	
		}	
	
		int k;
		numWeightIncreases++;
		int refCounter = 0;
		for(j=0; j<numSteps; j++)
		{
				// find the time of the jth firing event
			while(outputs[i][refCounter] != 1)
			{
				refCounter++;
				if(refCounter >= numSteps)
				{
					break;
				}
			}
			//unmatched spike in ref
			if(spikeDomAndPairs[i][numSteps+j] == -1)
			{
				refUnpairedSpikeNum++;
					//increase presynaptic weights to neuron i	
				int x, y;
				for(x=0; x<size; x++)
				{
					if(refCounter > *z)
					{
						for(y = refCounter - *z; y < refCounter; y++)
						{
							// this does not equally increase the presynaptic weights
							// if neuron x fires 4 times within the interval (optcounter -z, optcounter-1
							// then we will increase Wxi 4 times 
							if(optOutputs[x][y] == 1)
							{
						
								int rand1 = rand() % size;

								double q = epsilon * (rand() % 10);
								optNet[i][x] += q;
                                                                
								int b;
                                                                int xinhib = 0;

                                                                for(b=0; b<numInhibitory; b++)
                                                                {
                                                                        if(inhibIndices[b] == x)
                                                                        {
                                                                                xinhib = 1;
                                                                        }
                                                                }

                                                                if(xinhib == 0 && optNet[i][x] < 0)
                                                                {
                                                                        optNet[i][x] = 0;
                                                                }
                                                                if(xinhib == 1 && optNet[i][x] > 0)
                                                                {
                                                                        optNet[i][x] = 0;
                                                                }



								
								int rand2 = rand() % size;
				

								optNetControl[rand1][rand2] += q;



                                                                int rand2inhib = 0;
								int u;
                                                                for(u=0; u<numInhibitory; u++)
                                                                {
                                                                        if(inhibIndices[u] == rand2)
                                                                        {
                                                                                rand2inhib = 1;
                                                                        }
                                                                }
                                                                if(rand2inhib == 0 && optNetControl[rand1][rand2] < 0)
                                                                {
                                                                        optNetControl[rand1][rand2] = 0;
                                                                }
                                                                if(rand2inhib == 1 && optNetControl[rand1][rand2] > 0)
                                                                {
                                                                        optNetControl[rand1][rand2] = 0;
                                                                }

									
								break;
							}					
						}
					}
				}					
			}
			refCounter++;
		}


	//	fprintf(unmatchedspikesf, "refUnmatchedSpikeNum: %d\n",refUnpairedSpikeNum);		
		
	}

	fprintf(spikechangesf, "loopcount: %d\n", loopcount);
	for(i=0; i<size; i++)
	{
		for(j=0; j<numSteps; j++)
		{
			if(previousSpikeDomAndPairs[i][j] == -1 && spikeDomAndPairs[i][j] != -1)
			{
				fprintf(spikechangesf,"unmatched spike in opt eliminated | neuron: %d | spikenumber: %d\n",i,j);
			}
			if(previousSpikeDomAndPairs[i][j] != -1 && spikeDomAndPairs[i][j] == -1)
			{
				fprintf(spikechangesf, "new unmatched spike in opt | neuron: %d | spikenumber: %d\n", i, j);
			}	
		}
		for(j=numSteps; j<2*numSteps; j++)
		{
			if(previousSpikeDomAndPairs[i][j] == -1 && spikeDomAndPairs[i][j] != -1)
			{
                                fprintf(spikechangesf,"unmatched spike in ref eliminated | neuron: %d | spikenumber: %d\n",i,j-numSteps);
			}
			if(previousSpikeDomAndPairs[i][j] != -1 && spikeDomAndPairs[i][j] == -1)
			{
                                fprintf(spikechangesf, "new unmatched spike in ref | neuron: %d | spikenumber: %d\n", i, j-numSteps);
			}
		}
		
		for(j=0; j<2*numSteps; j++)
		{
			previousSpikeDomAndPairs[i][j] = spikeDomAndPairs[i][j];
		}
	}

	int optSpikeNum = 0;
	for(i=0; i<size; i++)
	{
		for(j=0; j<numSteps; j++)
		{
			optSpikeNum += optOutputs[i][j];
		}
	}

	
	printf("optSpikeNum: %d | original optSpikeNum %d | numSpikesRef: %d\n", optSpikeNum, originalSpikeNum, numSpikesOutputs);
	printf("optUnmatched: %d | refUnmatched %d \n", optUnpairedSpikeNum, refUnpairedSpikeNum);
	fprintf(unmatchedspikesf, "optUnmatched: %d | refUnmatched: %d \n", optUnpairedSpikeNum, refUnpairedSpikeNum);
	printf("loopy\n"); 


	int spikeNumDiff = numSpikesOutputs - optSpikeNum;

	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			
			optNet[i][j] += (rescale)*(spikeNumDiff);
                        optNetControl[i][j] += (rescale)*(spikeNumDiff);
			int u;

			int jInhib = 0;

			for(u=0; u<numInhibitory; u++)
			{
				if(inhibIndices[u] == j)
				{
						jInhib = 1;				
				}
			}
			if(jInhib == 1)
			{
				if(optNet[i][j] > 0)
				{
					optNet[i][j] = 0;
				}
				if(optNetControl[i][j] >0)
				{
					optNetControl[i][j] = 0;
				}
			}
			
			if(jInhib == 0)
			{
				if(optNet[i][j] < 0)
				{
					optNet[i][j] = 0;
				}
				if(optNetControl[i][j] < 0 )
				{			
					optNetControl[i][j] = 0;
				}
			}
			

		}
	}
	

	double optTotal = 0;
	double controlTotal = 0;
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			optTotal += optNet[i][j];
			controlTotal += optNetControl[i][j];
		}
	}

	printf("optTotal: %lf\n", optTotal);
	printf("controlTotal: %lf\n", controlTotal);


}


	fclose(errorProgressionf);

	fprintf(spikeDomAndPairsf, "after optimization!\n");
	for(i=0; i<size; i++)
	{
		for(j=0; j<2*numSteps; j++)
		{
			if(j<numSteps)
			{
				if(spikeDomAndPairs[i][j]==-1)
				{
					fprintf(spikeDomAndPairsf,"unmatched spike in opt at %d, %d\n", i, j);
				}
			}
			if(j >= numSteps)
			{
				if(spikeDomAndPairs[i][j] == -1)
				{
					fprintf(spikeDomAndPairsf, "unmatched spike in ref at %d, %d\n", i, j);
				}
			}
			
		}
	}



	fclose(unmatchedspikesf);
	fclose(spikeDomAndPairsf);
	fclose(spikechangesf);

	double absolute;



	int refSpikes, oldOptSpikes, newOptSpikes;
	refSpikes = 0;
	oldOptSpikes =0;
	newOptSpikes = 0;

	for(i=0; i<size; i++)
	{
		int countingSpikes = 0;
		for(j=0; j<numSteps; j++)
		{
			refSpikes += outputs[i][j];
			oldOptSpikes += rememberOptOutputs[i][j];
			newOptSpikes += optOutputs[i][j];
		}
	}	
	printf("refSpikes: %d | oldOptSpikes: %d | optSpikes: %d\n", refSpikes, oldOptSpikes, newOptSpikes);

	printf("Start Error = %d, Final Error = %d \n",startError,endError);

	FILE * fp;
	fp = fopen("newOptSpikes.txt","w+");
	for(i=0;i<size;i++)
	{
		for(j=0;j<numSteps;j++)
		{
			fprintf(fp,"%lf,",optOutputs[i][j]);
		}
		fprintf(fp,"\n");
	}	
	fclose(fp);

	FILE * fp1;
        fp1 = fopen("oldOptSpikes.txt","w+");
        for(i=0;i<size;i++)
        {
                for(j=0;j<numSteps;j++)
                {
                        fprintf(fp1,"%lf,",rememberOptOutputs[i][j]);
                }
                fprintf(fp1,"\n");
        }
        fclose(fp1);
	
	  FILE * fp2;
        fp2 = fopen("refSpikes.txt","w+");
        for(i=0;i<size;i++)
        {
                for(j=0;j<numSteps;j++)
                {
                        fprintf(fp2,"%lf,",outputs[i][j]);
                }
                fprintf(fp2,"\n");
        }
        fclose(fp2); 


	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
//			printf("optNet[%d][%d] = %lf\n", i, j, optNet[i][j]);
		}
	}



	FILE * fp3;
	fp3 = fopen("optweights.txt","w+");
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			fprintf(fp3,"%lf,",optNet[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fclose(fp3);

        FILE * fp4;
        fp4 = fopen("refweights.txt","w+");
        for(i=0;i<size;i++)
        {
                for(j=0;j<size;j++)
                {
                        fprintf(fp4,"%lf,",refNet[i][j]);
                }
                fprintf(fp4,"\n");
        }
        fclose(fp4);


        FILE * fp5;
        fp5 = fopen("controlweights.txt","w+");
        for(i=0;i<size;i++)
        {
                for(j=0;j<size;j++)
                {
                        fprintf(fp5,"%lf,",optNetControl[i][j]);
                }
                fprintf(fp5,"\n");
        }
        fclose(fp5);
	


	FILE * fs = fopen("syncurrentratio.txt","w+");
	FILE * fsp = fopen("synPercent.txt", "w+");
	
	for(i=0; i<size; i++)
	{
		fprintf(fs, "%lf,", synCurrentRatio[i]);
		fprintf(fsp, "%lf,", synPercent[i]);
	}


        simulate2(optOutputsCopy, optNetControl, inputs, refstatemat, synCurrentRatio, synPercent);

	FILE * fcontrol = fopen("controlSpikes.txt","w+");

	for(i=0; i<size; i++)
	{
		for(j=0; j<numSteps; j++)
		{
			fprintf(fcontrol, "%lf,", optOutputsCopy[i][j]);			
		}
		fprintf(fcontrol, "\n");
	}
	fclose(fcontrol);


	return 0;		
}


void loadNetwork(double ** network, int size)
{
	//Treat the network as a single dimensional array, 
	//loading in values like network[i*size + j] = ???
	   FILE * fp = fopen("weights.csv", "r");

        int currentrow, currentcolumn = 0;
        int counter = 0;

        int current = fgetc(fp);

        while(current != EOF)
        {
                char c = (char) current;

//		printf("%c ", c);
                if(c=='0')
                {
                        network[currentrow][currentcolumn] = 0;
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c == '1')
                {
                        network[currentrow][currentcolumn] = 1*weightScale;
	
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c=='2')
                {
                        network[currentrow][currentcolumn] = 2*weightScale;
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c == '3')
                {
                        network[currentrow][currentcolumn] = 3*weightScale;
	
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c=='4')
                {
                        network[currentrow][currentcolumn] = 4*weightScale;
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c == '5')
                {
                        network[currentrow][currentcolumn] = 5*weightScale;
	
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c=='6')
                {
                        network[currentrow][currentcolumn] = 6*weightScale;
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c == '7')
                {
                        network[currentrow][currentcolumn] = 7*weightScale;
	
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c=='8')
                {
                        network[currentrow][currentcolumn] = 8*weightScale;
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                if(c == '9')
                {
                        network[currentrow][currentcolumn] = 9*weightScale;
	
                        if(counter%size == 0)
                        {
                                currentrow = 0;
                                currentcolumn++;
                        }
                        else
                        {
                                currentrow++;
                        }
                        counter++;
                }
                
                current = fgetc(fp);
        }
        fclose(fp);
}

void loadSpikes(double ** inputSpikes)
{
	FILE * fp = fopen("spikes.csv","r");
	int currentrow = 0;
	int currentcolumn = 0;
	int counter = 0;
	
	int current = fgetc(fp);
	while(current != EOF && currentrow < size)
	{
		char c = (char) current;
		if(c=='0')
		{

			inputSpikes[currentrow][currentcolumn] = 0;
			if(counter%numSteps == 0)
			{
				currentcolumn = 0;
				currentrow++;
			}
			else
			{
				currentcolumn++;
			}
			counter++;
		}
		if(c=='1')
		{

			inputSpikes[currentrow][currentcolumn]= 1;
			if(counter%numSteps == 0)
			{
				currentcolumn = 0;
				currentrow++;
			}
			else
			{
				currentcolumn++;
			}
			counter++;
		}
		current = fgetc(fp);
	}
	fclose(fp);
}


//dynamic program for optimally pairing spikes from the reference and opt spike trains

double pairSpikes(double * optTrain, double * refTrain, double * maxSpikesAndPairs, int * a, int * maxSepDist, int ** d, int ** M)
{


	int i, j, k;
	int numSpikesO = 0;
	int numSpikesR = 0;

	
	//count the number of spikes in each spike train

	for(i=0; i<numSteps; i++)
	{
		if(optTrain[i] ==1.0)
		{
//			printf("still fine for optnet calculation, numSpikesO++");
			
			numSpikesO=numSpikesO+1;
		}		
		if(refTrain[i] == 1.0)
		{
			
			numSpikesR=numSpikesR+1;
		}
	}
//	printf("counted the number of spikes in each train");	
	
//	printf("numSpikesR: %d | numSpikesO: %d\n", numSpikesR, numSpikesO);

	int dummyspikes;
	double z;
	double q;
//	int Dz[numSpikesO];
//	double ans;

//	int ** d = (int **) malloc(numSpikesO*sizeof(int *));
//      for(i=0; i<numSpikesO; i++)
//        {
//                d[i] = (int *) malloc(numSpikesR*sizeof(int));
//
//        }

//        int ** M = (int **) malloc((numSpikesO+1)*sizeof(int *));
//        for(i=0; i<numSpikesO+1; i++)
//        {
//                M[i] = (int *) malloc((numSpikesR+1)*sizeof(int));
//        }

	
//	if(numSpikesR == 0)
//	{
//		return 0;
//	}



//if(numSpikesR > 0 && numSpikesO>0)
//{
	// ensure that optTrain is the spike train with more spikes
	
//	printf("optTrain and refTrain are set up now\n");
//	free(dummyspikes);
	



/*
	int ** d = (int **) malloc(numSpikesO*sizeof(int *));
	for(i=0; i<numSpikesO; i++)
	{
		d[i] = (int *) malloc(numSpikesR*sizeof(int));
	}



        int ** d = (int **) malloc(numSpikesO*sizeof(int *));
        for(i=0; i<numSpikesO; i++)
        {
                d[i] = (int *) malloc(numSpikesR*sizeof(int));

        }
//        printf("initialized d");

*/	
	



	//fill in values of d[i][j], where d[i][j] gives the temporal difference between the ith spike of optTrain and the jth spike of refTrain
	
	// oSpikeIndex and rSpikeIndex keep track of which spike I am looking at for a given spike train
	int oSpikeIndex= 0;
	int rSpikeIndex = 0;
	// to and tr represent the current time step that I am on within each spike train
	int tO=0;
	int tR=0;
	while(oSpikeIndex < numSpikesO)
	{
		//for this value of oSpikeIndex, reset rSpikeIndex to 0
		
		rSpikeIndex = 0;
		tR=0;
		//keep incrementing tO until the next spike in optTrain is discovered 
		while(optTrain[tO] != 1)
		{
			tO++;
			if(tO > numSteps-1)
			{
				
				break;
				
			}
		}
		
		//now tO is storing the value of the spike of interest
		//keep looking at the next spike in refTrain
		while(rSpikeIndex < numSpikesR)
		{
			//find the time of the next reTrain spike
			while(refTrain[tR] != 1)
			{
				tR++;
				if(tR > numSteps-1)
				{
					break;
				}
			}
		
			if(tO < numSteps && tR < numSteps)
			{
				d[oSpikeIndex][rSpikeIndex] = abs(tR-tO);			
				rSpikeIndex++;
				tR++;
			}
			else
			{
				break;
			}

		}

		// increment tO to prevent the same time step from being written to over and over again
		// increment oSpikeIndex because we have just accounted for this current spike of the spike train
		tO++;
		oSpikeIndex++;
			
	}

	
	// provide base cases for cost matrix, M
	
	
	for(i=0; i<numSpikesO+1; i++)
	{
		M[i][0] = *a * i;
	}
	for(j = 0; j<numSpikesR+1; j++)
	{
		M[0][j] = *a * j;
	}
	
	M[0][0] = 0;

	
	// compute cost matrix M
	
        for(i=1; i<numSpikesO+1; i++)
        {

                for(j=1; j<numSpikesR+1; j++)
                {
			int check = 0;
			while(check ==0)
			{
				
				int C1;
				if(d[i-1][j-1] <= *maxSepDist)
				{
                	        	C1 = d[i-1][j-1] + M[i-1][j-1];
				}
				else
				{
					C1 = 3* (*a) + M[i-1][j-1];	
				}
                	        int C2 = *a + M[i-1][j];
                	        int C3 = *a + M[i][j-1];
                	        int C4 = 2*(*a) + M[i-1][j-1];

                	        if(((C1<=C2 && C1<=C3) && C1<=C4)&&C1>=0)
                	        {
                	                M[i][j] = C1;
					check = 1;

                	        }
                	        else if(((C2<=C1 && C2<=C3) && C2 <= C4)&&C2>0)
                	        {
                	                M[i][j] = C2;
					check = 1;

                	        }
                	        else if(((C3<=C1 && C3 <= C2) && C3<=C4)&&C3>0)
                	        {
                	                M[i][j] = C3;
					check = 1;

                	        }
	        	       	else if(C4>0)
				{
                        
                	        	M[i][j] = C4;
					check = 1;

                	        }
			}
                }
        }


	int Dz[numSpikesO];
	int Dq[numSpikesR];

	for(i=0; i<numSpikesO; i++)
	{
		Dz[i] = 0;
	}
	for(j=0; j<numSpikesR; j++)
	{
		Dq[i] = 0;
	}


//	solving for D, the pairing of spike events

//	Trace back through M matrix
	i = numSpikesO;
	j = numSpikesR;

	while(i > 0 && j>0)
        {
		{
                        int C1;
			if(d[i-1][j-1] <= *maxSepDist)
			{
				C1 = d[i-1][j-1] + M[i-1][j-1];
			}
			else
			{
				C1 = 3*(*a) * ((numSteps + 1)) + M[i-1][j-1];
			}
                        int C2 = *a + M[i-1][j];
                        int C3 = *a + M[i][j-1];
              	        int C4 = 2*(*a) + M[i-1][j-1];
		
                        if((C1<=C2 && C1<=C3) && C1 <= C4)
                        {
                                Dz[i-1] = j-1;
				Dq[j-1] = i-1;
                                i--;
                                j--;
                        }
                        else if((C2<=C1 && C2<=C3) && C2<=C4)
                        {
                                Dz[i-1] = -1;
                                i--;
                        }
                        else if((C3<=C1 && C3 <= C2) && C3 <= C4)
                        {
                                Dq[j-1] = -1;
				j--;
                        }
			else
			{
				Dz[i-1] = -1;
				Dq[j-1] = -1;
				i--;
				j--;
			}   
		}
	}
	


	while(i>0)
	{
		Dz[i-1]= -1;
		i--;
	}
	while(j>0)
	{
		Dq[j-1] = -1;
		j--;
	}

	
	for(i=0; i<numSteps; i++)
	{
		if(i<numSpikesO)
		{
			maxSpikesAndPairs[i] = Dz[i];
		}
		else
		{
			maxSpikesAndPairs[i] = -5;
		}

	}
	
	for(j=0; j<numSteps; j++)
	{
		if(j<numSpikesR)
		{
			maxSpikesAndPairs[numSteps + j] = Dq[j];
		}
		else
		{
			maxSpikesAndPairs[numSteps + j] = -5;
		}
	}
		
	
	

	double ans = 0;

	ans = (double) (M[numSpikesO+1][numSpikesR+1]);
	
	return ans;
	

}


//	simulates network to produce output spike trains given a set of input stimulating currents and the weight matrix of the network

void simulate2(double ** outputs, double **  network, double ** inputs, double ** vstatemat, double * synCurrentRatio, double * synPercent)
{
        double vstate[size];
        int i, j;
        for(i=0; i<size; i++)
        {
                vstate[i] = 0;
        }
	double synContributions[size];
	double currContributions[size];
	for(i=0; i<size; i++)
	{
		synContributions[i] = 0;
		currContributions[i] = 0;
	}
        double temp = 0;
        double * currentmatmult = (double *)malloc(size*sizeof(double));
        for(i = 0; i<numSteps; i++)
        {
                for(j=0; j<size; j++)
                {
                        if(vstate[j]<-15)
			{
                                vstate[j]=0;
                        }
                        if(vstate[j]>voltageThreshold)
                        {

                                currentmatmult[j] = 1;
                        }
                        else
                        {
                                currentmatmult[j] = 0;
                        }
                }
                matrixMult(network, currentmatmult);
                for(j =0; j<size; j++)
                {


                        if(vstate[j]>voltageThreshold)
                        {
                                outputs[j][i]=1;
                                vstate[j] = 0;
				vstate[j] += Tau*step*(currentmatmult[j] + R* inputs[j][i]);
			//	synContributions[j] += abs(step*(currentmatmult[j]));
			//	currContributions[j] += abs(step*(R*inputs[j][i]));
                        }
                        else
                        {
                                outputs[j][i] = 0;
                                temp =  (step)*((-1)*vstate[j] + currentmatmult[j] + R*(inputs[j][i]));
                                temp = temp*Tau;


                                vstate[j] = vstate[j] + temp;
                                temp = 0;
				
				
                        }
                        synContributions[j] += sqrt(pow((step*(currentmatmult[j])),2));
                        currContributions[j] += sqrt(pow((step*(R*inputs[j][i])),2));

                        vstatemat[j][i] = vstate[j];
                }

	
        }
	
	for(i=0; i<size; i++)
	{
		if(currContributions[i] != 0)
		{
			synCurrentRatio[i] = synContributions[i] / currContributions[i]; 
		}
		else
		{
			synCurrentRatio[i] = -1;
		}
		
		synPercent[i] = synContributions[i] / (synContributions[i] + currContributions[i]);
	}
        free(currentmatmult);
}




void simulate1(double ** outputs, double **  network, double ** inputs, double ** vstatemat)
{
        double vstate[size];
        int i, j;
        for(i=0; i<size; i++)
        {
                vstate[i] = 0;
        }
        double temp = 0;
        double * currentmatmult = (double *)malloc(size*sizeof(double));
        for(i = 0; i<numSteps; i++)
        {
                for(j=0; j<size; j++)
                {
                        if(vstate[j]<-15){
                                vstate[j]=-15;
                        }
                        if(vstate[j]>voltageThreshold)
                        {

                                currentmatmult[j] = 1;
                        }
                        else
                        {
                                currentmatmult[j] = 0;
                        }
                }
                matrixMult(network, currentmatmult);
                for(j =0; j<size; j++)
                {


                        if(vstate[j]>voltageThreshold)
                        {
                                outputs[j][i]=1;
                                vstate[j] = 0;
				vstate[j] += step*(currentmatmult[j] + R* inputs[j][i]);
                        }
                        else
                        {
                                outputs[j][i] = 0;
                                temp =  (step)*((-1)*vstate[j] + currentmatmult[j] + R*(inputs[j][i]));
                                temp = temp*Tau;
                                vstate[j] = vstate[j] + temp;
                                temp = 0;
                        }
                        vstatemat[j][i] = vstate[j];
                }
        }
        free(currentmatmult);
}






void matrixMult(double ** A, double * B)
{
        int i, k;
        double row = 0;
        double entry = 0;
        double * C = (double *) malloc(size*sizeof(double));
        for(i = 0; i < size; i++)
        {
                double row = 0;
                for (k=0;k<size;k++)
                {
                        double entry = 0;
                        entry = (A[i][k])*(B[k]);
                        row = row + entry;
                }
                C[i] = row;
        }

        for(i=0; i<size; i++)
        {
                B[i] = C[i];
        }
        free(C);
	
}

