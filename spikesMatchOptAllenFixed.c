#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
typedef int bool;
#define true 1
#define false 0


//specifying global parameters

//membrane resistance
double R = 100000000;

//for generating random intiial weights

int scaler = 1000;

//number of iterations

int numIter = 100;

//mean input stimulating current

double inputMean = 0;

//standard deviation of input stimulating currents
double inputVar =  .00000000025;

//network size

int size = 383;

//voltage threshold in volts
double voltageThreshold = .03;

// reciprocal of membrane time constant in Hz
double Tau = 30.0;

// time step size (determined by 
double step = 0.025;

// number of time steps in the simulation
int numSteps = 5000;

//also for generating random numbers

double factor = .000001;


// maximum number of optimization trials
int maxLoopCount = 100;

//scaling input weighhts if reading in initial network structure
double weightScale = 3;

// maximum size of individual weight adjustment
double epsilon = .0000001;

// for synaptic rescaling
double rescale = .00000001;

// max number of timesteps back considered when determining 
//int Z = 5;


// FUNCTION DECLARATIONS


// simulates network with given input currents, producing output spike trains, also passes back voltage states for each neuron (knowledge of which is not used in optimization) the ratio of synaptic current contributions to external current contributions to voltage changes).
void simulate2(double ** outputs, double ** network, double ** inputs, double ** vstatemat, double * synCurrentRatio, double * synPercent);
void simulate1(double ** outputs, double ** network, double ** inputs, double ** vstatemat);

// premultiplies vector B with matrix A: AB
void matrixMult(double ** A, double * B);

// dynamic program that optimally pairs individual spikes given a reference spike train and a spike train from the current network. 
double pairSpikes(double * optTrain, double * refTrain, double * maxSpikesAndPairs, int * a, int * maxSepDist, int ** d, int **M);

// for loading in a preassigned initial network structure of specified size.
void loadNetwork(double ** network, int size);

// loads reference spike trains from biological recordings. 
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

	


	// initial specification of optNet. NOTE: other specifications can be used.

	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			optNet[i][j] = (rand() % scaler) * factor;
			
			if((rand() % 100) < 20)
			{
				optNet[i][j] = -optNet[i][j];
			}
			

//			if((rand() % 100) < 40)
//			{
//				optNet[i][j] = 0;
//			}
		


			// instantiate control network with same initial weights as optNet


			optNetControl[i][j] = optNet[i][j];
		}
	}




	// Gaussian Inputs; distribution generated from Box-Muller transform.

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
	

	//load in the objective outputs from file in directory

	loadSpikes(outputs);

	int numSpikesOutputs = 0;

	for(i=0; i<size; i++)
	{
		for(j=0; j<numSteps; j++)
		{
			numSpikesOutputs += outputs[i][j];
		}
	}
	

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


// maximum number of timesteps back that presynaptic firing events are considered as causally impacting a spiking event at time t.
int * z = (int *)malloc(sizeof(int)); 
*z = 20;

// upper bound on distance that spikes can be matched
int * maxSepDist = (int *) malloc(sizeof(int));
*maxSepDist = 15;

int optUnpairedSpikeNum = 0;
int refUnpairedSpikeNum = 0;

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
//		printf("optOutputs[%d][%d] = %lf \n", i, j, optOutputs[i][j]);
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




//begin optimization procedure



while(loopcount < maxLoopCount)
{

	
	//number of unpaired spikes in opt spike trains and those in ref.

	optUnpairedSpikeNum = 0;
	refUnpairedSpikeNum = 0;


	// as the procedure continues, reduce the maximum separation distance
	if(loopcount % 15 == 0 && *maxSepDist > 0)
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
	

	// use dynamic program to optimally pair spikes from the opt and ref spike trains. store result in spikeDomAndPairs

	for(i=0; i<size; i++)
	{

		toteError += pairSpikes(optOutputs[i], outputsCopy[i], spikeDomAndPairs[i], a, maxSepDist, d, M);//size;
		

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

	
	// make alterations to optNetwork weight matrix
	
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
					// reduce presynaptic weights to neuron i. Reduce a random synapse in control network.
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

								int rand2 = rand() % size;
								optNetControl[rand1][rand2] -= q;
								
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
								
								int rand2 = rand() % size;
				

								optNet[rand1][rand2] += q;									
								break;
							}					
						}
					}
				}					
			}
			refCounter++;
		}


			
		
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
			
		}
	}
	


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
	for(i=0;i<size;i++)
        {
                for(j=0;j<size;j++)
                {
			absolute = fabs(optNet[i][j]);		
                        if(absolute<0.01)
                        {
                                optNet[i][j] = 0.0;
                        }
                }
        }


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

	FILE * fp3;
	fp3 = fopen("weightsOut.txt","w+");
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			fprintf(fp3,"%lf,",optNet[i][j]);
		}
		fprintf(fp3,"\n");
	}
	fclose(fp3);


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



// can be used to load a size x size weight matrix saved as 'weights.csv' in the same directory

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


// dynamic string matching program

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
		
			numSpikesO=numSpikesO+1;
		}		
		if(refTrain[i] == 1.0)
		{
			
			numSpikesR=numSpikesR+1;
		}
	}

	int dummyspikes;
	double z;
	double q;


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


	
	for(i=0; i<numSpikesO+1; i++)
	{
		M[i][0] = *a * i;
	}
	for(j = 0; j<numSpikesR+1; j++)
	{
		M[0][j] = *a * j;
	}

	// base case for cost matrix
	M[0][0] = 0;


	// populate cost matrix
	//
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


	//Traceback to solve for D, the pairing of spike events
	
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


// simulates a network with provided weight matrix and input currents to generate output spike trains


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
//                                printf("currentmatmult[%d] = %lf\n", j, currentmatmult[j]);
                                outputs[j][i] = 0;
                                temp =  (step)*((-1)*vstate[j] + currentmatmult[j] + R*(inputs[j][i]));
                                temp = temp*Tau;
//                                printf("currentmatmult: %lf; R*inputs[%d][%d]: %lf, vstate[%d]=%lf", currentmatmult[j], i, j, R*inputs[j][i], j, vstate[j]);
//                        printf("Temp: %lf\n", temp);
                                vstate[j] = vstate[j] + temp;
                                temp = 0;
                        }
 //                       printf("vstate[%d] = %lf\n",j, vstate[j]);
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

