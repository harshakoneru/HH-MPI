/*
This code generates a small world network with N nodes with k average connections per node. It then simulates N stochastic HH neurons connected in such a topology using the 4th order Runge - Kutta method for t_max seconds. 
It calculates Mean ISI, Firing Correlation, Avg Firing correlation etc.
*/
#include <time.h>
#include <stdio.h>
#include <math.h> 
#include <stdlib.h>

//#define N 100 //Define the number of nodes in the network
#define PI 3.1415926536

const int N = 100;
const int k = 8;   //Define the mean degree of the network
const double timewell = 1.0;  // ????????????
const double time_max = 5000.0;   

double xna = 0.9; 
double xk = 0.6;
double rhona = 60.0; //um-2
double rhok = 18.0; //um-2
double S = 1.0; //Membrane patch size
double Nna;
double Nk;
double  GNa=120.0;
double	GK = 36.0;
double	GL=0.3;
double	VNa = 115.0;
//double VNa = 50.0; //mV (This is according to the paper)
double	VK = -12.0;
//double VK = -77.0; //mV (This is according to the paper)
double	VL=10.5995;
//double VL = -54.4; //mV (This is according to the paper)
double	c = 1.0;
int n_one_cycle = 1000;  // ????????????????
double D  = 0.001;
double firing_coherence = 0.0;
//float *voltagematrix, **tspike, *timematrix,*vjunk,*tjunk;

double KV( double C, double V,double m,double h,double n);
double Km( double V , double m );
double Kh( double V , double h );
double Kn( double V , double n );
double Kcoupling1(double a[N], double b[][N], int m, double k[N], int index);

int randomnumbergenerator(int max_num);
double gaussrandm();
double gaussrandh();
double gaussrandn();
double clusteringcoefficient(int a[][N], int m);
double meanshortest(int a[][N], int m);


int main()
{	

 int nodes = (int)N;
 double p = 0.05; //Rewiring probablity
 float randomnumber; //Float because we need to calculate probability
 int rewiringnumber, no_of_connections; /*rewiringnumber is the node to which a broken connection is rewired. 
 no_of_connections is the variable used to count the number of connections on each node */
 int max_randomnumber = 100000; //Maximum Random Number
 int degree[20];//Array containing number of neurons with n connections. For example, degree[2] will give number of neurons with 2 connections.
 srand(time(NULL));  //Seeding random number generator with current time
 int no_spikes[N];

//Defining and initializing the Matrix containing edges


double edgematrix[N][N];   // ????
double edgematrix1[N][N]; // ????

FILE *firingcoherenceid;
firingcoherenceid = fopen("FiringCoherenceVaryingDn0(100,8) (0.9; 0.6).txt","w"); // ????

FILE *fcdegreeid;
fcdegreeid = fopen("FCDegreeVaryingDsn0(100,8) (0.9; 0.6).txt","w");// ????

FILE *fisiid;
fisiid = fopen("ISI_no_connectionsn0(100,8) (0.9; 0.6).txt","w");// ????

FILE *fmeanisiid;
fmeanisiid = fopen("MeanISIn0(100,8) (0.9; 0.6).txt","w"); // ????

/*FILE *fspikecountid;
fspikecountid = fopen("No_of_spikes(100,8)(1.0,0.7).txt","w") // ????
*/

for(D = 0.1; D <= 0.2; D = D + 0.05)
{
	printf("D = %f\n",D);


	double tot_degree[12]; // ??????
	int count_degree[12];    // ???
	int i,j,l; //Looping variables. Not using k because it is already used above.
	for( i = 0; i < N; i++) // i ---> number of neurons
	{
	   no_spikes[i] = 0;  // 
	}

	for(i =0; i < 12; i++) // i --> degree of separation
	{ 
		tot_degree[i] = 0;  // ???
		count_degree[i] = 0;  // ???? , i ???
	}

	double FinalISI[N];
	double connections[N];


	for(i = 0; i < N; i++)
	{
		FinalISI[i] = 0.0;
		connections[i] = 0;
	}

	int trial = 0;

	for(trial = 0; trial < 10;  trial++)

	{

		printf("Trial %d\n",trial+1);

		//Initializing file

		 FILE *fV_vs_t_id;
		 fV_vs_t_id = fopen("V_vs_T.txt","w");

		 FILE *fspiketimeid;
		 fspiketimeid = fopen("Spiketime.txt","w");

		 FILE *yid;
		 yid = fopen("Y.txt","w");

		 FILE *fVTid;
		 fVTid = fopen("VT.txt","w");

		 FILE *fAdjacencymatrixId;
		 fAdjacencymatrixId = fopen("RewiredMatrix.txt", "w");

		 FILE *frandomnumberId;
		 frandomnumberId = fopen("Randomnumbers.txt", "w");

		 FILE *fNoConnectionsId;
		 fNoConnectionsId = fopen("NoOfConnections.txt","w");

		 FILE *fmnoiseId;
		 fmnoiseId = fopen("mnoise.txt","w");

		 FILE *fmvalueId;
		 fmvalueId = fopen("mvalues.txt","w");

		  FILE *fhvalueId;
		 fhvalueId = fopen("hvalues.txt","w");

		  FILE *fnvalueId;
		 fnvalueId = fopen("nvalues.txt","w");



  		printf("p = %f\n",p);
		int Done = 0; //Flag


		//Initializing degree array
		for(i = 0; i<20; i++)
		{
			degree[i] = 0;
		}
 
		//Initializing Edgematrix
		for(i = 0; i < N; i++)
 		{
			for (j = 0; j < N; j++)
  			{
  	  		edgematrix[i][j] = 0;	
  	  		edgematrix1[i][j] = 0;
  			}	
 		}
		Nna = rhona * S;
		Nk = rhok * S;
 		//Initializing V[j] Array
/*____________________________________________________________________________________________________________________________*/

  //Constructing a ring lattice with k/2 connections on either side of each node
 		for(i = 0; i < N; i++)
  		{  	
  			for(j = 1; j <= k/2; j++)
  			{

  			if((i+j) >= N)
  			{
  				edgematrix[i][abs(N - (i+j))] = 1; 
  				edgematrix[i][abs(i - j)] = 1;
  				edgematrix[abs(N - (i+j))][i] = 1; //Since for an undirectional graph, edgematrix[i][j] = edgematrix[j][i]
  				edgematrix[abs(i - j)][i] = 1;
  			}
  		        else if ((i-j) < 0)
  			{
  				edgematrix[i][abs(N + (i-j))] = 1;
  				edgematrix[i][abs(i+j)] = 1;
  				edgematrix[abs(N + (i-j))][i] = 1; //Since for an undirectional graph, edgematrix[i][j] = edgematrix[j][i]
  				edgematrix[abs(i+j)][i] = 1;
  			}
  			else
  			{
  				edgematrix[i][abs(i + j)] = 1;
  		    		edgematrix[i][abs(i - j)] = 1;
  		    		edgematrix[abs(i + j)][i] = 1; //Since for an undirectional graph, edgematrix[i][j] = edgematrix[j][i]
  		    		edgematrix[abs(i - j)][i] = 1;
  			}
  			}
		 }


		//printf("Finished printing. \n");
		//printf("Generic ring lattice with k/2 connections on either side of each node generated. \n");

		for(i = 0; i< N; i++)
		{
			for(j = 0; j<N; j++)
			{
			fprintf(fAdjacencymatrixId,"%f   ",edgematrix[i][j]);
			edgematrix1[i][j] = edgematrix[i][j];
			}
			fprintf(fAdjacencymatrixId, "\n");
		}

		fprintf(fAdjacencymatrixId, "\n\n\n\n");
/*______________________________________________________________________________________________________________________________________*/

		//Rewiring the system using the rewiring probability

		fprintf(frandomnumberId, "#N = %d, k = %d, p = %f \n", nodes, k , p);
		fprintf(frandomnumberId, "S.No i    j    r    re   \n");

		for(i = 0; i < N-1; i++)
		{
			for (j = i+1; j < N; j++)
			{
			if (edgematrix[i][j] == 1 )
			{
			randomnumber = randomnumbergenerator(max_randomnumber); //Generate a random number and divide it by max_random number.
			randomnumber = (float)randomnumber/(float)max_randomnumber; //If normalized value is less than p, remove connection
			if (randomnumber < p)
			{
				
			edgematrix[i][j] = 0;
			edgematrix[j][i] = 0;
			Done = 0;
			while (Done != 1)  // done ??????
			{
				rewiringnumber = randomnumbergenerator(N-1);
				if (rewiringnumber != i && edgematrix[i][rewiringnumber] == 0) //Eliminating self-loops& redundant connections
				{
				edgematrix[i][rewiringnumber] = 1;
				edgematrix[rewiringnumber][i] = 1; //Since for an undirectional graph, edgematrix[i][j] = edgematrix[j][i]
				fprintf(frandomnumberId, "%d    %d    %d    %f    %d\n",0,i,j,randomnumber,rewiringnumber);
				Done = 1;
						
				}
				else
				continue;
			}		
			}		
			}  
			else
			continue;	
			}
	//printf("Node %d rewired. \n", i+1);	
		}



/*________________________________________________________________________________________________________________________________*/
		//Writing into files and printing
		//printf("Starting to write adjacency matrix into the file. \n\n");
		fprintf(fAdjacencymatrixId, "#N = %d, k = %d, p = %f \n", nodes, k, p);
		fprintf(fNoConnectionsId, "#N = %d, k = %d, p = %f \n", nodes, k , p);
		for(i = 0; i < N; i++)
		{
	
		no_of_connections = 0;
		for(j = 0; j < N; j++)
		{
			if (edgematrix[i][j] == 1)
			{
			edgematrix[i][j] = 1;
			}
			fprintf(fAdjacencymatrixId, "%f   ", edgematrix[i][j]); //Printing Adjacency Matrix of neuron
			if(edgematrix[i][j] == 1)
			{
			no_of_connections = no_of_connections + 1;	//Calculating number of connections of current neuron
			}	
		}
		degree[no_of_connections] = degree[no_of_connections] + 1;
		connections[i] = connections[i] + no_of_connections;
		fprintf(fAdjacencymatrixId, "\n" );  
		fprintf(fNoConnectionsId, "%d    %d\n",i+1,no_of_connections);
		}
		//printf("Finished writing adjacency matrix into the file. \n\n");
		//printf("Finished printing number of connections per neuron. \n\n");

		for(i = 0; i < N; i++) // ????
		{
		connections[i] = connections[i]/2.0; // ???
		}

		fprintf(fAdjacencymatrixId, "\n\n\n\n");
		int countno = 0;
		for(i = 0; i < N-1; i++)
		{
			for(j = i+1; j<N; j++)
			{
			if(edgematrix[i][j] == 0 && edgematrix1[i][j] == 1)
			{
				countno++;
			}
			}
		}

		fprintf(fAdjacencymatrixId, "%d",countno);


/*______________________________________________________________________________________________________________________________________*/
		//Defining parameters related to V[j] and gating variables
		double sum_coupling = 0.0;
		double V[N];
		double m[N];
		double n[N];
		double h[N];
		double coupling[N];
		double alpham[N];
		double betam[N];
		double alphah[N];
		double betah[N];
		double alphan[N];
		double betan[N];
		double k1[N],m1[N],n1[N],h1[N],coupling1[N];
		double k2[N],m2[N],n2[N],h2[N],coupling2[N];
		double k3[N],m3[N],n3[N],h3[N],coupling3[N];
		double k4[N],m4[N],n4[N],h4[N];
		double dt = 0.01; //ms
		//double time_max = 100.0; //ms
		int n_steps = (int)time_max/dt;
		printf("No of time steps: %d\n",n_steps);
		double mnoise[N],nnoise[N],hnoise[N];

		 for(i = 0; i < N; i++)
		 {
		 	V[i] =  0.0;
		 	k1[i] = 0.0;
		 	k2[i] = 0.0;
		 	k3[i] = 0.0;
		 	k4[i] = 0.0;
		 	m1[i] = 0.0;
		 	m2[i] = 0.0;
		 	m3[i] = 0.0;
		 	m4[i] = 0.0;
		 	n1[i] = 0.0;
		 	n2[i] = 0.0;
		 	n3[i] = 0.0;
		 	n4[i] = 0.0;
		 	h1[i] = 0.0;
		 	h2[i] = 0.0;
		 	h3[i] = 0.0;
		 	h4[i] = 0.0;
		 	coupling1[i] = 0.0;
		 	coupling2[i] = 0.0;
		 	coupling3[i] = 0.0;
		 	mnoise[i] = 0.0;
		 	nnoise[i] = 0.0;
		 	hnoise[i] = 0.0;
 	
		 	alpham[i] = (0.1*(25-V[i]))/(exp((25 - V[i])/10) - 1);
			betam[i] = 4*exp(-V[i]/18);
			alphah[i] = 0.07*exp(-V[i]/20);
			betah[i] = 1/(exp((30-V[i])/10) + 1);
			alphan[i] = 0.01*(10-V[i])/(exp((10-V[i])/10) - 1);
			betan[i] = 0.125*exp(-V[i]/80);

			/*alpham[i] = (0.1*(V[i] + 40))/(-exp((-40 - V[i])/10) + 1);
			betam[i] = 4*exp((-V[i] - 65)/18);
			alphah[i] = 0.07*exp((-V[i]-65)/20);
			betah[i] = 1/(exp((-35-V[i])/10) + 1);
			alphan[i] = 0.01*(V[i] + 55)/(-exp((-55-V[i])/10) + 1);
			betan[i] = 0.125*exp((-V[i]-65)/80); */

			m[i] = (alpham[i]/(alpham[i] + betam[i]));
			h[i] = (alphah[i]/(alphah[i] + betah[i]));
			n[i] = (alphan[i]/(alphan[i] + betan[i]));
		 }
		 
		for (i = 0; i < N; i++)
		{
			sum_coupling = 0.0;
			for (j = 0; j < N; j++)
			{
				sum_coupling = sum_coupling + (edgematrix[i][j]*(V[j] - V[i]));
			}
			coupling[i] = sum_coupling;
		}

		//Starting calculation
		for(l = 0; l <= n_steps; l++)
		{
			//printf("%f\n",l*dt );
			fprintf(fV_vs_t_id, "%f\t",l*dt);
			fprintf(fmvalueId, "%f\n",l*dt);
			fprintf(fhvalueId, "%f\n",l*dt);
			fprintf(fnvalueId, "%f\n",l*dt);
			fprintf(fmnoiseId, "%f\t",l*dt);
			for(j = 0; j < N; j++)
			{
			alpham[j] = (0.1*(25-V[i]))/(exp((25 - V[i])/10) - 1);
			betam[j] = 4*exp(-V[i]/18);
			alphah[j] = 0.07*exp(-V[i]/20);
			betah[j] = 1/(exp((30-V[i])/10) + 1);
			alphan[j] = 0.01*(10-V[i])/(exp((10-V[i])/10) - 1);
			betan[j] = 0.125*exp(-V[i]/80);

			mnoise[j] =  sqrt(((2*alpham[j]*betam[j])/(Nna*xna*(alpham[j] + betam[j]))))*gaussrandm();
			fprintf(fmnoiseId, "%f\t",mnoise[j]);
			nnoise[j] =  sqrt(((2*alphan[j]*betan[j])/(Nk*xk*(alphan[j] + betan[j]))))*gaussrandn();
			hnoise[j] =  sqrt(((2*alphah[j]*betah[j])/(Nna*xna*(alphah[j] + betah[j]))))*gaussrandh();
			k1[j] = dt*KV(coupling[j],V[j],m[j],h[j],n[j]);
		   	m1[j] = dt*(Km(V[j],m[j]) + mnoise[j]);
		    h1[j] = dt*(Kh(V[j],h[j])+ hnoise[j]);
		    n1[j] = dt*(Kn(V[j],n[j])+ nnoise[j]);
		  	coupling1[j] = dt*Kcoupling(V,edgematrix,nodes,k1,j);
		}
		fprintf(fmnoiseId, "\n");
		for(j = 0; j <N; j++)
		{

		    k2[j] = dt*KV(coupling[j]+(0.5*coupling1[j]),V[j]+(0.5*k1[j]),m[j]+(0.5*m1[j]),h[j]+(0.5*h1[j]),n[j]+(0.5*n1[j]));
		    m2[j] = dt*(Km(V[j]+(0.5*k1[j]),m[j]+(0.5*m1[j]))+mnoise[j]);
		    h2[j] = dt*(Kh(V[j]+(0.5*k1[j]),h[j]+(0.5*h1[j]))+hnoise[j]);
		    n2[j] = dt*(Kn(V[j]+(0.5*k1[j]),n[j]+(0.5*n1[j]))+nnoise[j]);
		    coupling2[j] = dt*Kcoupling(V,edgematrix,nodes,k2,j);
		}

		for(j = 0; j<N; j++)
		{
	
		    k3[j] = dt*KV(coupling[j]+(0.5*coupling2[j]),V[j]+(0.5*k2[j]),m[j]+(0.5*m2[j]),h[j]+(0.5*h2[j]),n[j]+(0.5*n2[j]));
		    m3[j] = dt*(Km(V[j]+(0.5*k2[j]),m[j]+(0.5*m2[j]))+mnoise[j]);
		    h3[j] = dt*(Kh(V[j]+(0.5*k2[j]),h[j]+(0.5*h2[j]))+hnoise[j]);
		    n3[j] = dt*(Kn(V[j]+(0.5*k2[j]),n[j]+(0.5*n2[j]))+nnoise[j]);
		    coupling3[j] = dt*Kcoupling1(V,edgematrix,nodes,k3,j);
		}

		for (j = 0; j < N; j++)
		{

			k4[j] = dt*KV(coupling[j]+coupling3[j],V[j]+k3[j],m[j]+m3[j],h[j]+h3[j],n[j]+n3[j]);
		        m4[j] = dt*(Km(V[j]+k3[j],m[j]+m3[j])+mnoise[j]);
		   	h4[j] = dt*(Kh(V[j]+k3[j],h[j]+h3[j])+hnoise[j]);
		   	n4[j] = dt*(Kn(V[j]+k3[j],n[j]+n3[j])+nnoise[j]);
		}

		for (j = 0; j < N; j++)
		   	{
		    		V[j] = V[j] + ((k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6);
		    		m[j] = m[j] + ((m1[j] + 2*m2[j] + 2*m3[j] + m4[j])/6);
		   		h[j] = h[j] + ((h1[j] + 2*h2[j] + 2*h3[j] + h4[j])/6);
		   		n[j] = n[j] + ((n1[j] + 2*n2[j] + 2*n3[j] + n4[j])/6);
				fprintf(fV_vs_t_id, "%f\t",V[j]);
				fprintf(fVTid, "%f\t", l*dt);
				fprintf(fVTid, "%f\n", V[j]);
				fprintf(fmvalueId, "%f\t",m[j]);
				fprintf(fhvalueId, "%f\t",h[j]);
				fprintf(fnvalueId, "%f\t",n[j]);
			}	
			fprintf(fV_vs_t_id, "\n");
			fprintf(fmvalueId, "\n");
			fprintf(fhvalueId, "\n");
			fprintf(fnvalueId, "\n");
		for (i = 0; i < N; i++)
		{
			sum_coupling = 0.0;
			for (j = 0; j < N; j++)
			{
				sum_coupling = sum_coupling + (edgematrix[i][j]*(V[j] - V[i]));
			}
			coupling[i] = sum_coupling;
		}

	}


	fclose(fAdjacencymatrixId);
	fclose(fV_vs_t_id);
	fclose(fVTid);
	fclose(fNoConnectionsId);
	fclose(frandomnumberId);
	fclose(fmvalueId);
	fclose(fnvalueId);
	fclose(fhvalueId);
	//fclose(fMeanShortestPathId);
	//fclose(fClusteringCoeffId);
//________________________________________________________________________________________________________________________________________

//Calculation of ISI, spike times, firing correlation, average firing coherence etc.
	double f_correlation[N][N];
	float tmp0 = 0.0;
	float tmp1 = 0.0;
	int s = 0;
	int a = 0;
	int b =0;
	int c = 0;
	double z = 0.0;
	int smax = 0;
	char ignore[1024];
	double t1 = 0.0;
	double t2 = timewell;
	int no_cycles = (int)(time_max/timewell);
	double sum_total = 0.0;

	for(i = 0; i < N; i++)
	{
		//no_spikes[i] = 0;
		for(j = 0; j < N; j++)
		{
			f_correlation[i][j] = 0.0;
		}
	}

	//float *vjunk;
	//float *tjunk;
	//int countnotimes = 0;
	fVTid = fopen("VT.txt","r");
	//voltagematrix  = voltagematrix = (float *)(malloc(n_steps*sizeof(float)));
	float voltagematrix[n_steps];
	float timematrix[n_steps];
	for(i = 0; i < n_steps; i++)
	{
		voltagematrix[i] = 0.0;
		timematrix[i] = 0.0;
	}	
/*
	tspike = (float **)malloc((n_steps/n_one_cycle)*sizeof(float));
	for(c = 0; c < (n_steps); c++)
	{
	tspike[c] = (float *)malloc((n_steps/n_one_cycle) * sizeof(float));
	}
*/
	float	tspike[N][n_steps/n_one_cycle];
	for (i=0;i<N;i++)
		for (j=0;j< (n_steps/n_one_cycle);j++)
			tspike[i][j]=0;

//	timematrix = (float *)malloc((n_steps)*sizeof(float));

	//vjunk = (float *)malloc((n_steps)*sizeof(float));
	//tjunk = (float *)malloc((n_steps)*sizeof(float));

	for(i = 0; i < N; i++)
	{
		//countnotimes = 0;
		s = 0;
		tmp0 = 0.0;
		tmp1 = 0.0;

		for(j = 0; j < i; j++)
		{
			fgets(ignore, sizeof(ignore), fVTid);
		}
		for(l = 0; l < n_steps; l++)
		{
			fscanf(fVTid,"%f %f",&timematrix[l],&voltagematrix[l]); // reading from file
			for(a = 0; a < N; a++)
			{
				fgets(ignore, sizeof(ignore), fVTid);
			}
		}
		//printf("%f\n",voltagematrix[1500]);
		for(b = 1; b <= n_steps - 1; b++)
		{
			//printf("Here\n");
			//countnotimes = countnotimes + 1;
			tmp0 = voltagematrix[b] - voltagematrix[b+1];
			tmp1 = voltagematrix[b] - voltagematrix[b-1];
			if(tmp0 > 0 && tmp1 > 0)
			{
				if(voltagematrix[b] > 70.0)
				{
					tspike[i][s] = timematrix[b];
					if(floor(tspike[i][s]) == floor(tspike[i][s-1]) || floor(tspike[i][s]) == ceil(tspike[i][s-1]))
						{
						tspike[i][s] = 0.0;
						continue;
						}
					fprintf(fspiketimeid,"%f\t",tspike[i][s]);
					s = s + 1;
					no_spikes[i] = s;
					if(s > smax)
					{
						smax = s;
					}
				}
			
			}	
		}
		//printf("%d\n",countnotimes);
		fprintf(fspiketimeid,"\n");
		rewind(fVTid);
	//	printf("No of spikes of neuron %d: %d\n",i+1,no_spikes[i]);
	}

	int flag  = 0;
	for(i = 0; i < N; i++)
	{
		if(no_spikes[i] == 0)
			flag  = 1;		  
	}

	if(flag == 1)
	{
	trial = trial - 1;
	printf("Re-simulating, no spikes detected\n");
	continue;
	}

	double AvgISI[N];

	double ISI = 0.0;
	int count_isi = 0;
	for(i = 0; i < N; i++)
	{
		AvgISI[i] = 0.0;
	}

	for(i = 0; i < N; i++)
	{
		count_isi = 0;	
		for(j = 1; j < smax; j++)
		{
			if(tspike[i][j] == 0.0)
			{
				continue;		
			} 
			ISI = tspike[i][j] - tspike[i][j-1];
			AvgISI[i] =  AvgISI[i] + ISI;
			count_isi = count_isi+1;	
		}
		AvgISI[i] = AvgISI[i]/count_isi;
	}

	for(i = 0; i < N; i++)
	{
		FinalISI[i] = FinalISI[i] + AvgISI[i]; // ISI averaged over all neurons
	}

	fclose(fspiketimeid);

	int Y[N][no_cycles]; // 
	/*int **Y;
	Y = (int **)malloc(N*sizeof(float));
		for(c = 0; c < N; c++)
		{
		Y[c] = (int *)malloc(no_cycles * sizeof(int));
		}
	*/
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < no_cycles; j++)
		{
			Y[i][j] = 0;
		}
	}

	for(i = 0; i < N; i++)
	{
		for(a = 0; a < (n_steps/n_one_cycle) ; a++)  // n_steps/n_one_cycle is value of time bin
		{
			if (tspike[i][a]==0)
				continue;
			int bin_number = (int)(floor(tspike[i][a])/timewell); // timewell = n_steps/n_one_cycle 
			Y[i][bin_number] = 1;
		}
	}


	for(i = 0; i < N; i++)
	{
		for(j = 0; j < no_cycles; j++)
		{
			fprintf(yid,"%d\t",Y[i][j]);
		}
		fprintf(yid,"\n");
	}

	int sum_spikes = 0;
	for(i = 0; i < N-1; i++)
	{
		for(j = i+1; j < N;  j++)
		{
			for(l = 0; l < no_cycles; l++)
			{
				if(Y[i][l] == 1 && Y[j][l] == 1)
				{
					sum_spikes = sum_spikes + 1;
				}
			}	
			f_correlation[i][j] = (double)sum_spikes/sqrt(no_spikes[i]*no_spikes[j]);
			f_correlation[j][i] = (double)sum_spikes/sqrt(no_spikes[i]*no_spikes[j]);
			sum_spikes = 0;
		}
	}

	for(i = 0; i < N-1; i++)
	{
		for(j = i+1; j < N;  j++)
		{
			sum_total  = sum_total + f_correlation[i][j];
		}
	}
	sum_total = 2.0*sum_total/(N*(N-1));
	firing_coherence = firing_coherence + sum_total; // population coherence
	printf("Firing Coherence: %f\n", sum_total);

	fclose(yid);


	int distance[N][N]; //Array containing shortest paths between all pairs of nodes // FLOYD WARSHALL Algorithm implementation
	//int i,j,l;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			if (edgematrix[i][j] == 0)
			{
				if(i == j)
				{
					distance[i][j] = 0;
					distance[j][i] = 0;
				}
				else
				{
					distance[i][j] = 100; //100 is just a very large number. Symbolizes infinity.
					distance[j][i] = 100;
				}		
			}
			else
			{

				distance[j][i] = 1;
				distance[i][j] = 1;
			}
		}
	}

	//Calculating shortest paths
	for (l = 0; l < N; l++)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				if (distance[i][j] > (distance[i][l] + distance[l][j]))
				{
					distance[i][j] = distance[i][l] + distance[l][j];
				}
				else
					continue;
			}
		}
	}

	int max_dist = 0;
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
	//		printf("%d\t",distance[i][j]);
			if(distance[i][j] > max_dist)
				max_dist = distance[i][j];
		}
	//	printf("\n");
	}

	double fc_degree[max_dist]; // firing correlation by degree of separation.
	for(i = 0; i < max_dist; i++)
		fc_degree[i] = 0.0;

	double sum_fc_degree = 0.0;
	int count = 0;
	for(i = 1; i <= max_dist; i++)
	{
		sum_fc_degree = 0.0;
		count = 0;
		for(j = 0; j < N-1; j++)
		{
			for(l = i+1; l < N; l++)
			{
				if(distance[j][l] == i)
				{
					sum_fc_degree = sum_fc_degree + f_correlation[j][l];
					count  = count + 1;			
				}
			}
		}
		fc_degree[i-1] = sum_fc_degree/count;
	}


	for(i = 0; i < max_dist; i++)
	{
		printf("%d:  %f\n",i+1,fc_degree[i]);
		tot_degree[i] = tot_degree[i] + fc_degree[i];
		count_degree[i] = count_degree[i] + 1;
	}



	}
	//FILE *degreefcid;
	//degreefcid = fopen("FiringCoherenceDegree.txt", "w");


	firing_coherence  = firing_coherence/10.0;
	printf("Firing Coherence = %f\n", firing_coherence);
	fprintf(firingcoherenceid,"%f\t%f\n",D,firing_coherence);

	fprintf(fcdegreeid,"#D = %f\n",D);
	fprintf(fcdegreeid,"\n");
	for(i = 0; i < 12; i++)
	{
		tot_degree[i] = tot_degree[i]/count_degree[i];
		//fprintf(degreefcid, "%d\t%f\n",i+1,tot_degree[i]);
		fprintf(fcdegreeid, "%d\t%f\n",i+1,tot_degree[i]);
	}
	fprintf(fcdegreeid,"\n");

	fprintf(fisiid,"#D = %f\n",D);
	double mean = 0.0;
	double mean1 = 0.0;
	double stddev = 0.0;
	for(i = 0; i < N; i++)
	{
		FinalISI[i] = FinalISI[i]/10.0;
		mean  = mean + FinalISI[i];
		mean1 = mean1 + (FinalISI[i]*FinalISI[i]);	
		//connections[i] = connections[i]/2.0;
		// FinalISI, connections averaged over all trials. 
		fprintf(fisiid,"%d\t %f\t %d\t %f\n",i+1,connections[i],no_spikes[i],FinalISI[i]); 
	}
	fprintf(fisiid,"\n");

	mean = mean/N;
	mean1 = mean1/N;
	stddev = sqrt(mean1 - (mean*mean));
	fprintf(fmeanisiid,"#D = %f\n",D);
	fprintf(fmeanisiid,"Mean: %f\n",mean);
	fprintf(fmeanisiid,"Standard Deviation: %f\n",stddev);
	fprintf(fmeanisiid,"\n");
	//fclose(degreefcid);

	}

	fclose(fcdegreeid);
	fclose(firingcoherenceid);
	fclose(fisiid);
	fclose(fmeanisiid);

	return 0;
	}


	//End of main function

	/*Function to generate random number between 0 and max_num
	  Random number is generated using current time*/
	int randomnumbergenerator(int max_num)
	{
	    int result = 0;
	    result = (rand()%(max_num+1));
	    return result;
	}

	//Calculating the clustering coefficient.
	/*Clustering coefficient is given by:
	3*(number of triangles)/(number of triads) */
	double clusteringcoefficient(int a[][N], int m)
	{
	int i,j,l = 0;
	int count = 0; 
	double result = 0.00; //clustering coefficient to be returned
	int no_triads = 0; //Number of triads
	int no_triangles = 0; //No of triangles

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < N; j++)
		{
			if(a[i][j] == 1)
			{
				count++;
			}
		}
	}

	no_triads = count*(count-1)/2;
	no_triads = count*(count-1)/2;
	printf("Number of triads = %d\n", no_triads);
	for (i = 0; i < N-2; i++)
	{
		for(j = i+1; j < N-1; j++)
		{
			for(l = j+1; l < N; l++)
			{
				if (a[i][j] == 1 && a[j][l] == 1 && a[l][i] == 1)
				{
					no_triangles++;
				}
			}
		}
	}	
	result = (3*((float)no_triangles)/((float)no_triads));
	return result;
	}

	//Calculating Mean Shortest Path
	/*You calculate the shortest path between all pairs of nodes using the Floyd-Warshall Algorithm
	and then use the formula l(G) = (2/n.(n-1))*(sum of all unique shortest path lengths) */

	//Initializing distance matrix with Weights
	//All the weights are 1 in our case
	

	double meanshortest(int a[][N], int m) // This function is identical to the portion earlier for calculating shortest path
	{
	int i,j,l = 0;
	int distance[N][N]; //Array containing shortest paths between all pairs of nodes
	double result = 0.00;
	int sumofshortestpaths = 0; //Value used to calculate mean shortest path
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			if (a[i][j] == 0)
			{
				if(i == j)
				{
					distance[i][j] = 0;
					distance[j][i] = 0;
				}
				else
				{
					distance[i][j] = 100; //100 is just a very large number. Symbolizes infinity.
					distance[j][i] = 100;
				}		
			}
			else
			{

				distance[j][i] = 1;
				distance[i][j] = 1;
			}
		}
	}

	//Calculating shortest paths
	for (l = 0; l < N; l++)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				if (distance[i][j] > (distance[i][l] + distance[l][j]))
				{
					distance[i][j] = distance[i][l] + distance[l][j];
				}
				else
					continue;
			}
		}
	}

	//Calculating mean shortest path using the formula l(G) = (2/n.(n-1))*(sum of all unique shortest path lengths)
	for (i = 0; i < N; i++)
	{
		for (j = i+1; j < N; j++)
		{
			sumofshortestpaths = sumofshortestpaths + distance[i][j];
		}
	}
	result = 2*(float)sumofshortestpaths/(N*(N-1));
	return result;
	}	

/* ______________________________________________________________________________________________________________________________________*/
/* ______________________________________________________________________________________________________________________________________*/

	double KV( double C, double V,double m,double h,double n )
	{
		double p ;
		int sum = 0;
		p =  (((D*C) - ((GNa*xna*(m*m*m)*h*(V-VNa)) + (GK*xk*(n*n*n*n)*(V-VK))+(GL*(V-VL))))/c);
		return (p) ;
	}

	double Km( double V , double m )
	{
		double q;
		q =  ((((0.1*(25-V))/(exp((25-V)/10) - 1))*(1 - m)) - (m*(4*exp(-V/18))));
		return (q);
	}

	double Kh( double V , double h )
	{
		double r;
		r = ((0.07*exp(-V/20))*(1-h)) - ((1/(exp((30-V)/10) + 1))*h);
		return (r);
	}

	double Kn( double V , double n )
	{
		double s;
		s = (0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80));
		return (s);
	}


double Kcoupling1(double a[N], double b[][N], int m, double k[N], int index)
	{
		//a[N] is the V[j] array at a particular time
		//k[N] is the runge kutta coefficient matrix
		//b[][N] is the edgematrix
		//m is the runge-kutta step flag
		double sum_coupling = 0.0;
		int i = 0;
		for(i = 0; i<N; i++)
		{
			sum_coupling = sum_coupling + (b[index][i]*((a[i] + k[i]) - (a[index] + k[index])));
		}
		return sum_coupling;
	}




/* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */
	double gaussrandm()
	{
  		static double V1, V2, S;
  		static int phase = 0;
  		double X;

  		if(phase == 0) {
    	do {
      			double U1 = (double)rand() / RAND_MAX;
      			double U2 = (double)rand() / RAND_MAX;

      			V1 = 2 * U1 - 1;
      			V2 = 2 * U2 - 1;
      			S = V1 * V1 + V2 * V2;
      		} 
      	while(S >= 1 || S == 0);

    	X = V1 * sqrt(-2 * log(S) / S);
  		} else
    		X = V2 * sqrt(-2 * log(S) / S);

  			phase = 1 - phase;

  		return X;
	}

	double gaussrandh()
	{
  		static double V1, V2, S;
  		static int phase = 0;
  		double X;

  		if(phase == 0) {
    	do {
      			double U1 = (double)rand() / RAND_MAX;
      			double U2 = (double)rand() / RAND_MAX;

      			V1 = 2 * U1 - 1;
      			V2 = 2 * U2 - 1;
      			S = V1 * V1 + V2 * V2;
      		} 
      	while(S >= 1 || S == 0);

    	X = V1 * sqrt(-2 * log(S) / S);
  		} else
    		X = V2 * sqrt(-2 * log(S) / S);

  			phase = 1 - phase;

  		return X;
	}

	double gaussrandn()
	{
  		static double V1, V2, S;
  		static int phase = 0;
  		double X;

  		if(phase == 0) {
    	do {
      			double U1 = (double)rand() / RAND_MAX;
      			double U2 = (double)rand() / RAND_MAX;

      			V1 = 2 * U1 - 1;
      			V2 = 2 * U2 - 1;
      			S = V1 * V1 + V2 * V2;
      		} 
      	while(S >= 1 || S == 0);

    	X = V1 * sqrt(-2 * log(S) / S);
  		} else
    		X = V2 * sqrt(-2 * log(S) / S);

  			phase = 1 - phase;

  		return X;
	}

