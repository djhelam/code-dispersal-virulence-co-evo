/*
	Copyright (C) 2024  Jhelam N. Deshpande
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

//============================================================================


#include <cstdlib>			//standard C library
#include <iostream>			//standard input output function library
#include <fstream>			//file stream library
#include <numeric>		
#include <string>
#include <vector>
#include<cmath>			//need to raise numbers to powers
#include <gsl/gsl_rng.h>     //random number generatir gsl library      
#include <gsl/gsl_randist.h> //gsl random distribution library
#include <algorithm>	
#include <math.h>			 //standard math library
using namespace std;

//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------Class defining the individuals
class TInd{       
public:
	TInd();
	bool infection_state;			//stores the infectious state of an individual 0 if susceptible 1 if infected
	double dispersal_probability;	//stores the host dispersal genotype
	double resistance; 				//stores the host resistance genotype
	double virulence;				//stores the parasite virulence
};
//--------------------------------------------------------------------------------------------Constructor for class TInd
TInd::TInd(){     
	infection_state=0;
	dispersal_probability=0.0;
	resistance=0.0;
	virulence=0.0;
}

//______________________________________________________________________________________________________________________
//--------------------------------------------------------------------------------------------Class defining the patches
class TPatch      
{
public:
	TPatch();
	vector <TInd> females;      //females in a patch
	vector <TInd> newfemales;   //vector to store new disperses of new born females in a patch
	vector<int> neighbours;		//stores the neighbours of a given patch
	double measured_dispersal;	//stores measred dispersal probability in a patch
	double measured_resistance;	//stores measured resistance in a patch
	double measured_virulence;	//stores measured virulence in a patch	
	double measured_transmission_rate;	//stores measured transmission rate in a patch
	int host_extinction; // Parasite induced host extinctions before recolonaization is possible
};
//------------------------------------------------------------------------------------------Constructor for class TPatch
TPatch::TPatch(){     
	females.clear();
	newfemales.clear();
	measured_dispersal=0.0;
	measured_resistance=0.0;
	measured_virulence=0.0;
	measured_transmission_rate=0.0;
	host_extinction=0;
}

//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------------Declare global variables
//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------------------Model parameters 
const int NUMBER_OF_PARAMETERS=24;
int No;         				//Initial number of individuals per patch
double INITIAL_PREVALENCE;		//Initial prevalence in the range core
double LAMBDA;      			//Growth rate or mean female host fecundity
double ALPHA;     				//Host intraspecific competition coefficient
int BURN_IN_TIME;				//Time steps before invasion begins
int REPLICATES;					//number of replicates to run
double DISPERSAL_PROBABILITY;  	// Dispersal probability if dispersal does not evolve
double DISPERSAL_MORTALITY;  	// Cost of dispersal
double EXTINCTION_PROBABILITY;  //probability of random patch extinction
double RESISTANCE;				//set host resistance if it does not evolve
double RESISTANCE_COST;			//cost to resistance
double STANDING_GENETIC_VARIATION_HOST_DISPERSAL;		//standing genetic variation for host traits
double STANDING_GENETIC_VARIATION_HOST_RESISTANCE;		//standing genetic variation for host traits
double STANDING_GENETIC_VARIATION_PARASITE;	//standing genetic variation for parasite traits
double MUTATION_RATE_HOST_DISPERSAL;		//mutation rates of host traits
double MUTATION_RATE_HOST_RESISTANCE;		//mutation rates of host traits
double MUTATION_RATE_PARASITE;	//mutation rates of parasite traits
double MUTATION_EFFECT_SD_HOST;	//stanfdard deviation of mutation effects host
double MUTATION_EFFECT_SD_PARASITE;//standard deviation of mutation effects parasite
double VIRULENCE;			    //parasite virulence if no evolution
double SEARCH_EFFICIENCY;		//search efficiency parameter for transmission
double SLOPE_VIRULENCE;			//slope of the line that relates virulence to transmission
double MAXIMUM_TRANSMISSION;	//maximum transmission rate
bool RANDOM_NETWORK;			//stores whether a given network is fixed or random

//______________________________________________________________________________________________________________________
//--------------------------------------------------------------------------------------------------Landscape properties
const int NUMBER_OF_PATCHES=100;	//Defining number of patches
TPatch world[NUMBER_OF_PATCHES];	//Creating patches

//______________________________________________________________________________________________________________________
const gsl_rng *gBaseRand;	//Seed for random number generator

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);

	srand(randSeed);
	unsigned long r = rand();
	gsl_rng_set(gBaseRand, r);
}

//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------------Simplifications

//-------------------------------------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

//---------------------------------------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

//-----------------------------------------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}

double logit(double p)
{
	return log(p/(1.0-p));
}

double inverse_logit(double x)
{
	return 1.0/(1.0+exp(-x));
}


const int RS = 100;                 // random seed

//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------Reading and setting parameters
void input_parameters()  
{
	string para[NUMBER_OF_PARAMETERS];							//stores model parameter values
	string line;	//string stores each line of the parameter input file
	ifstream myfile ("input.txt");				//read the parmeter input file
	int count=0;
	if (myfile.is_open())
	{
		while ( getline (myfile,line))
		{
			if(count%2==1)
			{
				para[count/2]=line;				//store only numeric values from the input file
			}
			count++;

		}
		myfile.close();
	}
	else cout << "Unable to open file";
	No = (int) std::atof(para[0].c_str());				//sets initial population size per patch
	INITIAL_PREVALENCE=std::atof(para[1].c_str());		//set initial prevalence in range core
	LAMBDA=std::atof(para[2].c_str());					//sets mean fecundity of the female
	ALPHA=std::atof(para[3].c_str());			//sets intra-specific competition coefficient of the Beverton-Holt model
	BURN_IN_TIME= (int) std::atof(para[4].c_str());//sets the number of time steps before the beginning of invasion
	REPLICATES=(int) std::atof(para[5].c_str());	    //sets number of replicate simulations that are run
	DISPERSAL_PROBABILITY=std::atof(para[6].c_str());	//sets dispersal probability if it is not genetically encoded
	DISPERSAL_MORTALITY=std::atof(para[7].c_str());		//sets dispersal costs
	EXTINCTION_PROBABILITY=std::atof(para[8].c_str());	//sets probability of random patch extinction
	RESISTANCE=std::atof(para[9].c_str());				//sets resistance if it does not evolve
	RESISTANCE_COST=std::atof(para[10].c_str());				//sets cost to resistance
	STANDING_GENETIC_VARIATION_HOST_DISPERSAL=std::atof(para[11].c_str());//sets standing genetic variation dispersal
	STANDING_GENETIC_VARIATION_HOST_RESISTANCE=std::atof(para[12].c_str());//sets standing genetic variation resistance
	STANDING_GENETIC_VARIATION_PARASITE=std::atof(para[13].c_str());//sets standing genetic variation for virulence
	MUTATION_RATE_HOST_DISPERSAL=std::atof(para[14].c_str());	//sets standing mutation rate for host dispersal
	MUTATION_RATE_HOST_RESISTANCE=std::atof(para[15].c_str());	//sets standing mutation rate for host resistance
	MUTATION_RATE_PARASITE=std::atof(para[16].c_str());//sets mutation rate for parasite traits
	MUTATION_EFFECT_SD_HOST =std::atof(para[17].c_str());	//sets standing mutation size for host traits
	MUTATION_EFFECT_SD_PARASITE=std::atof(para[18].c_str());//sets mutation size for parasite traits
	VIRULENCE=std::atof(para[19].c_str());	//sets disease transmission rate
	SEARCH_EFFICIENCY=std::atof(para[20].c_str()); //sets search efficiency of transmission
	SLOPE_VIRULENCE=std::atof(para[21].c_str()); //sets how quickly transmission increases with virulence
	MAXIMUM_TRANSMISSION=std::atof(para[22].c_str());	//sets the maximum transmission rate possible
	RANDOM_NETWORK=(bool) std::atof(para[23].c_str());

}

//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------------Initialise the landscape
void initialise_landscape()
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)				//loop through all patches
	{

		world[x].females.clear();	    //clear females in the patch	
		world[x].newfemales.clear();	//clear new females in patch
		world[x].neighbours.clear();	//clear neighbours of patch
		world[x].measured_dispersal=0.0;
		world[x].measured_virulence=0.0;
		world[x].measured_resistance=0.0;
		world[x].measured_transmission_rate=0.0;
	}
	for(int x=0;x<NUMBER_OF_PATCHES;x++) 		//loop through
	{
		for(int n=0;n<No;n++)	//create No new individuals 
		{
			TInd newind;	//create a new individual
				//if there is standing genetic variation in host  dispersal trait
			if(STANDING_GENETIC_VARIATION_HOST_DISPERSAL>0)
					//initialise dispersal probability with standing genetic variation
				newind.dispersal_probability=ran()*STANDING_GENETIC_VARIATION_HOST_DISPERSAL;
				//if there is no standing genetic variation for dispersal trait
			else 
					//initialise with dispersal parameter
				newind.dispersal_probability=DISPERSAL_PROBABILITY;
				//if there is standing genetic variation for resistance
			if(STANDING_GENETIC_VARIATION_HOST_RESISTANCE>0)
					//initialise resistance trait with standing genetic variation
				newind.resistance=ran()*STANDING_GENETIC_VARIATION_HOST_RESISTANCE;	
				//if there is no standing genetic variation for host resistance
			else
					//initialise with external model parameter
				newind.resistance=RESISTANCE;
				//individuals are infected with a probability INITIAL_PREVALENCE
			if(ran()<INITIAL_PREVALENCE)
			{
					newind.infection_state=1;	//set infection state to infected
					//if there is standing genetic variation for parasite virulence
					if(STANDING_GENETIC_VARIATION_PARASITE>0)
						//initialise parasite virulence with standing genetic variation
						newind.virulence=ran()*STANDING_GENETIC_VARIATION_PARASITE;	
					else newind.virulence=VIRULENCE;
				}
				else
				{
				newind.infection_state=0;	//set infection state to Susceptible
				newind.virulence=0;			//set virulence to 0
			}
			world[x].females.push_back(newind);
		}
	}
}

//______________________________________________________________________________________________________________________
//---------------------------------------------------------------------------Input adjacency matrix defining a landscape
void input_adjacency_matrix(int r) 	//r is the replicate number
{
	string filename; 		//stores name of input file
	filename="../matrices/adjacency_matrix_"+to_string(r)+".txt";	//read new adjacency matrix for every replicate
	//filename="adjacency_matrix_"+to_string(r)+".txt";
	ifstream adj;
	adj.open(filename.c_str()); //open adjacency matrix input file
	int adjacency_matrix[NUMBER_OF_PATCHES][NUMBER_OF_PATCHES];	//stores the adjacency matrix
	if(adj.is_open())	//if the file is open
	{
		for(int x=0;x<NUMBER_OF_PATCHES;x++)	
		{
			for(int y=0;y<NUMBER_OF_PATCHES;y++)
			{
				adj>>adjacency_matrix[x][y];	//store the adjacency matrix
			}
		}
		adj.close();
	}

	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all patches x
	{
		for(int y=0;y<NUMBER_OF_PATCHES;y++)	//go through all elements in the row x of the adjacency matrix
		{
			if(adjacency_matrix[x][y]==1)	//if patch y is a neightbor of patch x
			{
				world[x].neighbours.push_back(y);	//append patch y to the vector of patch x's neighbours
			}
		}
	}
}



//______________________________________________________________________________________________________________________
//--------------------------------------------------------------------Deciding coordinates of new patch while dispersing
int decide_patch(int x) 
{ 
	int new_x;
	if(world[x].neighbours.size()>0)
	{
		int pos=floor(ran()*world[x].neighbours.size());	//choose index of the new location of the disperser
		new_x=world[x].neighbours.at(pos);
	}
	else new_x=-1000;
	return new_x;	//return new patch of disperser
}

//______________________________________________________________________________________________________________________
//---------------------------------------------------------------------------------------------------mutation procedures
double mutate_dispersal(double d)
{
	//add mutation to dispersal trait with a given probability
	double mutation_effect=gauss(MUTATION_EFFECT_SD_HOST);
	if(ran()<MUTATION_RATE_HOST_DISPERSAL)
		d=inverse_logit(logit(d)+mutation_effect); //mutation drawn from a normal distribution
	return (d);
}

double mutate_resistance(double r)
{
	//add mutation to resistance trait with a given probability
	if(ran()<MUTATION_RATE_HOST_RESISTANCE)
		r=r+gauss(MUTATION_EFFECT_SD_HOST); //mutation drawn from a normal distribution
	if(r<0)
		r=0;
	if(r>1)
		r=1;
	return (r);
}

double mutate_virulence(double v)
{
	//add mutation to parasite virulence with gven probability
	double mutation_effect=gauss(MUTATION_EFFECT_SD_PARASITE);;
	if(ran()<MUTATION_RATE_PARASITE)
		v=inverse_logit(logit(v)+mutation_effect); //mutation drawn from a normal distribution
	return v;
}

//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------------------equations
//______________________________________________________________________________________________________________________
//--------------------------------------------------------------------------------------Density regulation Beverton-Holt
double density_regulation(double population_size)
{
	return (1/(1+(ALPHA*population_size)));
}
//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------------------Virulence
double virulence_calculation(double virulence)
{
	return 1-virulence;
}
//------------------------------------------------------------------------------------------------------------Resistance
double resistance_cost_calculation(double resistance)
{
	return 1-resistance*RESISTANCE_COST;
}
//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------Calculate transmission probability
double transmission_probability(double S,double beta)
{
	return 1-exp(-SEARCH_EFFICIENCY*beta/(1.0+SEARCH_EFFICIENCY*S));
}
double virulence_transmission_tradeoff(double virulence)
{
	return  MAXIMUM_TRANSMISSION*pow(virulence,SLOPE_VIRULENCE);	//virulence transmission trade-off function
	//
	//return SLOPE_VIRULENCE*virulence;
}

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------------------output landscape
void output_landscape(ofstream& op,int r,int t)
{
	//go through all the patches
	for(int x=0;x<NUMBER_OF_PATCHES;x++)
	{
		int count_susceptible=0;	//stores number of susceptible
		int count_infected=0;
		for(int f=0;f<world[x].females.size();f++)	//goes through all individuals
		{
			if(world[x].females.at(f).infection_state==0)	//count susceptible
				count_susceptible++;
			else if(world[x].females.at(f).infection_state==1)	//count infected
				count_infected++;
		}
		op<<r<<" "<<t<<" "<<x<<" "<<world[x].females.size()<<" "<<count_susceptible<<" "<<count_infected<<" "<<world[x].measured_dispersal<<" "<<world[x].measured_resistance<<" "<<world[x].measured_virulence<<" "<<world[x].measured_transmission_rate<<" "<<world[x].host_extinction<< endl;	//output patch properties
		
	}
}
//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------------------output genotypes
void output_genotypes(ofstream& op1, int t
	, int r)
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)
	{
		for(int f=0;f<world[x].females.size();f++)
		{
			//if(world[x].females.at(f).infection_state==1)
				op1<<r<<" "<<t<<" "<<x<<" "<<world[x].females.at(f).infection_state<<" "<<world[x].females.at(f).dispersal_probability<<" "<<world[x].females.at(f).virulence<<" "<<world[x].females.at(f).resistance<<endl;
		}
	}
}


//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------life cycle procedures

//______________________________________________________________________________________________________________________
//---------------------------------------------------------------------------------------------------dispersal procedure
void disperse()
{
	//clear the newfemales vector
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all patches
	{
		world[x].newfemales.clear();	//clear newfemales vector so it can store dispersers
	}
	//dispersal
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all patches
	{
		int count_dispersers=0;
		int population_size=world[x].females.size();
		for(int f=0;f<world[x].females.size();f++)
		{
			//disperse with the genetically encoded dispersal probability
			if(ran()<world[x].females.at(f).dispersal_probability)
			{
				//dispersal mortalilty
				if(ran()>DISPERSAL_MORTALITY)	//if the individual does not die while dispersing
				{
					int new_patch=decide_patch(x);//randomly choose one of 8 nearest neighbours to disperse
					//store this individuals in the newfemales vector of its target patch
					if(new_patch!=-1000)	//if the patch actually has neighbours
						world[new_patch].newfemales.push_back(world[x].females.at(f)); //disperse to one of the neighbours
				}
				world[x].females.erase(world[x].females.begin()+f);	//remove this female from its old patch
				//if there is no patch surrounding the patch the disperser dies
				//the next female is now at the position in the females vector where the dispersed/dead female was
				f--;
				count_dispersers++;
			}
		}
		if(population_size>0){
			world[x].measured_dispersal=double(count_dispersers)/double(population_size);
		}
		else
			world[x].measured_dispersal=0;
	}

	//add the dispersers to their target patches
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all the patches
	{
		if(world[x].newfemales.size()>0)	 //if there are dispersers that are arriving in the patch
		{
			for(int f=0;f<world[x].newfemales.size();f++)
			{
				world[x].females.push_back(world[x].newfemales.at(f));//add these dispersers to the new patch
			}
		}
		world[x].newfemales.clear();	//clear the newfemales vector
	}
}

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------------reproduction procedure
void reproduce()
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)			//loop through all the patches
	{
		world[x].newfemales.clear();	//clear all individuals in the patch
		double measure_virulence=0.0;	//track the measured virulemce in the patch
		int count_infected=0;
		int population_size=world[x].females.size();
		for(int f=0;f<world[x].females.size();f++)
		{
			//mean number of offspring calculated from the Beverton-Holt model
			//its fecundity is reduced according to the virulence of the parasite it bears	
			//fecundity is reduced according to cost of resistance		
			double mean_offspring=LAMBDA*resistance_cost_calculation(world[x].females.at(f).resistance)*virulence_calculation(world[x].females.at(f).virulence)*density_regulation(double(population_size));
			if(world[x].females.at(f).infection_state==1)//if the individual is infected
			{
				measure_virulence=measure_virulence+world[x].females.at(f).virulence;	//measure virulence
				count_infected++;
			}
			int number_of_offspring=poisson(mean_offspring);//number of offspring drawn from a Poisson distribution
			for(int b=0;b<number_of_offspring;b++)
			{
				TInd newind;	//create new individual
				//inherit the mother's allele for dispersal probability
				newind.dispersal_probability=mutate_dispersal(world[x].females.at(f).dispersal_probability);
				//inherit the mother's allele for resistance
				newind.resistance=mutate_resistance(world[x].females.at(f).resistance);
				newind.infection_state=0; //all individuals are born susceptible
				newind.virulence=0;	//virulence acting on susceptible individuals is 0
				world[x].newfemales.push_back(newind);
			}

		}

		if(world[x].newfemales.size()==0 && world[x].females.size()>0)
		{
			world[x].host_extinction=1;

		}
		else world[x].host_extinction=0;

		if(count_infected!=0)
			world[x].measured_virulence=measure_virulence/double(count_infected);	//store measured virulence as a patch property
		else world[x].measured_virulence=0;

	}

} 
//______________________________________________________________________________________________________________________
//----------------------------------------------------------------------------------------disease transmission procedure
void transmission()
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)		//go through all the patches
	{
		vector<double> all_virulence;	//stores the genotypic virulence value of the parasite infecting each infected
		for(int f=0;f<world[x].females.size();f++)	//go through all females in the parent generation
		{
			if(world[x].females.at(f).infection_state==1)//if they are infected
			{
				all_virulence.push_back(world[x].females.at(f).virulence); //store the genotypic value of the virulence infecting that infected individual 
			}

		}
		//if there are infected individuals in the parental generation
		int count_infected=0;
		int count_resistant=0;
		if(all_virulence.size()>0)
		{
			//go through all the newborn susceptible offspring
			for(int nf=0;nf<world[x].newfemales.size();nf++)
			{
				vector<double> contact_virulence; //this vector stores the virulence of a given parasite strain in contact with a susceptible
				for(int av=0;av<all_virulence.size();av++)
				{
					//check whether this susceptible is in contact with a given parasite strain
					if(ran()<transmission_probability(world[x].newfemales.size(),virulence_transmission_tradeoff(all_virulence.at(av))))
					{
						contact_virulence.push_back(all_virulence.at(av));	//store virulence of parasite strain in contact with a given newborn susceptible individual
					}		
				}
				if(contact_virulence.size()>0)	//if the given newborn is in contact with any parasite strain
				{
					if(ran()<1-world[x].newfemales.at(nf).resistance)	//if the the newborn female is not resistant
					{
						world[x].newfemales.at(nf).infection_state=1;	//the newborn is infected
						//we now determine the starin it is infected by by drawibg one strain from those the newborn is in contact with
						int pos=floor(ran()*contact_virulence.size());	//draw a strain at random to infect the given susceptible
						world[x].newfemales.at(nf).virulence=mutate_virulence(contact_virulence.at(pos));	//assign parasite strain, it may mutate while replicating within the host
						count_infected++;
					}
					else
					{
						count_resistant++;
					}
				}
			}
			if(count_resistant+count_infected!=0)
			{
			//stores fraction of resistant as a patch property
				world[x].measured_resistance=double(count_resistant)/double(count_resistant+count_infected);
			}
			else
			{
				world[x].measured_resistance=0;
			}

		}
	}
}




//______________________________________________________________________________________________________________________
//-------------------------------------------------------------------------------------------------------death procedure
void death()
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//loop through all patches
	{

		world[x].females.clear();	//parent generation dies
		world[x].females=world[x].newfemales;	//offspring replace their parents
	}
}

//______________________________________________________________________________________________________________________
//------------------------------------------------------------------------------------------add random patch extinctions
void patch_extinction()	//random patch extinction
{
	for(int x=0;x<NUMBER_OF_PATCHES;x++)	//go through all the patches
	{
		if(ran()<EXTINCTION_PROBABILITY)	//the patch is cleared wih a probability EXTINCTION_PROB
		{
			world[x].females.clear();
		}
		
	}
}

int main()
{
	ofstream op;
	op.open("metapopulation.txt"); 	//output the population size, prevalence.etc
	op <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"N S I measured_dispersal measured_resistance measured_virulence measured_transmission host_extinction";
	op<<endl;
	ofstream op1;
	op1.open("genotypes.txt");	
	op1 <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"infection_state dispersal_probability virulence resistance";
	op1<<endl;
	specify_rng(RS);	//set the seed of the the random number generator
	input_parameters();	//input and set model parameters

	for(int r=0;r<REPLICATES;r++)	//go through different replicates
	{
		initialise_landscape();		//initialise the landscape without specifying connectivity
		if(RANDOM_NETWORK==0)		//specify connectivity of landscape for a fixed graph
			input_adjacency_matrix(0);
		if(RANDOM_NETWORK==1)		//specify connectivity of a landscape for a random graph
			input_adjacency_matrix(r);	
		int t=0;	//set time step to 0
		do
		{
  if(t>BURN_IN_TIME-100)	//output landscape only in the last 100 time steps
				output_landscape(op,r,t);	//output the landscape
			//host life cycle begins
			disperse();	//natal dispersal
			reproduce();	//reproduction
			transmission();	//parasite transmission to offspring generation
			death();	//death of parent generation
			patch_extinction();	 //random patch extinction
			t++;	//increase time step counter
			if(t==BURN_IN_TIME)
				output_genotypes(op1,t,r);
		}
		while(t<BURN_IN_TIME);

	}
	op.close();	//close output file
	op1.close();
	return 0;
}