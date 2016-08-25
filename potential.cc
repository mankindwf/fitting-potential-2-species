#include <fstream>       
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <infer/state.h>
#include <infer/objfunction.h>
#include <infer/pso.h>
#define LBOX 3.61*3

//Para La interaccion HCl

//Prueba

using namespace std;
using namespace pso;



typedef std::vector<double> config;//config : vector de doubles (Lista de to pociciones) 3*numero_atomos
typedef std::vector<std::string> species_config; //species_config : vector de strings (Lista de species)

//Cada cuanto mostrar 

int n_show = 10;

// Sobrecargamos el constructor de la clase Minimizer
class MyMinimizer: public Minimizer
{
 public:
  MyMinimizer(int nparams, int nparticles): Minimizer(nparams, nparticles) { }

  void OnIteration(int n, double objval, const State & current) const
  {
  if(n%n_show==0){
  //Lo que muestra el código al ejecutarce. ( Numero de Paso=n, Valor de la función Error= objval,Parametros=current)
  cout << n/n_show << " Error: " << objval << "  Params: " << current << endl;
  }}
};


//Creamos la clase Potential, es un conjunto de potenciales en este caso para 3 interacciones (H-H, Cl-Cl, H-Cl)

class Potential{
public:
	//El numero de parametros sera 6.
	virtual double operator()(const config &conf, const species_config &species, const State &param)=0;
	int N;
};

// En este caso todos interactuan con LJ
// H-H
double PairPotentialHH(double r, double epsilon, double sigma){
	return 4.0*epsilon*(pow(sigma/r,12) - pow(sigma/r,6) );}

// Cl-Cl
double PairPotentialHCl(double r, double epsilon, double sigma){
	return 4.0*epsilon*(pow(sigma/r,12) - pow(sigma/r,6) );}

// H-Cl
double PairPotentialClCl(double r, double epsilon, double sigma){
	return 4.0*epsilon*(pow(sigma/r,12) - pow(sigma/r,6) );}


double Energia(const config &conf, const species_config &species, const State &param){
		double LHALF = 0.5*LBOX;
		int natoms = conf.size()/3;
		double U = 0.0e0;
		for(int i=0;i<natoms-1;++i){for(int j=i+1;j<natoms;++j){
    		double dr[3], r2=0.0e0;
    		for (int q=0;q<3;++q) dr[q] = conf[j*3+q]- conf[i*3+q];
    		for (int q=0;q<3;++q){
     		double * dd = dr+q;
     		if (*dd >= LHALF) (*dd) -= LBOX;
     		else if (*dd < -LHALF) (*dd) += LBOX;
     		r2 += (dr[q]*dr[q]);
    		}//Agregar el if para usar diferente pairpotential en cada combinacion de especie
		if(species[i]==species[j]){
			//Potencial Para H-H 
			if(species[i]=="H"){U += PairPotentialHH(sqrt(r2), param[0], param[1]);}
			//Potencial para Cl-Cl
			else{		    U += PairPotentialClCl(sqrt(r2), param[2], param[3]);}}
			//Potencial H-Cl
		else{		    
				    U += PairPotentialHCl(sqrt(r2), param[4], param[5]);}
    		}}
 	return U;}

// 
class LJ: public Potential{
public:
	LJ() { N = 6; }
	double operator()(const config &conf, const species_config &species, const State &param)
	{return Energia(conf, species, param);}
};

void LeerPos(char *data_pos, vector <config> &all_atoms, vector <species_config> &all_species){
		ifstream dataxyz;
		dataxyz.open(data_pos);
		int i=0;
		while( 1 ){
			int num_atoms;
			dataxyz >> num_atoms;
			if( dataxyz.eof() ){break;}
			config atoms(3*num_atoms);
			species_config species(num_atoms);
			//Lee Formato xyz  (Specie x y z )
		for(int k=0; k<num_atoms; k++){dataxyz >> species [k] >> atoms[3*k] >> atoms[3*k+1] >> atoms[3*k+2];
			}
			all_species.push_back(species);
			all_atoms.push_back(atoms);
			i++;}}


class Error: public ObjectiveFunction{
public:
	Error(char *data_pos, char *data_ene, Potential &Pot): ObjectiveFunction(Pot.N), potential(Pot)
               {
		//Lectura de posiciones y especies
		LeerPos(data_pos, all_atoms, all_species);
		//Lectura de Energias
		ifstream dataene;
		dataene.open(data_ene);
		double each_energy;
		while(!dataene.eof()){  dataene >> each_energy;
					energies.push_back(each_energy);}
		}

	double operator()(const State &param)const{	
		double total=0.0;		
		int num_config = all_atoms.size();
		//La funcion Error, es la diferencia de la energia potencial con las energias para cada configuracion
		for(int i=0; i<num_config; i++){
				total += pow(potential(all_atoms[i], all_species[i], param) - energies[i],2);}
		return total;
	}

private:
	vector <config> all_atoms;
	vector <species_config> all_species;
	config energies;
        Potential & potential;
	
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//					 	INICIO	 						////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv){
MPI_Init(&argc, &argv);


vector <config> all_atoms;
vector <species_config> all_species;
LeerPos("atoms.xyz", all_atoms, all_species);
int num_config = all_atoms.size();

LJ Pot;
State param({1,3.4,0.5,4,0.7,4.5});//Todos los cores usan la misma semilla.
for(int i=0; i<num_config; i++){cout << Pot(all_atoms[i], all_species[i], param) << endl;}

Error MinLJ("atoms.xyz","energias_MultiEspecie.dat",Pot);
State s({1,3.4,0.5,4,1,4.3});
cout << MinLJ(s) << endl;

MyMinimizer min(6, 100);
min.Minimize(MinLJ, s, 1e-20);

MPI_Finalize();
return 0;}
