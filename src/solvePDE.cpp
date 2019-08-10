//Importing compiler libraries
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <math.h>
#include <chrono>
#include <fstream>
//Importing external libraries
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
#include <iomanip>
#include <nlopt.hpp>


//Parameters set at compilation time & compiler macros
#define MESH_SIZE 				10
#define maxC					0.4
#define power					2
#define PLATE_WIDTH       			0.005
#define SINK_WIDTH 				0.002
#define Ts 					293
#define Km 					65.0
#define Kp					0.2
#define TEMP_SIZE				(MESH_SIZE*MESH_SIZE)
#define MESH_SIZE2				((2*MESH_SIZE)-1)
#define MESH_WIDTH 				(PLATE_WIDTH/(MESH_SIZE-1))
#define Origq             			20000000
#define q                 			(-(Origq*MESH_WIDTH*MESH_WIDTH))
#define metal              			((2*MESH_SIZE)/5) // Number of Temperature points in the heat sink (we set their link as 100% metal)



using namespace Eigen;
using namespace std;

void computeK(const double* V,double K[MESH_SIZE2][MESH_SIZE]);
void generateSystem(double K[MESH_SIZE2][MESH_SIZE], SparseMatrix<double, ColMajor> &M,VectorXd &bb);
void dump2matlab(string filename, VectorXd sol, const double* v);
void loadV(string filename, double K[MESH_SIZE2*MESH_SIZE]);

int main(int argc, char *argv[])
{	
if (argc > 1) {
	string file = argv[1];
	file = file + ".txt";
	double vv[MESH_SIZE2*MESH_SIZE] = {0};
	loadV(file,vv);
	double KMesh [MESH_SIZE2][MESH_SIZE] = {0};
	computeK(vv,KMesh);
	//Setting up the system of equations
	SparseMatrix<double, ColMajor> M(TEMP_SIZE,TEMP_SIZE);
	VectorXd bb(TEMP_SIZE), Temp(TEMP_SIZE);
 	generateSystem(KMesh,M,bb);

	//Building and calling the solver		
	SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver;
	solver.analyzePattern(M);
	solver.compute(M);
	Temp = solver.solve(bb);
	file = argv[1];
	file = file + ".m";
	dump2matlab(file,Temp,vv);

/*cout << "K = [";
	for(int i = 0; i < MESH_SIZE*MESH_SIZE2; ++i) {
		if(i%MESH_SIZE != 0 ){cout << ", ";}
    cout <<  KMesh[i/MESH_SIZE][i%MESH_SIZE];
		if(i%MESH_SIZE == MESH_SIZE -1){cout << endl;}
  }
  cout << "];" << endl;*/
}else{
	cout << "Please provide a file to read." << endl;
}

}

void loadV(string filename,double K[MESH_SIZE2*MESH_SIZE]){
	int line = 0;
	std::string value;
  ifstream myfile (filename);
  if (myfile.is_open())
    {
        while (! myfile.eof() )
        {
            getline (myfile,value);
            K[line] = atof(value.c_str());
						line++;
        }
        myfile.close();
    }
    else cout << "Unable to open file";

}

//Takes the V matrix and builds the corresponding heat transfer coefficient matrix
void computeK(const double* V,double K[MESH_SIZE2][MESH_SIZE]){
 for (int i=0; i<MESH_SIZE2-1; i=i+2){
		for (int j=0; j<MESH_SIZE-1; j++){
        K[i][j]= Kp +  (pow(V[i*MESH_SIZE+j],power)) *(Km-Kp);
		}
	}
	for (int j=0; j<MESH_SIZE-1; j++){
    K[MESH_SIZE2-1][j]= (Kp +  (pow(V[(MESH_SIZE2-1)*MESH_SIZE+j],power)) *(Km-Kp)) / 2.0;
	}
	for (int i=1; i<MESH_SIZE2; i=i+2){
		for (int j=0; j<MESH_SIZE-1; j++){
        K[i][j]= Kp +  (pow(V[i*MESH_SIZE+j],power)) *(Km-Kp);
		}
		K[i][MESH_SIZE-1]= (Kp +  (pow(V[i*MESH_SIZE+MESH_SIZE-1],power)) *(Km-Kp)) / 2.0;
	}
}


//Builds the finite differences scheme system from the data.
void generateSystem(double K[MESH_SIZE2][MESH_SIZE], SparseMatrix<double, ColMajor> &M, VectorXd &bb){
		double qvalue =  q;
		int t_index;
		//Top border
		M.insert(0,0) =  -(K[0][0]+K[1][0]);//corner
		M.insert(0,1) =  K[0][0];//right
		M.insert(0,MESH_SIZE) =  K[1][0];//down
		bb(0) = qvalue;

		for (int j = 1; j < MESH_SIZE-1; j++){
						M.insert(j,j) =  -(K[0][j-1]+K[0][j]+K[1][j]);
						M.insert(j,j-1) =  K[0][j-1];//left
						M.insert(j,j+1) =  K[0][j];//right
						M.insert(j,j+MESH_SIZE) =  K[1][j];//down
						bb(j) = qvalue;
		}
		M.insert(MESH_SIZE-1,MESH_SIZE-1) =  1.0;//sink corner
		bb(MESH_SIZE-1) = Ts;
		//equation for the rows with heatsink points	
		for(int i = 1; i < metal; i++){
				//Left Border
						t_index = (i*MESH_SIZE);
						M.insert(t_index,t_index) =  -(K[(2*i)+1][0]+K[(2*i)-1][0]+K[2*i][0]);
						M.insert(t_index,t_index-MESH_SIZE) =  K[(2*i)-1][0];//up
						M.insert(t_index,t_index+1) =  K[2*i][0];//right
						M.insert(t_index,t_index+MESH_SIZE) =  K[(2*i)+1][0];//down
						bb(t_index) = qvalue;
				//Interior points
				for (int j = 1; j < MESH_SIZE-1; j++){
						t_index = (i*MESH_SIZE)+j;
						M.insert(t_index,t_index) =  -(K[(2*i)+1][j]+K[(2*i)-1][j]+K[2*i][j]+K[2*i][j-1]);
						M.insert(t_index,t_index-MESH_SIZE) =  K[(2*i)-1][j];//up
						M.insert(t_index,t_index-1) =  K[2*i][j-1];//left
						M.insert(t_index,t_index+1) =  K[2*i][j];//right
						M.insert(t_index,t_index+MESH_SIZE) =  K[(2*i)+1][j];//down
						bb(t_index) = qvalue;
				}
				//heat sink points
				M.insert((i*MESH_SIZE)+MESH_SIZE-1,(i*MESH_SIZE)+MESH_SIZE-1) =  1.0;
				bb((i*MESH_SIZE)+MESH_SIZE-1) = Ts;
		}
		//equation for the rows without heatsink points	
		for(int i = metal; i < MESH_SIZE -1; i++){
				//Left Border
						t_index = (i*MESH_SIZE);
						M.insert(t_index,t_index) =  -(K[(2*i)+1][0]+K[(2*i)-1][0]+K[2*i][0]);
						M.insert(t_index,t_index-MESH_SIZE) =  K[(2*i)-1][0];//up
						M.insert(t_index,t_index+1) =  K[2*i][0];//right
						M.insert(t_index,t_index+MESH_SIZE) =  K[(2*i)+1][0];//down
						bb(t_index) = qvalue;
				//Interior points
				for (int j = 1; j < MESH_SIZE-1; j++){
						t_index = (i*MESH_SIZE)+j;
						M.insert(t_index,t_index) =  -(K[(2*i)+1][j]+K[(2*i)-1][j]+K[2*i][j]+K[2*i][j-1]);
						M.insert(t_index,t_index-MESH_SIZE) =  K[(2*i)-1][j];//up
						M.insert(t_index,t_index-1) =  K[2*i][j-1];//left
						M.insert(t_index,t_index+1) =  K[2*i][j];//right
						M.insert(t_index,t_index+MESH_SIZE) =  K[(2*i)+1][j];//down
						bb(t_index) = qvalue;
				}
				//right border points
						t_index = (i*MESH_SIZE)+MESH_SIZE-1;
						M.insert(t_index,t_index) =  -(K[(2*i)+1][MESH_SIZE-1]+K[(2*i)-1][MESH_SIZE-1]+K[2*i][MESH_SIZE-2]);
						M.insert(t_index,t_index-MESH_SIZE) =  K[(2*i)-1][MESH_SIZE-1];//up
						M.insert(t_index,t_index-1) =  K[2*i][MESH_SIZE-2];//left
						M.insert(t_index,t_index+MESH_SIZE) =  K[(2*i)+1][MESH_SIZE-1];//down
						bb(t_index) = qvalue;
		}
		//Bottom border

		int i = MESH_SIZE-1;
		//j = 0
		t_index = (i*MESH_SIZE);
		M.insert(t_index,t_index) =  -(K[(2*i)-1][0]+K[2*i][0]);
		M.insert(t_index,t_index-MESH_SIZE) =  K[(2*i)-1][0];//up
		M.insert(t_index,t_index+1) =  K[2*i][0];//right
		bb(t_index) = qvalue;
		for (int j = 1; j < MESH_SIZE-1; j++){
						t_index = (i*MESH_SIZE)+j;
						M.insert(t_index,t_index) =  -(K[(2*i)-1][j]+K[2*i][j]+K[2*i][j-1]);
						M.insert(t_index,t_index-MESH_SIZE) =  K[(2*i)-1][j];//up
						M.insert(t_index,t_index-1) =  K[2*i][j-1];//left
						M.insert(t_index,t_index+1) =  K[2*i][j];//right
						bb(t_index) = qvalue;
		}
		int j = MESH_SIZE-1;
		t_index = (i*MESH_SIZE)+j;
		M.insert(t_index,t_index) =  -(K[(2*i)-1][j]+K[2*i][j-1]);//corner
		M.insert(t_index,t_index-MESH_SIZE) =  K[(2*i)-1][j];//up
		M.insert(t_index,t_index-1) =  K[2*i][j-1];//left

		bb(t_index) = qvalue;
}




//Generates a matlab file with a given name able to plot a hetmap for the given temperature vector
void dump2matlab(string filename, VectorXd sol,const double* v){
	ofstream temps;
  temps.open (filename);
	temps << "M = [";
	for (int i=MESH_SIZE-1; i >= 0; i--){
		for (int j=MESH_SIZE-1; j>=0; j--){
        temps << std::setprecision(5) << sol.coeff(i*MESH_SIZE + j) << " ";
		}
		for (int j=0; j<MESH_SIZE; j++){
        temps << std::setprecision(5) << sol.coeff(i*MESH_SIZE + j) << " ";
		}
		temps <<  endl;
	}
	for (int i=0; i<MESH_SIZE; i++){
		for (int j=MESH_SIZE-1; j>=0; j--){
        temps << std::setprecision(5) << sol.coeff(i*MESH_SIZE + j) << " ";
		}
		for (int j=0; j<MESH_SIZE; j++){
        temps << std::setprecision(5) << sol.coeff(i*MESH_SIZE + j) << " ";
		}
		temps <<  endl;
	}
	temps << "];" << endl;
	temps << "mean2(M)" << endl;
	temps << "figure;" << endl;
	temps << "heatmap(M);" << endl;
	
	temps << "v = [";
	for (int i=MESH_SIZE2-1; i >= 0; i--){
		if(i%2 == 0){
			for (int j=MESH_SIZE-2; j>=0; j--){
        temps << v[i*MESH_SIZE + j] << " ";
			}
			temps << "0 ";
			for (int j=0; j<MESH_SIZE; j++){
        temps << v[i*MESH_SIZE + j] << " ";
			}
			temps <<  endl;
		}else{
			for (int j=MESH_SIZE-1; j>=0; j--){
        temps << v[i*MESH_SIZE + j] << " ";
			}
			for (int j=0; j<MESH_SIZE; j++){
        temps << v[i*MESH_SIZE + j] << " ";
			}
			temps <<  endl;
		}		
	}
	for (int j=MESH_SIZE-1; j>=0; j--){
		temps << "0 ";
	}
	for (int j=0; j<MESH_SIZE; j++){
  	temps << "0 ";
	}
	temps <<  endl;
	for (int i=0; i<MESH_SIZE2; i++){
		if(i%2 == 0){
			for (int j=MESH_SIZE-2; j>=0; j--){
        temps << v[i*MESH_SIZE + j] << " ";
			}
			temps << "0 ";
			for (int j=0; j<MESH_SIZE; j++){
        temps << v[i*MESH_SIZE + j] << " ";
			}
			temps <<  endl;
		}else{
			for (int j=MESH_SIZE-1; j>=0; j--){
        temps << v[i*MESH_SIZE + j] << " ";
			}
			for (int j=0; j<MESH_SIZE; j++){
        temps << v[i*MESH_SIZE + j] << " ";
			}
			temps <<  endl;
		}	
	}
	temps << "];" << endl;
	/*temps << "v = [";
				for(int i = 0; i < MESH_SIZE*MESH_SIZE2; ++i) {
						if(i%MESH_SIZE != 0 ){temps << ", ";}
            temps <<  v[i];
						if(i%MESH_SIZE == MESH_SIZE -1){temps << endl;}
        }
        temps << "];" << endl;*/
	temps << "figure;" << endl;
	temps << "plotVT(v);" << endl;

	temps.close();
}




