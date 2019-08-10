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
#define PLATE_WIDTH      			0.005
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

void initV(vector<double> &V);
void fillbounds(vector<double>& lb,vector<double>& hb);
void computeK(const double* V,double K[MESH_SIZE2][MESH_SIZE]);
double Vcon(const double* v);
double V0(); 
double J(VectorXd T);
SparseMatrix<double, ColMajor> derMrespVtimesT(const double *v, VectorXd T);
void generateSystem(double K[MESH_SIZE2][MESH_SIZE], SparseMatrix<double, ColMajor> &M,VectorXd &bb);
void dump2matlab(string filename, VectorXd sol, const double* v);
double mySumconstraint(unsigned n, const double *v, double *grad, void *data);
double myfunc(unsigned n, const double *v, double *grad, void *my_func_data);

int main()
{	
	nlopt::opt opt(nlopt::LD_MMA, MESH_SIZE2*MESH_SIZE); //Create the optimizer object with all the Vs as decision variables
	vector<double> lb(MESH_SIZE2*MESH_SIZE, 0); //Initializing the lower bound vector
	vector<double> ub(MESH_SIZE2*MESH_SIZE, 1); //Initializing the upper bound vector
	fillbounds(lb,ub);//Setting the lower and upper bounds for all the points

	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	opt.set_min_objective(myfunc, NULL);//Setting the objective function + gradient
	opt.add_inequality_constraint(mySumconstraint, NULL, 1e-8);//Setting the inequality constraint + gradient
	opt.set_xtol_rel(1e-4); //Stop criterium
	opt.set_ftol_rel(1e-4);
	std::vector<double> v(MESH_SIZE2*MESH_SIZE);//Initial value for the result vector.
	initV(v);
#ifndef NDEBUG
cout << "Initial v = [" << v[0];
for(int i = 1; i < MESH_SIZE*MESH_SIZE2; ++i) {
	if(i%MESH_SIZE != 0 ){cout << ", ";}
	cout <<  v[i];
	if(i%MESH_SIZE == MESH_SIZE -1){cout << endl;}
}
cout << ']' << endl;
#endif

	double minf;
	try{
		auto t_start = std::chrono::high_resolution_clock::now();
  	  	nlopt::result result = opt.optimize(v, minf);
		auto t_end = std::chrono::high_resolution_clock::now();
		double	elapsed_time=std::chrono::duration<double>(t_end-t_start).count(); 
		//Dumping the results
		std::cout << "found minimum = "<< std::setprecision(10) << minf << std::endl;
		std::cout << "Time taken: " << elapsed_time     << " seconds" << std::endl;
			
		double vv[MESH_SIZE2*MESH_SIZE] = {0};			
		for(int i = 0; i < MESH_SIZE*MESH_SIZE2; ++i) {
			if(v[i] != 0){vv[i] = v[i];}
		}

		double KMesh [MESH_SIZE2][MESH_SIZE] = {0};
		computeK(vv,KMesh);
		//Setting up the system of equations
		SparseMatrix<double, ColMajor> M(TEMP_SIZE,TEMP_SIZE);
		VectorXd bb(TEMP_SIZE), Temp(TEMP_SIZE);
		generateSystem(KMesh,M,bb);
#ifndef NDEBUG
	cout << endl << "Final A matrix is below " << endl;
	cout << MatrixXd(M) << endl;
#endif
		//Building and calling the solver		
		SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver;
		solver.analyzePattern(M);
		solver.compute(M);
		Temp = solver.solve(bb);
		string s = "temps";
		s = s + to_string(MESH_SIZE)  + ".m";

		dump2matlab(s,Temp,vv);
			
	}
	catch(std::exception &e) {
	    std::cout << "nlopt failed: " << e.what() << std::endl;
	}


}


//Sets the V matrix as 1 for the points in the heat sink and as the highest possible constant value for the rest of the pointa
void initV(vector<double>& V){
	double fill = V0();
	int i = 0;
		while (i<2*(metal-1)){
			for (int j=0; j<MESH_SIZE-1; j++){
		      V[i*MESH_SIZE+j] = fill;
			}
			if(i%2 != 0){
				V[i*MESH_SIZE+MESH_SIZE-1] = 1;
			}
			i++;
		}
		while (i<MESH_SIZE2){
			for (int j=0; j<MESH_SIZE-1; j++){
		      V[i*MESH_SIZE+j] = fill;
			}
			if(i%2 != 0){
				V[i*MESH_SIZE+MESH_SIZE-1] = fill;
			}
			i++;
		}
}

//Sets the bounds for the V values
//TODO: Tidy this by making use of the default values.

void fillbounds(vector<double>& lb,vector<double>& hb){
	int i = 0;
		while (i<2*(metal-1)){
			for (int j=0; j<MESH_SIZE-1; j++){
		      lb.at(i*MESH_SIZE+j) = 0;
					hb.at(i*MESH_SIZE+j) = 1;
			
			}
			if(i%2 != 0){
				lb[i*MESH_SIZE+MESH_SIZE-1] = 1;
				hb[i*MESH_SIZE+MESH_SIZE-1] = 1;
			}else{
				lb[i*MESH_SIZE+MESH_SIZE-1] = 0;
				hb[i*MESH_SIZE+MESH_SIZE-1] = 0;
			}
			i++;
		}
		while (i<MESH_SIZE2){
			for (int j=0; j<MESH_SIZE-1; j++){
		      lb[i*MESH_SIZE+j] = 0;
					hb[i*MESH_SIZE+j] = 1;
			}
			if(i%2 != 0){
				lb[i*MESH_SIZE+MESH_SIZE-1] = 0;
				hb[i*MESH_SIZE+MESH_SIZE-1] = 1;
			}else{
				lb[i*MESH_SIZE+MESH_SIZE-1] = 0;
				hb[i*MESH_SIZE+MESH_SIZE-1] = 0;
			}
			i++;
		}
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

//computes the overall concentration of metal in the plate
double Vcon(const double* V){
	double con = 0;
  double tot = (MESH_SIZE2-1)*(MESH_SIZE-1) + 0.5*(MESH_SIZE2+MESH_SIZE);
	for (int j=0; j<MESH_SIZE-1; j++){
        	con += V[j];	
		con += V[(MESH_SIZE2-1)*MESH_SIZE+j]/2;
	}
	for (int i=2; i<MESH_SIZE2-2; i=i+2){
		for (int j=0; j<MESH_SIZE-1; j++){
			con += V[i*MESH_SIZE+j];
		}
	}
	for (int i=1; i<MESH_SIZE2-1; i=i+2){
		con += (V[i*MESH_SIZE]);
		con += (V[i*MESH_SIZE+MESH_SIZE-1])/2;
		for (int j=1; j<MESH_SIZE-1; j++){
			con += V[i*MESH_SIZE+j];
		}
	}

cout << "V concentration = " << con/tot << endl;
  return con/tot;
}

//Computes the maximum average concentration of metal for the points outside the heatsink.
double V0(){
	double tot= (MESH_SIZE2-1)*(MESH_SIZE-1) + 0.5*(MESH_SIZE2+MESH_SIZE);
  double available = maxC*tot;
	available = available-((metal-1)/2.0);
	return available/(tot-((metal-1)/2));
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



//Function for the sumation constraint on the Vs
double mySumconstraint(unsigned n, const double *v, double *grad, void *data)
{
    if (grad) {/*Code the gradient for the restrictions*/
				int i = 0;
			for (int j=0; j<MESH_SIZE-1; j++){
					grad[i*MESH_SIZE+j] = 1;
				}
				grad[i*MESH_SIZE+MESH_SIZE-1] = 0;
				i++;
			while (i<2*(metal-1)){
				for (int j=1; j<MESH_SIZE-1; j++){
					grad[i*MESH_SIZE+j] = 1;
				}
				if(i%2 != 0){
					grad[i*MESH_SIZE] = 0.5;
					grad[i*MESH_SIZE+MESH_SIZE-1] = 0;
				}else{
					grad[i*MESH_SIZE] = 1;
					grad[i*MESH_SIZE+MESH_SIZE-1] = 0;
				}
				i++;
			}	
			while (i<MESH_SIZE2-1){
				
				for (int j=1; j<MESH_SIZE-1; j++){
		      grad[i*MESH_SIZE+j] = 1;
				}
				if(i%2 != 0){
					grad[i*MESH_SIZE] = 1;
					grad[i*MESH_SIZE+MESH_SIZE-1] = 0.5;
				}else{
					grad[i*MESH_SIZE] = 1;
					grad[i*MESH_SIZE+MESH_SIZE-1] = 0;
				}
				i++;
			}
			for (int j=0; j<MESH_SIZE-1; j++){
					grad[i*MESH_SIZE+j] = 0.5;
				}
				grad[i*MESH_SIZE+MESH_SIZE-1] = 0;
    }
		/*Restrictions Constraint(v) <= 0*/
    return Vcon(v) - 0.4;
 } 


//Computes the derivative of the PDE system matrix with respect to v times T

//TODO: check if the use of MESH_WIDTH here is correct


SparseMatrix<double, ColMajor> derMrespVtimesT(const double *v, VectorXd T){
//For each partial derivative of the system matrix
//we find the relevant T indexes,
//the relevant coefficients and we multiply them
		SparseMatrix<double, ColMajor> derv(TEMP_SIZE,MESH_SIZE2*MESH_SIZE);
		int T1,T2;
    //double k = power*(Km-Kp)/(MESH_WIDTH*MESH_WIDTH); actually correct, but we simplify

		double k = power*(Km-Kp);
		int vi = 0;
		while (vi < MESH_SIZE-2){ //Vs in the top row that dont touch the heat sink
				T1= vi;
				T2= T1 + 1;
				derv.insert(T1,vi) =  k*pow(v[vi],power-1)*(T(T2)-T(T1));
				derv.insert(T2,vi) =  k*pow(v[vi],power-1)*(T(T1)-T(T2));
				vi++;
		}
		//LastV in the top row, touches the heatsink 
		T1= vi;
		T2= T1 + 1;
		derv.insert(T1,vi) = k*pow(v[vi],power-1)*(T(T2)-T(T1));
		vi++;
		//Now we skip the fictional point
		vi++;
		while ((vi/MESH_SIZE) < 2*metal){//While there is a heat sink at the end of the row
			if ((vi/MESH_SIZE)%2 == 0){//horizontal Vs
				while (vi%MESH_SIZE < MESH_SIZE-2){ //Vs not touching the heatsink
					T1= (((vi/MESH_SIZE)/2)*MESH_SIZE) + (vi%MESH_SIZE);
					T2= T1 + 1;
					derv.insert(T1,vi) =  k*pow(v[vi],power-1)*(T(T2)-T(T1));
					derv.insert(T2,vi) =  k*pow(v[vi],power-1)*(T(T1)-T(T2));
					vi++;
				}
				//Second to last V in the row, touches the heatsink 
				T1= (((vi/MESH_SIZE)/2)*MESH_SIZE) + (vi%MESH_SIZE);
				T2= T1 + 1;
				derv.insert(T1,vi) =  k*pow(v[vi],power-1)*(T(T2)-T(T1));
				vi++;
				//Last V of the row does not actually exist
				vi++;
			}else{//Vertical Vs
				while (vi%MESH_SIZE < MESH_SIZE-1){ //Vs not touching the heatsink
					T1= (((vi/MESH_SIZE)/2)*MESH_SIZE) + (vi%MESH_SIZE);
					T2= T1 + MESH_SIZE;
					derv.insert(T1,vi) =  k*pow(v[vi],power-1)*(T(T2)-T(T1));
					derv.insert(T2,vi) =  k*pow(v[vi],power-1)*(T(T1)-T(T2)); 
					vi++;
				}
				//Last V of the row is between heat sinks, so it does not affect anything
				//Unless it is at the end of the heat sink
				if((vi/MESH_SIZE)==(2*metal)-1){
					T1= (((vi/MESH_SIZE)/2)*MESH_SIZE) + (vi%MESH_SIZE);
					T2= T1 + MESH_SIZE;
					derv.insert(T2,vi) =  (k*pow(v[vi],power-1)*(T(T1)-T(T2)))/2.0;
				}
				vi++;
			}
		}	
		while ((vi/MESH_SIZE)<MESH_SIZE2-1){//For the next rows
			if ((vi/MESH_SIZE)%2 == 0){//horizontal Vs
				while (vi%MESH_SIZE < MESH_SIZE-1){ //All Vs not touching the heatsink
					T1= (((vi/MESH_SIZE)/2)*MESH_SIZE) + (vi%MESH_SIZE);
					T2= T1 + 1;
					derv.insert(T1,vi) =  k*pow(v[vi],power-1)*(T(T2)-T(T1));
					derv.insert(T2,vi) =  k*pow(v[vi],power-1)*(T(T1)-T(T2));
					vi++;
				}
				//Last V of the row does not actually exist
				vi++;
			}else{//Vertical Vs
				while (vi%MESH_SIZE < MESH_SIZE-1){ //All Vs not touching the heatsink
					T1= (((vi/MESH_SIZE)/2)*MESH_SIZE) + (vi%MESH_SIZE);
					T2= T1 + MESH_SIZE;
					derv.insert(T1,vi) =  k*pow(v[vi],power-1)*(T(T2)-T(T1));
					derv.insert(T2,vi) =  k*pow(v[vi],power-1)*(T(T1)-T(T2));
					vi++;
				}
				//Last V actually exists for this
				T1= (((vi/MESH_SIZE)/2)*MESH_SIZE) + (vi%MESH_SIZE);
				T2= T1 + MESH_SIZE;
				derv.insert(T1,vi) =  (k*pow(v[vi],power-1)*(T(T2)-T(T1)))/2.0;
				derv.insert(T2,vi) =  (k*pow(v[vi],power-1)*(T(T1)-T(T2)))/2.0;
				vi++;
			}
		}
		//Last row
		while (vi%MESH_SIZE < MESH_SIZE-1){ //All Vs not touching the heatsink
					T1= (((vi/MESH_SIZE)/2)*MESH_SIZE) + (vi%MESH_SIZE);
					T2= T1 + 1;
					derv.insert(T1,vi) =  (k*pow(v[vi],power-1)*(T(T2)-T(T1)))/2.0;
					derv.insert(T2,vi) =  (k*pow(v[vi],power-1)*(T(T1)-T(T2)))/2.0;
					vi++;
				}
				//Last V of the row does not actually exist
				vi++;
		return derv;
		
}

//Objective function of the optimization
double myfunc(unsigned n, const double *v, double *grad, void *my_func_data)
{
		/*T computation: We compute the temperatures for the V values using the Eigen library sparse LU solver*/
			double KMesh [MESH_SIZE2][MESH_SIZE] = {0};
			computeK(v,KMesh);
			//Setting up the system of equations
			SparseMatrix<double, ColMajor> M(TEMP_SIZE,TEMP_SIZE);
	  		VectorXd bb(TEMP_SIZE), T(TEMP_SIZE), derGt(TEMP_SIZE), Lambda(TEMP_SIZE);
 			generateSystem(KMesh,M,bb);

			SparseMatrix<double, ColMajor> MT(TEMP_SIZE,TEMP_SIZE); //The transposed matrix will be useful for the gradient
			MT = SparseMatrix<double, ColMajor>(M.transpose()); // / (MESH_WIDTH*MESH_WIDTH);

			//Building and calling the solver		
	 		SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver;
			solver.analyzePattern(M);
			solver.compute(M);
			T = solver.solve(bb);
#ifndef NFDDEBUG
cout << "Computing finite differences..." <<endl;
SparseMatrix<double, ColMajor> Mv(TEMP_SIZE,MESH_SIZE2*MESH_SIZE);
double g [MESH_SIZE2*MESH_SIZE] = {0};
for (int index = 0; index < MESH_SIZE2*MESH_SIZE; index++){
			
double KMeshP [MESH_SIZE2][MESH_SIZE] = {0};
double KMeshN [MESH_SIZE2][MESH_SIZE] = {0};
double vp [MESH_SIZE2*MESH_SIZE] = {0};
double vn [MESH_SIZE2*MESH_SIZE] = {0};
double delta = 0.001;
for(int t = 0; t < MESH_SIZE2*MESH_SIZE; t++){vp[t] = v[t];vn[t] = v[t]; }
vp[index] = v[index] + delta;
vn[index] = v[index] - delta;
computeK(vp,KMeshP);
computeK(vn,KMeshN);
SparseMatrix<double, ColMajor> MP(TEMP_SIZE,TEMP_SIZE);
SparseMatrix<double, ColMajor> MN(TEMP_SIZE,TEMP_SIZE);
VectorXd bp(TEMP_SIZE), TP(TEMP_SIZE),bn(TEMP_SIZE), TN(TEMP_SIZE),TT(TEMP_SIZE);
generateSystem(KMeshP,MP,bp);
generateSystem(KMeshN,MN,bn);

SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solverP,solverN;
solverP.analyzePattern(MP);
solverN.analyzePattern(MN);
solverP.compute(MP);
solverN.compute(MN);
TP = solverP.solve(bp);
TN = solverN.solve(bn);
g[index] = ((J(TP)-J(TN))/(2*delta));
M = (MP - MN)/(2*delta);
cout << "FD with respect to v[" << index/MESH_SIZE << "," << index%MESH_SIZE << "] = " << v[index] << endl;
TT = M*T;
for(int t = 0; t < TEMP_SIZE; t++){Mv.insert(t,index) = TT.coeff(t);}
}
cout << "Computed finite diferences." <<endl;
#endif

#ifndef NDEBUG
	cout << endl << "A matrix is below " << endl;
	cout << MatrixXd(M) << endl;
#endif

#ifndef NDEBUG
					cout << "T = [";
					for (int i=0; i<MESH_SIZE; i++){
						for (int j=0; j<MESH_SIZE; j++){
        			cout << std::setprecision(5) << T((i*MESH_SIZE)+j) << " ";
						}
					cout <<  endl;
					}
					cout << "]" << endl;
#endif

		if (grad) {/*derivative of G(v) using the adjoint method*/
			SparseMatrix<double, ColMajor> MvT(TEMP_SIZE,MESH_SIZE2*MESH_SIZE);
			MvT= derMrespVtimesT(v,T);
			//for (int d = 0; d < TEMP_SIZE; d++){derGt(d) = T(d); }
			for (int d = 0; d < TEMP_SIZE; d++){derGt(d) = T(d) - 293.0; }
			//derivative of the objective function with respect to the temperature
			//Computation of lambda
			SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solverL;
			solverL.analyzePattern(MT);
			solverL.compute(MT);
			Lambda = solverL.solve(derGt);
			bb = (-Lambda.transpose())*MvT;
#ifndef NFDDEBUG
			for (int row = 0; row < TEMP_SIZE; row++){
				for(int col = 0; col < MESH_SIZE*MESH_SIZE2; col++){
					if(abs(MvT.coeff(row,col)-Mv.coeff(row,col))>0.001){
						cout << "Disparitiy in (" << row << "," << col << ") =" << abs(MvT.coeff(row,col)-Mv.coeff(row,col)) << endl;
					}
				}
			}
			cout << MvT-Mv << endl;
			cout << endl;
#endif
		for (int j=0; j<MESH_SIZE2*MESH_SIZE; j++){
			grad[j] = bb(j);
    }
	}



#ifndef NDEBUG
cout  << "v= {";
for(int t = 0; t < MESH_SIZE2*MESH_SIZE-1; t++){
	
cout  << v[t] << ", ";
}
cout << v[MESH_SIZE2*MESH_SIZE-1] << "}" << endl;

cout  << "Adjoint Grad(cost)v= {";
for(int t = 0; t < MESH_SIZE2*MESH_SIZE-1; t++){
	
cout  << grad[t] << ", ";
}
cout  << grad[MESH_SIZE2*MESH_SIZE-1] << "}" << endl << endl;
#endif
#ifndef NFDDEBUG
cout  << "FinDiff Grad(cost)v= {";
for(int t = 0; t < MESH_SIZE2*MESH_SIZE-1; t++){
	
cout  << g[t] << ", ";
}
cout  << g[MESH_SIZE2*MESH_SIZE-1] << "}" << endl << endl;

cout  << "Gradient error= {";
for(int t = 0; t < MESH_SIZE2*MESH_SIZE-1; t++){
	
cout  << grad[t]-g[t] << ", ";
}
cout  << grad[MESH_SIZE2*MESH_SIZE-1]-g[MESH_SIZE2*MESH_SIZE-1] << "}" << endl << endl;
#endif
/*Cost function computation*/
cout << "cost = " << J(T) << endl;

    return J(T); 
}

//Computation of the cost function from the temperatures
double J(VectorXd T){
		double J = 0;
		for (int i = 0; i < TEMP_SIZE; i++){
			//J += (T(i))*(T(i));
			J += (T(i)-293)*(T(i)-293);
		}
		return J/2.0;
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
	temps << "set(gca,'visible','off')" << endl;
	temps << "xlim([0 1]);" << endl;
	temps << "ylim([0 1]);" << endl;
	temps << "plotVT(v);" << endl;
	temps << "saveas(gcf,'temps"<< MESH_SIZE << ".png');" << endl;

	temps.close();
}




