/*
	V-Scheme and FMG MultiGrid Method
		
	Assignment 2: ME670 - Advanced CFD
	Nirmal S. [234103107]

*/

#include <stdio.h>
#include <math.h>
#include <string.h>

void read_inputs(char file_name[50], int *N_grid, char *method, int *nu1, int *nu2, double *k, double *sigma,
				double *epsilon, int *max_iter, double *J_weight, char *MG_method, int *nu0, int *FMG_level);
int index_level(int Levels, int lvl);
int	num_grid(int Levels, int lvl);
void relaxation_methods(int M, int lvl, int Levels, double v[M], double f[M], double sigma, int nu, char method, double J_weight, double h);
void restriction_methods(int M, int lvl, int Levels, double f[M]);
void prolongation_methods(int M, int lvl, int Levels, double v[M], double v_temp[M]);
void post_processing_methods(int M, int lvl, int Levels, double h, double sigma, double f[M], double v[M], double u[M], double *error, double *res);
void v_cycle(int M, int Levels, int lvl, double v[M], double f[M], double residual[M], double u[M],
		    double C, double k, double sigma, int nu1, int nu2, double h, char method, double J_weight);
				
int main(){
	
	char ip_file_name[50] = "input.txt";	// File name of the txt for inputs
	char op_file_name[50];
	
	// *** Variables that are read from the input documents *** 
	int N_grid, nu1, nu2, nu0, max_iter, FMG_level;
	double J_weight, k, sigma, epsilon;
	char method, MG_method;
	
	// *** Read the input file and extract the details of problem *** 
	read_inputs(ip_file_name, &N_grid, &method, &nu1, &nu2, &k, &sigma,
				&epsilon, &max_iter, &J_weight, &MG_method, &nu0, &FMG_level);
	
	
	// *** Defining parameters and declaring f and v variables ***
	int Levels = log(N_grid)/log(2);
	int M = index_level(Levels, 0); // Array Total Length
	double f[M], v[M], residual[M], u[M];
	for(int i=0; i<M; i++){
		f[i] = 0;
		v[i] = 0;
		residual[i] = 0;
		u[i] = 0;
	}
	double h = 1/(float)N_grid;
	double C = pow(M_PI,2)*pow(k,2) + sigma;	
	double error, res;
	
	char text_method[20];
	if (method == 'G') strcpy(text_method, "Gauss-Seidel");
	else if (method == 'J' && J_weight == 1) strcpy(text_method, "Jacobi");
	else if (method == 'J' && J_weight != 1) strcpy(text_method, "Weighted-Jacobi");
	
	
	if(MG_method=='V'){
		
		sprintf(op_file_name, "Output_k-%.0f %s .txt", k, text_method);
		
		FILE *fp = fopen(op_file_name, "w");
	
		if(fp==NULL){
			printf("Could not find %s\n", op_file_name);
		}
		
		printf("\n\n\t\t\t ######  Result V-Cycle #####\n\n");
		printf("\t Grid_Size: %d\t\t k: %.0f\t\t sigma: %.0f\n", N_grid, k, sigma);
		printf("\t nu1: %d\t nu2: %d\t Relaxation_Method: %s\t", nu1, nu2, text_method);
		if(method == 'J' && J_weight != 1) printf("Weight: %.4f\n", J_weight);
		else printf("\n");
		
		fprintf(fp, "\n\n\t\t\t ######  Result V-Cycle #####\n\n");
		fprintf(fp, "\t Grid_Size: %d\t\t k: %.0f\t\t sigma: %.0f\n", N_grid, k, sigma);
		fprintf(fp, "\t nu1: %d\t nu2: %d\t Relaxation_Method: %s\t", nu1, nu2, text_method);
		if(method == 'J' && J_weight != 1) fprintf(fp, "Weight: %.4f\n", J_weight);
		else fprintf(fp, "\n");
		
		int lvl, iter=0;
		
		do{
			lvl = 1;
			res = 0;
			error = 0;
			
			iter += 1;
			
			v_cycle(M, Levels, lvl, v, f, residual, u,
					C, k, sigma, nu1, nu2, h, method, J_weight);
			
			post_processing_methods(M, lvl, Levels, lvl*h, sigma, f, v, u, &error, &res);
			
			printf("\n\t Iteration:%d \t\tError: %.4f  \t\tResidual: %.8f", iter, error, res);
			fprintf(fp, "\n\t Iteration:%d \t\tError: %.4f  \t\tResidual: %.8f", iter, error, res);
			
			
			/*
			printf("\n\n*** u *** \n");
			for(int i=0; i<M; i++) {
				if(i==3 || i==8 || i==17|| i==34|| i==67|| i==132|| i==261|| i==518|| i==1031|| i==2056) printf("\n");
				printf("\t%.4f", u[i]);
			}
					
			printf("\n\n*** v *** \n");
			for(int i=17; i<M; i++) {
				if(i==3 || i==8 || i==17|| i==34|| i==67|| i==132|| i==261|| i==518|| i==1031|| i==2056) printf("\n");
				printf("\t%.4f", v[i]);
			}
			*/
			
		}while(res>epsilon && iter<max_iter);
		
		printf("\n");
		fclose(fp);
		
	}
	
	if(MG_method=='F'){
		
		sprintf(op_file_name, "Output_FMG Level-%d.txt", FMG_level);
		FILE *fp = fopen(op_file_name, "w");
	
		if(fp==NULL){
			printf("Could not find %s\n", op_file_name);
		}
		
		printf("\n\n\t\t ######  Result Full-Multigrid #####\n\n");
		printf("\t Grid_Size: %d\t\t k: %.0f\t\t sigma: %.0f\n", N_grid, k, sigma);
		printf("\t nu1: %d\t nu2: %d\t Relaxation_Method: %s\t", nu1, nu2, text_method);
		if(method == 'J' && J_weight != 1) printf("Weight: %.4f\n", J_weight);
		else printf("\n");
		printf("\t nu0: %d\t\t FMG_level: %d\n", nu0, FMG_level);
		
		fprintf(fp, "\n\n\t\t ######  Result Full-Multigrid #####\n\n");
		fprintf(fp, "\t Grid_Size: %d\t\t k: %.0f\t\t sigma: %.0f\n", N_grid, k, sigma);
		fprintf(fp, "\t nu1: %d\t nu2: %d\t Relaxation_Method: %s\t", nu1, nu2, text_method);
		if(method == 'J' && J_weight != 1) fprintf(fp, "Weight: %.4f\n", J_weight);
		else fprintf(fp, "\n");
		fprintf(fp, "\t nu0: %d\t\t FMG_level: %d\n", nu0, FMG_level);
		
		
		// ** Relaxing at the level we want to start Full MultiGrid
		double v_temp[M];	// temp variable used in V-cycle while going up the levels
		for(int i=0; i<M; i++) v_temp[i] = 0;
		int h1 = h*pow(2, FMG_level-1);
		relaxation_methods(M, FMG_level, Levels, v, f, sigma, nu1, method, J_weight, h1);
		
		for(int FMG_lvl_i=FMG_level; FMG_lvl_i>1; FMG_lvl_i--){
			
			// *** Prolongate the 'v' values to a finer grid
			prolongation_methods(M, FMG_lvl_i, Levels, v, v_temp);
			
			for(int i=0; i<nu0; i++){
				
				// *** Apply V-Cycle from the finer Grid 'nu0' times
				v_cycle(M, FMG_lvl_i, 1, v, f, residual, u,
						C, k, sigma, nu1, nu2, h, method, J_weight);
			}
		}
		
			
		post_processing_methods(M, 1, Levels, h, sigma, f, v, u, &error, &res);
		printf("\n\t After 1 Full-MultiGrid Cycle:");
		printf("\n\t\t Error: %.4f  \t\tResidual: %.8f", error, res);
		fprintf(fp, "\n\t After 1 Full-MultiGrid Cycle:");
		fprintf(fp, "\n\t\t Error: %.4f  \t\tResidual: %.8f", error, res);
			
			
		// Further Performing V Cycles until residual is below epsilon
		int lvl, iter=0;
		do{
			lvl = 1;
			res = 0;
			error = 0;
			
			iter += 1;
			
			v_cycle(M, Levels, lvl, v, f, residual, u,
					C, k, sigma, nu1, nu2, h, method, J_weight);
			
			post_processing_methods(M, lvl, Levels, lvl*h, sigma, f, v, u, &error, &res);
			
		}while(res>epsilon && iter<max_iter);
		
		printf("\n\n\t Number of V-Cycle Iteration Needed after FMG: %d\n\n", iter);
		fprintf(fp, "\n\n\t Number of V-Cycle Iteration Needed after FMG: %d", iter);
		
			
		printf("\n");
		fclose(fp);
			
	}
	
		
	return 0;
}


int index_level(int Levels, int lvl){return pow(2, Levels-lvl+1) - 2 + (Levels-lvl);}

int	num_grid(int Levels, int lvl){return pow(2,Levels-lvl+1) + 1;}

void read_inputs(char file_name[50], int *N_grid, char *method, int *nu1, int *nu2, double *k, double *sigma,
				double *epsilon, int *max_iter, double *J_weight, char *MG_method, int *nu0, int *FMG_level){
	
	FILE *fp = fopen(file_name, "r");
	
	if(fp==NULL){
		printf("Could not find %s\n", file_name);
		return;
	}
	
	char temp[10];
	
	fscanf(fp, "\n\t##### Inputs to Multigrid Problem #####\n\n\n");
	
	fscanf(fp, "\t*** Grid Inputs ***\n");
	fscanf(fp, "\t\tGrid Size(n): %d\n\n", N_grid);
	
	fscanf(fp, "\t*** Relaxation Method ***\n");
	fscanf(fp, "\t\tMethod(J-Weighted_Jacobi, G-Gauss_Seidel): %c\n", method);
	fscanf(fp, "\t\tJacobi_Weight(w): %s\n\n", &temp);
	
	fscanf(fp, "\t*** V-Cycle Inputs ***\n");
	fscanf(fp, "\t\tnu1: %d\n", nu1);
	fscanf(fp, "\t\tnu2: %d\n\n", nu2);

	fscanf(fp, "\t*** Initial Value Constants ***\n");
	fscanf(fp, "\t\tk: %lf\n", k);
	fscanf(fp, "\t\tsigma: %lf\n\n", sigma);

	fscanf(fp, "\t*** Termination Condition ***\n");
	fscanf(fp, "\t\tEpsilon[1x10^value]: %lf\n", epsilon);
	fscanf(fp, "\t\tMax Iteration: %d\n\n", max_iter);
	
	fscanf(fp, "\t*** Multigrid Method ***\n");
	fscanf(fp, "\t\tMethod(V-Cycle(V), Full-Multigrid(F)): %c\n", MG_method);
	fscanf(fp, "\t\tnu0: %d", nu0);
	fscanf(fp, "\t\tFMG_level: %d", FMG_level);
	
	char *slashPos = strchr(temp, '/');
	if(slashPos != NULL){
		*slashPos = '\0';
		double numerator = atof(temp);
        double denominator = atof(slashPos + 1);
		*J_weight = (double)numerator/(double)denominator;
	}
	else *J_weight = atof(temp);
	
	*epsilon = pow(10,*epsilon);
	
	
	fclose(fp);
	return;
}

void relaxation_methods(int M, int lvl, int Levels, double v[M], double f[M], double sigma, int nu, char method, double J_weight, double h){
	
	int index_lvl = index_level(Levels, lvl);
	int m = index_lvl + num_grid(Levels, lvl);
	double h2 = pow(h,2);
	
	if(method == 'J') {
		int n = pow(2,Levels-lvl+1) + 1;
		double v_new[n];
		for(int i=0; i<n; i++) v_new[i] = v[index_lvl+i];
	
		for(int iter=0; iter<nu; iter++){
			for(int i=1; i<n-1; i++){
				v_new[i] = (1-J_weight)*v_new[i] + J_weight*(h2*f[index_lvl+i] + v[index_lvl+i-1] + v[index_lvl+i+1])/(2+sigma*h2);			
			}
			for(int i=0; i<n; i++)  v[index_lvl+i] = v_new[i];
		}
	}
	
	if(method == 'G') {
		for(int iter=0; iter<nu; iter++){
			for(int i=index_lvl+1; i<m-1; i++){
				v[i] = (h2*f[i] + v[i-1] + v[i+1])/(2+sigma*h2);
			}
		}
	}
	
	return;
}

void restriction_methods(int M, int lvl, int Levels, double f[M]){
	
	int index_lvl = index_level(Levels, lvl);
	int index_lvl_nxt = index_level(Levels, lvl+1);
	int n = num_grid(Levels, lvl);
	
	for(int i=1; i<n/2; i++){
		f[index_lvl_nxt+i] = 0.25*(f[index_lvl+2*i-1] + 2*f[index_lvl+2*i] + f[index_lvl+2*i+1]);
	}
	return;
}

void prolongation_methods(int M, int lvl, int Levels, double v[M], double v_temp[M]){
	
	int index_lvl = index_level(Levels, lvl);
	int index_lvl_prev = index_level(Levels, lvl-1);
	int n = num_grid(Levels, lvl-1);
		
	for(int i=0; i<n/2; i++){
		v_temp[index_lvl_prev+2*i] = v[index_lvl+i];
		v_temp[index_lvl_prev+2*i+1] = 0.5*(v[index_lvl+i]+v[index_lvl+i+1]);
	}	
	
	
	return;
}

void post_processing_methods(int M, int lvl, int Levels, double h, double sigma, double f[M], double v[M], double u[M], double *error, double *res){
	
	int index_lvl = index_level(Levels, lvl);
	int n = num_grid(Levels, lvl);
	double h2 = pow(h,2);
	
	for(int i=index_lvl+1; i<index_lvl+n-1; i++){
		*res += pow(f[i] - 1/h2*(- v[i-1] - v[i+1] + (2+sigma*h2)*v[i]), 2);
		*error += pow(u[i] - v[i], 2);
	}
	*res = sqrt(*res);
	*error = sqrt(*error);
	
	return;
}

void v_cycle(int M, int Levels, int lvl, double v[M], double f[M], double residual[M], double u[M],
		    double C, double k, double sigma, int nu1, int nu2, double h, char method, double J_weight){
	
	double v_temp[M];	// temp variable used in V-cycle while going up the levels
	for(int i=0; i<M; i++) v_temp[i] = 0;
	for(int i=0; i<index_level(Levels, 1); i++) v[i] = 0;
	
	int index_lvl = index_level(Levels, lvl);
	int n = num_grid(Levels, lvl);
	double h1, h2;
	
	// *** Initialising 'f' (RHS) array and exact solution 'u' ***
	for(int i=index_lvl; i<index_lvl+n; i++) f[i] = C*sin((float)k*M_PI*(i-index_lvl)*h);
	for(int i=index_lvl; i<index_lvl+n; i++) u[i] = C/(pow(M_PI,2)*pow(k,2) + sigma)*sin((float)k*M_PI*(i-index_lvl)*h);
	
	
	// *** Moving towards coarser grids ***
	for(int lvl_i=lvl; lvl_i<Levels; lvl_i++){
		
		h1 = h*pow(2,lvl_i-1);
		h2 = pow(h1,2);
		
		// Applying relaxation method to next level
		relaxation_methods(M, lvl_i, Levels, v, f, sigma, nu1, method, J_weight, h1);

		// Calculating the residual
		index_lvl = index_level(Levels, lvl_i);
		n = num_grid(Levels, lvl_i);
		for(int i=index_lvl+1; i<index_lvl+n-1; i++)	residual[i] = f[i] - 1/h2*(- v[i-1] - v[i+1] + (2+sigma*h2)*v[i]);
		
		// Taking residual to next level
		restriction_methods(M, lvl_i, Levels, residual);
		index_lvl = index_level(Levels, lvl_i+1);
		n = num_grid(Levels, lvl_i+1);
		
		// Storing residual back into 'f' at the next level
		for(int i=index_lvl+1; i<index_lvl+n-1; i++)	f[i] = residual[i];
			
	}
	
	
	// *** Solve the equation at the coarsest grid ***
	h1 = h*pow(2,Levels-1);
	relaxation_methods(M, Levels, Levels, v, f, sigma, nu1, method, J_weight, h1);
	
	
	// *** Moving towards finer grids ***
	for(int lvl_i=Levels; lvl_i>lvl; lvl_i--){
			
		// Interpolate to the previous level
		prolongation_methods(M, lvl_i, Levels, v, v_temp);	
		
		// Recalculate the 'v' value with the interpolated value
		index_lvl = index_level(Levels, lvl_i-1);
		n = num_grid(Levels, lvl_i-1);
		for(int i=index_lvl+1; i<index_lvl+n-1; i++)	v[i] += v_temp[i];
		
		// Relaxing the equation at the previous level
		h1 = h*pow(2,lvl_i-2);
		relaxation_methods(M, lvl_i-1, Levels, v, f, sigma, nu2, method, J_weight, h1);
		
	}

	return;
}

















