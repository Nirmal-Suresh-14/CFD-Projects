/*
	Conjugate Gradient and PreConditioned Conjugate Gradient
		
	Assignment 3: ME670 - Advanced CFD
	Nirmal S. [234103107]

*/

#include <stdio.h>
#include <math.h>
#include <string.h>

void read_inputs(char file_name[50], int *M, int *N, double T_ip[4], int *method, double *epsilon, int *max_iter); 
void initialise_values(int M, int N, int L, double T[L], double An[L], double Aw[L], double Ap[L], double Ae[L], double As[L], double RHS[L], double T_ip[4]);
void op_file_name(int method, char op_iter_residual[50], char op_file_name[50], char op_diagonal[50]);
void calculate_preconditioner(int M, int N, int L, double An[L], double Aw[L], double Ap[L], double Ae[L], double As[L], double Ln[L], 
		double Lnw[L], double Lw[L], double Lp[L], double Up[L], double Ue[L], double Use[L], double Us[L], double method);
void calculate_Minv_r(int M, int N, int L, double r[L], double Mr[L], double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], 
		double Ue[L], double Use[L], double Us[L]);
void initialise_r_d(int M, int N, int L, double r[L], double Mr[L], double d[L], double T[L], double An[L], double Aw[L], double Ap[L], 
		double Ae[L], double As[L], double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], 
		double Ue[L], double Use[L], double Us[L], double RHS[L], int method, double *residual, double *residual_rMr);
void Conjugate_Gradient(int M, int N, int L, double T[L], double An[L], double Aw[L], double Ap[L], double Ae[L], double As[L], 
		double RHS[L], double r[L], double d[L], double Ad[L], double *residual);
void PreConditioned_CG(int M, int N, int L, double T[L], double An[L], double Aw[L], double Ap[L], double Ae[L], double As[L], 
		double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], double Ue[L], double Use[L], double Us[L],
		double RHS[L], double r[L], double Mr[L], double d[L], double Ad[L], double *residual, double *residual_rMr);
double dot_product(int L, double a[L], double b[L]);
void print_iter_res(char op_iter_residual[50], int iter, int max_iter, double iter_res[max_iter][2]);
void print_temp(char op_temp[50], int N, int L, double T[L], double T_ip[4]);
void print_temp_analytic(char op_temp[50], int M, int N);
void print_diag(char op_diagonal[50], int M, int L, double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], double Ue[L], double Use[L], double Us[L]);


int main(){
	
	char ip_file_name[50] = "input.txt";	// File name of the txt for inputs
	char op_iter_residual[50], op_temp[50], op_diagonal[50];
	
	
	// *** Variables that are read from the input documents *** 
	int M, N, method, max_iter;
	double epsilon, T_ip[4];
	
	
	// *** Read the input file and extract the details of problem *** 
	read_inputs(ip_file_name, &M, &N, T_ip, &method, &epsilon, &max_iter);
	//printf("M, N %d, %d, T {%.4f, %.4f, %.4f, %.4f}, Method: %d, Epsilon: %.6f, Max_iter: %d", M, N, T_ip[0], T_ip[1], T_ip[2], T_ip[3], method, epsilon, max_iter);
	
	// *** Creating the required arrays for computation ***
	int L = M*N+1;
	double T[L], An[L], Aw[L], Ap[L], Ae[L], As[L], RHS[L], r[L], Mr[L], d[L], Ad[L];
	double Ln[L], Lnw[L], Lw[L], Lp[L], Up[L], Ue[L], Use[L], Us[L];
	double residual, residual_rMr, iter_res[max_iter][2];
	
	
	// *** Initialising the arrays and output file names ***
	initialise_values(M, N, L, T, An, Aw, Ap, Ae, As, RHS, T_ip);
	op_file_name(method, op_iter_residual, op_temp, op_diagonal);
	
	// *** Calculate the Preconditioner based on the method input ***
	calculate_preconditioner(M, N, L, An, Aw, Ap, Ae, As, Ln, Lnw, Lw, Lp, Up, Ue, Use, Us, method);
	
	// ***	Initialise 'r' and 'd' ***
	initialise_r_d(M, N, L, r, Mr, d, T, An, Aw, Ap, Ae, As, Ln, Lnw, Lw, Lp, Up, Ue, Use, Us, RHS, method, &residual, &residual_rMr);
	
	// *** Running Iterations using the method selected ***
	int iter;
	for(iter=1; iter<max_iter; iter++){
		
		// ***** Applying the Conjugate / PreConditioned CG Method in the iteration
		if(method==1)	Conjugate_Gradient(M, N, L, T, An, Aw, Ap, Ae, As, RHS, r, d, Ad, &residual);
		if(method==2)	PreConditioned_CG(M, N, L, T, An, Aw, Ap, Ae, As, Ln, Lnw, Lw, Lp, Up, Ue, Use, Us, RHS, r, Mr, d, Ad, &residual, &residual_rMr);
		if(method==3||method==4){
			printf("\n\t **** Not coded for ILU and SIP ****\n\n");
			break;
		}
		if(method==5)	break;
		
		if(iter%10==0)	printf("\t\tIter: %d \t Residual: %.6f\n", iter, residual);
		
		iter_res[iter][0] = (float)iter;
		iter_res[iter][1] = residual;
		
		if(residual <= epsilon) break;
	}
	
	printf("\n\t\tFINAL Iter: %d \t Residual: %.6f\n\n", iter, residual);
	
	
	
	// Writing the output data into csv files
	if(method==3){
		print_diag(op_diagonal, M, L, Ln, Lnw, Lw, Lp, Up, Ue, Use, Us);
	}
	if(method!=5){
		print_iter_res(op_iter_residual, iter, max_iter, iter_res);
		print_temp(op_temp, N, L, T, T_ip);
	}
	if(method==5) print_temp_analytic(op_temp, M, N);
	
		
	return 0;
}


void read_inputs(char file_name[50], int *M, int *N, double T_ip[4], int *method, double *epsilon, int *max_iter){
	
	FILE *fp = fopen(file_name, "r");
	
	if(fp==NULL){
		printf("Could not find %s\n", file_name);
		return;
	}
	
	char temp[10];
	
	fscanf(fp, "\n\t##### Inputs to Conjugate Gradient #####\n\n\n");
	
	fscanf(fp, "\t*** Control Volume Inputs ***\n");
	fscanf(fp, "\t\tx-direction(n): %d\n", N);
	fscanf(fp, "\t\ty-direction(m): %d\n\n", M);
	
	fscanf(fp, "\t*** Boundary Conditions ***\n");
	fscanf(fp, "\t\tTemp Top: %lf\n", &T_ip[0]);
	fscanf(fp, "\t\tTemp Bottom: %lf\n", &T_ip[1]);
	fscanf(fp, "\t\tTemp Left: %lf\n", &T_ip[2]);
	fscanf(fp, "\t\tTemp Right: %lf\n\n", &T_ip[3]);
	
	fscanf(fp, "\t*** Solution Method ***\n");
	fscanf(fp, "\t\t1. Conjugate Gradient\n");
	fscanf(fp, "\t\t2. Preconditioned CG, Jacobi\n");
	fscanf(fp, "\t\t3. Preconditioned CG, ILU\n");
	fscanf(fp, "\t\t4. Preconditioned CG, SIP\n");
	fscanf(fp, "\t\t5. Analytical Solution\n");
	fscanf(fp, "\t\tSelection (1, 2, 3, 4, 5): %d\n\n", method);
	
	fscanf(fp, "\t*** Termination Condition ***\n");
	fscanf(fp, "\t\tEpsilon[1x10^value]: %lf\n", epsilon);
	fscanf(fp, "\t\tMax Iteration: %d", max_iter);
	
	*epsilon = pow(10,*epsilon);
	
	
	fclose(fp);
	return;
}

void initialise_values(int M, int N, int L, double T[L], double An[L], double Aw[L], double Ap[L], double Ae[L], double As[L], double RHS[L], double T_ip[4]){
	
	float beta_2 = pow((float)M/N, 2);
	
	for(int i=1; i<=L; i++){
		
		// Default Values
		An[i] = -beta_2;
		Aw[i] = -1.0;
		Ap[i] = 2*(1 + beta_2);
		Ae[i] = -1.0;
		As[i] = -beta_2;
		T[i] = 0.0;
		RHS[i] = 0.0;
		
		// Top Values
		if(i>=1 && i<=N){
			An[i] = 0;
			Ap[i] += beta_2;
			RHS[i] += 2*T_ip[0];
		}
		
		// Bottom Values
		if(i>=(M-1)*N+1 && i<=M*N){
			As[i] = 0;
			Ap[i] += beta_2;
			RHS[i] += 2*T_ip[1];
		}
		
		// Left Values
		if(i%N==1){
			Aw[i] = 0;
			Ap[i] += 1;
			RHS[i] += 2*T_ip[2];
		}
		
		// Right Values
		if(i%N==0){
			Ae[i] = 0;
			Ap[i] += 1;
			RHS[i] += 2*T_ip[3];
		}
	}
	return;
}

void op_file_name(int method, char *op_iter_residual, char *op_temp, char *op_diagonal){
	
	if(method==1){
		strcpy(op_iter_residual, "1.1 Iter_Residual (CG).csv");
		strcpy(op_temp, "1.2 Temp (CG).csv");
	}
		
	if(method==2){
		strcpy(op_iter_residual, "2.1 Iter_Residual (PCG_Jacobi).csv");
		strcpy(op_temp, "2.2 Temp (PCG_Jacobi).csv");
	}
		
	if(method==3){
		strcpy(op_iter_residual, "3.1 Iter_Residual (PCG_ILU).csv");
		strcpy(op_temp, "3.2 Temp (PCG_ILU).csv");
		strcpy(op_diagonal, "3.3 Diagonal Values.csv");
	}
		
	if(method==4){
		strcpy(op_iter_residual, "4.1 Iter_Residual (PCG_SIP).csv");
		strcpy(op_temp, "4.2 Temp (PCG_SIP).csv");
	}
			
	if(method==5){
		strcpy(op_temp, "5. Temp (Analytical).csv");
	}
	
	return;
}

void calculate_preconditioner(int M, int N, int L, double An[L], double Aw[L], double Ap[L], double Ae[L], double As[L], 
		double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], double Ue[L], double Use[L], double Us[L], double method){
	
	for(int i=1; i<L; i++){
		Ln[i] = 0.0;
		Lnw[i] = 0.0;
		Lw[i] = 0.0;
		Lp[i] = 0.0;
		Up[i] = 1.0;
		Ue[i] = 0.0;
		Use[i] = 0.0;
		Us[i] = 0.0;
	}
	
	if(method==2){
		for(int i=1; i<L; i++)		Lp[i] = Ap[i];
	}
	
	if(method==3){
		for(int i=1+N; i<L; i++)	Ln[i] = An[i];
		for(int i=2; i<L; i++)		Lw[i] = Aw[i];
		Lp[1] = Ap[1];
		for(int i=2; i<L; i++){
			if((i-1)<=L-1)	Ue[i-1] = Ae[i-1]/Lp[i-1];
			if((i-1)<=L-N)	Us[i-1] = As[i-1]/Lp[i-1];
			if(i<=N)	Lp[i]	= Ap[i] - Lw[i]*Ue[i-1];
			else		Lp[i]	= Ap[i] - Lw[i]*Ue[i-1] - Ln[i]*Us[i-N];
		}
	}	
		
		
	
	return;
}

void calculate_Minv_r(int M, int N, int L, double r[L], double Mr[L], double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], 
		double Ue[L], double Use[L], double Us[L]){
		      
		      
	// Minv_r = Uinv * Linv * r	(taking Linv * r = y)
	
	double y[L];	// Intermediate temp array to store y
	double sum;
	
	for(int i=1; i<L; i++){		// Fwd Substitution
		sum = Ln[i]*y[i-N] + Lnw[i]*y[i-N+1] + Lw[i]*y[i-1];
		y[i] = r[i]/Lp[i];
	}
	
	for(int i=L-1; i>0; i--){	// Bwd Substitution
		sum = Ue[i]*Mr[i+1] + Use[i]*Mr[i+N-1] + Us[i]*Mr[i+N];
		Mr[i] = r[i] - sum;
	}
	
	return;		      
}

void initialise_r_d(int M, int N, int L, double r[L], double Mr[L], double d[L], double T[L], double An[L], double Aw[L], 
		double Ap[L], double Ae[L], double As[L], double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], 
		double Ue[L], double Use[L], double Us[L], double RHS[L], int method, double *residual, double *residual_rMr){
	
	if(method==1){		// For Conjugate Gradient Method
		for(int i=1; i<L; i++){
			r[i] = RHS[i] - (An[i]*T[i-N] + Aw[i]*T[i-1] + Ap[i]*T[i] + Ae[i]*T[i+1] + As[i]*T[i+N]);
			d[i] = r[i];
		}
		*residual = sqrt(dot_product(L, r, r));
	}
	
	
	
	else{		// For Preconditioned CG Method
		
		
		for(int i=1; i<L; i++)	r[i] = RHS[i] - (An[i]*T[i-N] + Aw[i]*T[i-1] + Ap[i]*T[i] + Ae[i]*T[i+1] + As[i]*T[i+N]);
		calculate_Minv_r(M, N, L, r, Mr, Ln, Lnw, Lw, Lp, Up, Ue, Use, Us);
		
		for(int i=1; i<L; i++)	d[i] = Mr[i];		
		
		*residual = sqrt(dot_product(L, r, r));
		*residual_rMr = sqrt(dot_product(L, Mr, r));
	}
	
	
	
	return;	
}

void Conjugate_Gradient(int M, int N, int L, double T[L], double An[L], double Aw[L], double Ap[L], double Ae[L], double As[L], 
		double RHS[L], double r[L], double d[L], double Ad[L], double *residual){
	
	// Calculating vector A d
	for(int i=1; i<L; i++)	Ad[i] = An[i]*d[i-N] + Aw[i]*d[i-1] + Ap[i]*d[i] + Ae[i]*d[i+1] + As[i]*d[i+N];
	
	// 1. Calculating alpha(i)
	double alpha = pow(*residual,2)/dot_product(L, d, Ad);
	
	// 2. Calculating T(i+1)
	for(int i=1; i<=L; i++)	T[i] += alpha*d[i];
	
	// 3. Calculating r(i+1)
	for(int i=1; i<L; i++)	r[i] -= alpha*Ad[i];
	
	// 4. Calculating Beta(i+1)
	double beta = dot_product(L, r, r)/ pow(*residual,2);
	
	// 5. Calculating d(i+1)
	for(int i=1; i<L; i++)	d[i] = r[i] + beta*d[i];
	
	// 6. Calculating Residual value residual(i+1)
	*residual = sqrt(dot_product(L, r, r));
	
	return;
}

void PreConditioned_CG(int M, int N, int L, double T[L], double An[L], double Aw[L], double Ap[L], double Ae[L], double As[L], 
		double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], double Ue[L], double Use[L], double Us[L],
		double RHS[L], double r[L], double Mr[L], double d[L], double Ad[L], double *residual, double *residual_rMr){
		
	// Calculating vector A d
	for(int i=1; i<L; i++)	Ad[i] = An[i]*d[i-N] + Aw[i]*d[i-1] + Ap[i]*d[i] + Ae[i]*d[i+1] + As[i]*d[i+N];
	
	// 1. Calculating alpha(i)
	double alpha = pow(*residual_rMr,2)/dot_product(L, d, Ad);
		
	// 2. Calculating T(i+1)
	for(int i=1; i<=L; i++)	T[i] += alpha*d[i];
	
	// 3. Calculating r(i+1)
	for(int i=1; i<L; i++)	r[i] -= alpha*Ad[i];
	
	// 4. Calculating Beta(i+1)
	calculate_Minv_r(M, N, L, r, Mr, Ln, Lnw, Lw, Lp, Up, Ue, Use, Us);
	double beta = dot_product(L, Mr, r)/ pow(*residual_rMr,2);
	*residual_rMr = sqrt(dot_product(L, Mr, r));
	
	// 5. Calculating d(i+1)
	for(int i=1; i<L; i++)	d[i] = Mr[i] + beta*d[i];
	
	// 6. Calculating Residual value residual(i+1)
	*residual = sqrt(dot_product(L, r, r));
	
	return;
}

double dot_product(int L, double a[L], double b[L]){
	double value = 0.0;
	for(int i=1; i<L; i++)	value += a[i]*b[i];
	return value;
}

void print_iter_res(char op_iter_residual[50], int iter, int max_iter, double iter_res[max_iter][2]){
	
	FILE *fp = fopen(op_iter_residual, "w");
	
	if(fp==NULL){
		printf("Could not find %s\n", op_iter_residual);
		return;
	}
	
	fprintf(fp, "Iteration, Residual\n");
	for(int i=1; i<=iter; i++)	fprintf(fp, "%d, %.7f\n", (int)iter_res[i][0], iter_res[i][1]);
	
	fclose(fp);
	
	return;
}

void print_temp(char op_temp[50], int N, int L, double T[L], double T_ip[4]){
	
	FILE *fp = fopen(op_temp, "w");
	
	if(fp==NULL){
		printf("Could not find %s\n", op_temp);
		return;
	}

	for(int i=1; i<=N+1; i++)	fprintf(fp, "%.7f, ", T_ip[0]);
	for(int i=1; i<L-N-1; i++){
		if(i%N==1)	fprintf(fp, "\n%.7f, ", T_ip[2]);
		fprintf(fp, "%.7f, ", 0.25*(T[i] + T[i+1] + T[i+N] + T[i+N+1]));
		if(i%N==N-1){
			fprintf(fp, "%.7f", T_ip[3]);
			i += 1;
			}
	}
	fprintf(fp, "\n");
	for(int i=1; i<=N+1; i++)	fprintf(fp, "%.7f, ", T_ip[1]);
	
	fclose(fp);
	
	return;
}

void print_temp_analytic(char op_temp[50], int M, int N){
	
	FILE *fp = fopen(op_temp, "w");
	
	if(fp==NULL){
		printf("Could not find %s\n", op_temp);
		return;
	}
	
	double dx = (double)1/N;
	double dy = (double)1/M;
	
	double temp;
	
	for(int i=1; i<=M+1; i++){
		for(int j=1; j<=N+1; j++){
			temp = 0.0;
			for(int n=1; n<200; n++)	temp += (pow(-1, n+1) + 1)/n*sin(n*M_PI*(j-1)*dx)*sinh(n*M_PI*(1 - (i-1)*dy))/sinh(n*M_PI);
			temp = 2/M_PI*temp;
			fprintf(fp, "%.7f, ", temp);
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	
	return;
}

void print_diag(char op_diagonal[50], int M, int L, double Ln[L], double Lnw[L], double Lw[L], double Lp[L], double Up[L], double Ue[L], double Use[L], double Us[L]){
	
	FILE *fp = fopen(op_diagonal, "w");
	
	if(fp==NULL){
		printf("Could not find %s\n", op_diagonal);
		return;
	}
	
	fprintf(fp, "n, Ln, Lw, Lp, Up, Ue, Us\n");
	
	for(int i=1; i<L; i++){
		fprintf(fp, "%d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", i, Ln[i], Lw[i], Lp[i], Up[i], Ue[i], Us[i]);
	}
	
	fclose(fp);
	
	
	return;
}



