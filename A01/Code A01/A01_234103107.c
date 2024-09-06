/*
	2D Lid Driven Cavity using SIMPLE Algorithm
		
	Assignment 1: ME670 - Advanced CFD
	Nirmal S. [234103107]

*/

#include<stdio.h>
#include<math.h>

void read_Inputs(char file_name[50], double *Re, double *x_length, double *x_no_divisions, double *y_length, double *y_no_divisions, 
				double *u_left_value, double *u_right_value, double *u_top_value, double *u_bottom_value,
				double *v_left_value, double *v_right_value, double *v_top_value, double *v_bottom_value,
				double *u_initialise, double *v_initialise, double *p_initialise,
				double *vel_alpha, double *pressure_alpha,
				double *epsilon, double *max_iter, char *method);
void initialise_Values(int m, int n, double u[m+1][n], double u_star[m+1][n], double d_e[m+1][n],
									 double v[m][n+1], double v_star[m][n+1], double d_n[m][n+1],
									 double p[m+1][n+1], double p_star[m+1][n+1],
									 double pc[m+1][n+1], double b[m+1][n+1], double u_initialise, double v_initialise, double p_initialise);
void pc_initialise(int m, int n, double pc[m+1][n+1]);
void print_Array(int m, int n, double Arr[m][n]);
void u_Boundary_Conditions(int m, int n, double u[m+1][n], double u_left_value, double u_right_value, double u_top_value, double u_bottom_value);
void v_Boundary_Conditions(int m, int n, double v[m][n+1], double v_left_value, double v_right_value, double v_top_value, double v_bottom_value);
void p_Boundary_Conditions(int m, int n, double p[m+1][n+1]);
double max(double a, double b, double c);
void solver_SIMPLE_Algorithm(int m, int n, double u[m+1][n], double u_star[m+1][n], double d_e[m+1][n],
										   double v[m][n+1], double v_star[m][n+1], double d_n[m][n+1],
										   double p[m+1][n+1], double p_star[m+1][n+1],
										   double pc[m+1][n+1], double b[m+1][n+1], 
										   double vel_alpha, double pressure_alpha, double epsilon, double max_iter,
										   double dx, double dy, double Re,
										   double u_left_value, double u_right_value, double u_top_value, double u_bottom_value,
										   double v_left_value, double v_right_value, double v_top_value, double v_bottom_value, char method);
void print_Results(double Re, char u_op_file[50], char v_op_file[50], char tecplot_op_file[50],  double dx, double dy,
					int m, int n, double u_final[m][n], double v_final[m][n], double stream_final[m][n], double vorticity_final[m][n]);
void calculate_Collocated_Grid(int m, int n, double u[m+1][n], double v[m][n+1], double p[m+1][n+1], double u_final[m][n], double v_final[m][n], double p_final[m][n]);
void calculate_Stream_Vorticity(int m, int n, double stream_final[m][n], double vorticity_final[m][n], double u_final[m][n], double v_final[m][n], double dx, double dy);

void main(){
	
	char ip_file_name[50] = "input.txt";	// File name of the txt for inputs
	char u_op_file[50], v_op_file[50], tecplot_op_file[50], method;
	
	// *** Variables that are read from the input documents *** 
	double Re, x_length, x_no_divisions, y_length, y_no_divisions, 
			u_left_value, u_right_value, u_top_value, u_bottom_value,
			v_left_value, v_right_value, v_top_value, v_bottom_value,
			u_initialise, v_initialise, p_initialise,
			vel_alpha, pressure_alpha,
			epsilon, max_iter;
	
	
	// *** Read the input file and extract the details of problem *** 
	read_Inputs(ip_file_name, &Re, &x_length, &x_no_divisions, &y_length, &y_no_divisions, 
				&u_left_value, &u_right_value, &u_top_value, &u_bottom_value,
				&v_left_value, &v_right_value, &v_top_value, &v_bottom_value,
				&u_initialise, &v_initialise, &p_initialise,
				&vel_alpha, &pressure_alpha,
				&epsilon, &max_iter, &method);
	epsilon = pow(10,epsilon);
	


	// *** Calculations and 2D array creation *** 
	double dx = x_length/x_no_divisions;	// Division length
	double dy = y_length/y_no_divisions;
	int n = x_no_divisions + 1;	// No. of points
	int m = y_no_divisions + 1;
	
	// Final Collocated Variables
	double u_final[m][n], v_final[m][n], p_final[m][n], stream_final[m][n], vorticity_final[m][n];
	//Staggered Grid
	double u[m+1][n], u_star[m+1][n], d_e[m+1][n],
		   v[m][n+1], v_star[m][n+1], d_n[m][n+1],
		   p[m+1][n+1], p_star[m+1][n+1],
		   pc[m+1][n+1], b[m+1][n+1]; 
	
	
	// *** Initialisation *** 
	initialise_Values(m, n, u, u_star, d_e, v, v_star, d_n, p, p_star,
					  pc, b, u_initialise, v_initialise, p_initialise);

	
	// *** Apply Boundary Conditions *** 
	u_Boundary_Conditions(m, n, u, u_left_value, u_right_value, u_top_value, u_bottom_value);
	u_Boundary_Conditions(m, n, u_star, u_left_value, u_right_value, u_top_value, u_bottom_value);
	v_Boundary_Conditions(m, n, v, v_left_value, v_right_value, v_top_value, v_bottom_value);
	v_Boundary_Conditions(m, n, v_star, v_left_value, v_right_value, v_top_value, v_bottom_value);
	p_Boundary_Conditions(m, n, p);
	
	// *** Solve using SIMPLE Algorithm ***
	solver_SIMPLE_Algorithm(m, n, u, u_star, d_e, v, v_star, d_n, p, p_star, pc, b, vel_alpha, pressure_alpha, epsilon, max_iter, dx, dy, Re,
							u_left_value, u_right_value, u_top_value, u_bottom_value, v_left_value, v_right_value, v_top_value, v_bottom_value, method);

	// *** Reverting back to Collocated Grid
	calculate_Collocated_Grid(m, n, u, v, p, u_final, v_final, p_final);
	
	// *** Calculating Stream-Function and Vorticity ***
	calculate_Stream_Vorticity(m, n, stream_final, vorticity_final, u_final, v_final, dx, dy);
	
	// *** Calculate min stream function point ***
	double min = 0;
	int i_min, j_min;
	for(int i=1; i<m-1; i++){
		for(int j=1; j<n-1; j++){
			if(stream_final[i][j] < min){
				min = stream_final[i][j];
				i_min = i;
				j_min = j;
			}
		}
	}
	
	printf("Stream min Values\n Stream: %.4f,\tx_val:%.3f,\ty_val:%.3f,\tvorticity: %.4f", min, j_min*dx, (m-i_min)*dy, vorticity_final[i_min][j_min]);
	
	// *** Print Results ***
	print_Results(Re, u_op_file, v_op_file, tecplot_op_file, dx, dy,
					m, n, u_final, v_final, stream_final, vorticity_final);
	
	
	return;
}

void read_Inputs(char file_name[50], double *Re, double *x_length, double *x_no_divisions, double *y_length, double *y_no_divisions, 
				double *u_left_value, double *u_right_value, double *u_top_value, double *u_bottom_value,
				double *v_left_value, double *v_right_value, double *v_top_value, double *v_bottom_value,
				double *u_initialise, double *v_initialise, double *p_initialise,
				double *vel_alpha, double *pressure_alpha,
				double *epsilon, double *max_iter, char *method){
	
	FILE *fp = fopen(file_name, "r");
	
	if(fp==NULL){
		printf("Could not find %s\n", file_name);
		return;
	}
	
	fscanf(fp, "\n\t##### Inputs to 2D Lid Driven Cavity #####\n\n");
	
	fscanf(fp, "\t*** Material Constants ***\n");
	fscanf(fp, "\t\tRe: %lf\n\n", Re);
	
	fscanf(fp, "\t*** Grid Inputs ***\n");
	fscanf(fp, "\t\tx-length: %lf\n", x_length);
	fscanf(fp, "\t\tx-no of divisions(no.): %lf\n", x_no_divisions);
	fscanf(fp, "\t\ty-length: %lf\n", y_length);
	fscanf(fp, "\t\ty-no of divisions(no.): %lf\n\n", y_no_divisions);
	
	fscanf(fp, "\t*** Boundary Conditions ***\n");
	fscanf(fp, "\t\tU-Velocity\n");
	fscanf(fp, "\t\t\tLeft Value: %lf\n", u_left_value);
	fscanf(fp, "\t\t\tRight Value: %lf\n", u_right_value);
	fscanf(fp, "\t\t\tTop Value: %lf\n", u_top_value);
	fscanf(fp, "\t\t\tBottom Value: %lf\n\n", u_bottom_value);
	fscanf(fp, "\t\tV-Velocity\n");
	fscanf(fp, "\t\t\tLeft Value: %lf\n", v_left_value);
	fscanf(fp, "\t\t\tRight Value: %lf\n", v_right_value);
	fscanf(fp, "\t\t\tTop Value: %lf\n", v_top_value);
	fscanf(fp, "\t\t\tBottom Value: %lf\n\n", v_bottom_value);

	fscanf(fp, "\t*** Initialisation Values ***\n");
	fscanf(fp, "\t\tU: %lf\n", u_initialise);
	fscanf(fp, "\t\tV: %lf\n", v_initialise);
	fscanf(fp, "\t\tP: %lf\n\n", p_initialise);

	fscanf(fp, "\t*** Under Relaxation Values ***\n");
	fscanf(fp, "\t\tVelocity: %lf\n", vel_alpha);
	fscanf(fp, "\t\tPressure: %lf\n\n", pressure_alpha);

	fscanf(fp, "\t*** Termination Condition ***\n");
	fscanf(fp, "\t\tEpsilon[1x10^value]: %lf\n", epsilon);
	fscanf(fp, "\t\tMax Iteration: %lf\n\n", max_iter);
	
	fscanf(fp, "\t*** Method to solve Ax=B ***\n");
	fscanf(fp, "\t\tJacobi(J)/Gauss Seidel(G): %s\n", method);
	
	fclose(fp);
	return;
}

void initialise_Values(int m, int n, double u[m+1][n], double u_star[m+1][n], double d_e[m+1][n],
					  double v[m][n+1], double v_star[m][n+1], double d_n[m][n+1],
					  double p[m+1][n+1], double p_star[m+1][n+1],
					  double pc[m+1][n+1], double b[m+1][n+1], double u_initialise, double v_initialise, double p_initialise){
	
	for(int i=0; i<m+1; i++){
		for(int j=0; j<n+1; j++){
			
			if(i<m && j<n){
				u[i][j] = u_initialise;
				u_star[i][j] = u_initialise;
				d_e[i][j] = u_initialise;
				
				v[i][j] = v_initialise;
				v_star[i][j] = v_initialise;
				d_n[i][j] = v_initialise;

				p[i][j] = p_initialise;
				p_star[i][j] = p_initialise;
				pc[i][j] = 0.0;
				}
				
			else if(i==m && j<n){
				u[i][j] = u_initialise;
				u_star[i][j] = u_initialise;
				d_e[i][j] = u_initialise;

				p[i][j] = p_initialise;
				p_star[i][j] = p_initialise;
				pc[i][j] = 0.0;
				}
				
			else if(j==n && i<m){				
				v[i][j] = v_initialise;
				v_star[i][j] = v_initialise;
				d_n[i][j] = v_initialise;

				p[i][j] = p_initialise;
				p_star[i][j] = p_initialise;
				pc[i][j] = 0.0;
				}
			
			else{
				p[i][j] = p_initialise;
				p_star[i][j] = p_initialise;
				pc[i][j] = 0.0;
				}
		}
	}
	return;
}

void print_Array(int m, int n, double Arr[m][n]){
	printf("\nThe Array is: (%d, %d)\n", m, n);
	for(int i=0; i<m; i++){
		printf("\t");
		for(int j=0; j<n; j++){
			printf("\t%.3f",Arr[i][j]);
		}
		printf("\n");
	}
	return;
}

void u_Boundary_Conditions(int m, int n, double u[m+1][n], double u_left_value, double u_right_value, double u_top_value, double u_bottom_value){
	for(int i=0; i<m+1; i++){
		u[i][0] = u_left_value;
		u[i][n-1] = u_right_value; 
	}
	for(int j=0; j<n; j++){
		u[0][j] = u_top_value*2 - u[1][j];
		u[m][j] = u_bottom_value*2 - u[m-1][j];
	}
	return;		
}

void v_Boundary_Conditions(int m, int n, double v[m][n+1], double v_left_value, double v_right_value, double v_top_value, double v_bottom_value){
	for(int i=0; i<m; i++){
		v[i][0] = v_left_value*2 - v[i][1];
		v[i][n] = v_right_value*2 - v[i][n-1]; 
	}
	for(int j=0; j<n+1; j++){
		v[0][j] = v_top_value;
		v[m-1][j] = v_bottom_value;
	}
	return;		
}

void p_Boundary_Conditions(int m, int n, double p[m+1][n+1]){
	for(int i=0; i<m+1; i++){
		p[i][0] = p[i][1];
		p[i][n] = p[i][n-1];
	}
	for(int j=0; j<n+1; j++){
		p[0][j] = p[1][j];
		p[m][j] = p[m-1][j];
	}
	return;		
}

double max(double a, double b, double c){
	
	if(a>=b){
		if(a>=c){
			return a;
		}
		else return c;
	}
	else if(b>=c){
		return b;
	}
	else return c;
}

void solver_SIMPLE_Algorithm(int m, int n, double u[m+1][n], double u_star[m+1][n], double d_e[m+1][n],
										   double v[m][n+1], double v_star[m][n+1], double d_n[m][n+1],
										   double p[m+1][n+1], double p_star[m+1][n+1],
										   double pc[m+1][n+1], double b[m+1][n+1],
										   double vel_alpha, double pressure_alpha, double epsilon, double max_iter,
										   double dx, double dy, double Re,
										   double u_left_value, double u_right_value, double u_top_value, double u_bottom_value,
										   double v_left_value, double v_right_value, double v_top_value, double v_bottom_value, char method){
	
	int iter = 0;
	double error_u, error_v;
	
	double Fe, Fw, Fn, Fs, De, Dw, Dn, Ds, aE, aW, aN, aS, aP;
	
	do{
		
		// *** STEP 1 : Solving x and y momentum equation on starred variables ***
		
		// 1.1 X-Momentum Equation Interior:
		for(int i=1; i<m; i++){
			for(int j=1; j<n-1; j++){
				Fe = (u[i][j+1] + u[i][j])/2.0*dy*1;
				Fw = (u[i][j] + u[i][j-1])/2.0*dy*1;
				Fn = (v[i-1][j+1] + v[i-1][j])/2.0*dx*1;
				Fs = (v[i][j+1] + v[i][j])/2.0*dx*1;
				
				De = (1/Re)*(dy*1/dx);
				Dw = (1/Re)*(dy*1/dx);
				Dn = (1/Re)*(dx*1/dy);
				Ds = (1/Re)*(dx*1/dy);
				
				aE = max(-Fe, De - Fe/2.0, 0);
				aW = max( Fw, Dw + Fw/2.0, 0);
				aN = max(-Fn, Dn - Fn/2.0, 0);
				aS = max( Fs, Ds + Fs/2.0, 0);
				
				aP = aE + aW + aN + aS + Fe - Fw + Fn - Fs;
				
				d_e[i][j] = dy*1/aP;
				
				if(method =='G') u_star[i][j] = (aE*u_star[i][j+1] + aW*u_star[i][j-1] + aN*u_star[i-1][j] + aS*u_star[i+1][j])/aP - d_e[i][j]*(p[i][j+1] - p[i][j]);
				if(method =='J') u_star[i][j] = (aE*u[i][j+1] + aW*u[i][j-1] + aN*u[i-1][j] + aS*u[i+1][j])/aP - d_e[i][j]*(p[i][j+1] - p[i][j]);
				
			}
		}

		// 1.2 u-Boundary Conditions
		u_Boundary_Conditions(m, n, u_star, u_left_value, u_right_value, u_top_value, u_bottom_value);
		
		// 1.3 Y-Momentum Equation Interior:
		for(int i=1; i<m-1; i++){
			for(int j=1; j<n; j++){
				Fe = (u[i+1][j] + u[i][j])/2.0*dy*1;
				Fw = (u[i+1][j-1] + u[i][j-1])/2.0*dy*1;
				Fn = (v[i][j] + v[i-1][j])/2.0*dx*1;
				Fs = (v[i+1][j] + v[i][j])/2.0*dx*1;
				
				De = (1/Re)*(dy*1/dx);
				Dw = (1/Re)*(dy*1/dx);
				Dn = (1/Re)*(dx*1/dy);
				Ds = (1/Re)*(dx*1/dy);
				
				aE = max(-Fe, De - Fe/2.0, 0);
				aW = max( Fw, Dw + Fw/2.0, 0);
				aN = max(-Fn, Dn - Fn/2.0, 0);
				aS = max( Fs, Ds + Fs/2.0, 0);
				
				aP = aE + aW + aN + aS + Fe - Fw + Fn - Fs;
				
				d_n[i][j] = dx*1/aP;
				
				if(method =='G') v_star[i][j] = (aE*v_star[i][j+1] + aW*v_star[i][j-1] + aN*v_star[i-1][j] + aS*v_star[i+1][j])/aP - d_n[i][j]*(p[i][j] - p[i+1][j]);
				if(method =='J') v_star[i][j] = (aE*v[i][j+1] + aW*v[i][j-1] + aN*v[i-1][j] + aS*v[i+1][j])/aP - d_n[i][j]*(p[i][j] - p[i+1][j]);
			}
		}


		// 1.4 v-Boundary Conditions
		v_Boundary_Conditions(m, n, v_star, v_left_value, v_right_value, v_top_value, v_bottom_value);

		
		
		// *** STEP 2 : Solving Pressure Correction Equation ***
		
		// 2.1 Initialising pressure correction array to zero
		for(int i=0; i<m+1; i++){
			for(int j=0; j<n+1; j++){
				pc[i][j] = 0.0;
			}
		}

		// 2.2 Pressure Correction Interior
		for(int i=1; i<m; i++){
			for(int j=1; j<n; j++){
				aE = d_e[i][j]*dy*1;
				aW = d_e[i][j-1]*dy*1;
				aN = d_n[i-1][j]*dx*1;
				aS = d_n[i][j]*dx*1;
				
				aP = aE + aW + aN + aS;
				
				b[i][j] =  (u_star[i][j] - u_star[i][j-1])*dy*1 + (v_star[i-1][j] - v_star[i][j])*dx*1;
				
				pc[i][j] = (aE*pc[i][j+1] + aW*pc[i][j-1] + aN*pc[i-1][j] + aS*pc[i+1][j] - b[i][j])/aP;			
			}
		}
		
		
		
		// *** STEP 3 : Correct Pressure and Velocities ***

		error_u = 0.0;
		error_v = 0.0;

		// 3.1 Correcting Pressure Field
		for(int i=1; i<m; i++){
			for(int j=1; j<n; j++){
				p[i][j] = p[i][j] + pressure_alpha*pc[i][j];
			}
		}

		// 3.2 p-Boundary Conditions
		p_Boundary_Conditions(m, n, p);


		// 3.3 Correcting u-velocity
		for(int i=1; i<m; i++){
			for(int j=1; j<n-1; j++){
				error_u += pow(u[i][j] - (u_star[i][j] - vel_alpha*d_e[i][j]*(pc[i][j+1] - pc[i][j])), 2);
				u[i][j] = u_star[i][j] - vel_alpha*d_e[i][j]*(pc[i][j+1] - pc[i][j]);
				u_star[i][j] = u[i][j];
			}
		}

		// 3.4 u-Boundary Conditions
		u_Boundary_Conditions(m, n, u, u_left_value, u_right_value, u_top_value, u_bottom_value);
		u_Boundary_Conditions(m, n, u_star, u_left_value, u_right_value, u_top_value, u_bottom_value);

		// 3.5 Correcting v-velocity
		for(int i=1; i<m-1; i++){
			for(int j=1; j<n; j++){
				error_v += pow(v[i][j] - (v_star[i][j] - vel_alpha*d_n[i][j]*(pc[i][j] - pc[i+1][j])), 2);
				v[i][j] = v_star[i][j] - vel_alpha*d_n[i][j]*(pc[i][j] - pc[i+1][j]);
				v_star[i][j] = v[i][j];
			}
		}

		// 3.6 v-Boundary Conditions
		v_Boundary_Conditions(m, n, v, v_left_value, v_right_value, v_top_value, v_bottom_value);
		v_Boundary_Conditions(m, n, v_star, v_left_value, v_right_value, v_top_value, v_bottom_value);
		
		
		// *** STEP 4 : Print Error ***
		error_u = sqrt(error_u);
		error_v = sqrt(error_v);
		if(iter%100==0) printf("  Iteration %d:\tError_u: %.6lf \tError_v: %.6lf\n", iter, error_u, error_v);
		
		
		// Increase Iteration
		iter += 1;
		
		
		if(iter==max_iter) printf("\n\n *** Did not Converge :( *** \n\n");
		
		if(u[2][2]!=u[2][2]){
			printf("\n\n *** Nan Error Occured ***\n\n");
			break;		
		}
		
		if(error_u<=epsilon && error_v<=epsilon) printf("\n\n *** Converged at Iteration No: %d*** \n\n", iter);
		
	}while(error_u>epsilon && error_v>epsilon && iter<max_iter);
	


	return;
}

void print_Results(double Re, char u_op_file[50], char v_op_file[50], char tecplot_op_file[50], double dx, double dy,
					int m, int n, double u_final[m][n], double v_final[m][n], double stream_final[m][n], double vorticity_final[m][n]){
	
	if(u_final[2][2]!=u_final[2][2]){
		printf("\n\n *** Not Printing the Values ***\n\n");
		return;
	}
	
	sprintf(u_op_file, "U_Re-%d_vertical_mid.csv", (int)Re);
	sprintf(v_op_file, "V_Re-%d_horizontal_mid.csv", (int)Re);
	sprintf(tecplot_op_file, "Tecplot_Re-%d.plt", (int)Re);

	// Writing U-vertical Midsection file
	FILE *fp = fopen(u_op_file, "w");
	if(fp==NULL){
		printf("Could not create %s\n", u_op_file);
		return;
	}
	fprintf(fp, "U\n");
	if(n%2==0) for(int i=0; i<m; i++) fprintf(fp, "%.5f\n", (u_final[i][(int)(n/2.0-1)] + u_final[i][(int)(n/2.0)])/2.0);
	else for(int i=0; i<m; i++) fprintf(fp, "%.5f\n", u_final[i][(int)((n+1)/2.0)]);
	fclose(fp);
	
	// Writing V-vertical Midsection file
	FILE *fp2 = fopen(v_op_file, "w");
	if(fp2==NULL){
		printf("Could not create %s\n", v_op_file);
		return;
	}
	fprintf(fp2, "U\n");
	if(m%2==0) for(int j=0; j<n; j++) fprintf(fp2, "%.5f\n", (v_final[(int)(m/2.0-1)][j] + v_final[(int)(m/2.0)][j])/2.0);
	else for(int j=0; j<n; j++) fprintf(fp2, "%.5f\n", v_final[(int)((m+1)/2.0)][j]);
	fclose(fp2);
	
	// Writing TECPLOT file
	FILE *fp3 = fopen(tecplot_op_file, "w");
	if(fp3==NULL){
		printf("Could not create %s\n", v_op_file);
		return;
	}
	fprintf(fp3, "ZONE I = %d, J = %d\n", m, n);
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			fprintf(fp3, "%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", (double)j*dx, (double)(m-i-1)*dy, u_final[i][j], v_final[i][j], stream_final[i][j], vorticity_final[i][j]);
		}
	}
	fclose(fp3);
	
	return;
	
}

void calculate_Collocated_Grid(int m, int n, double u[m+1][n], double v[m][n+1], double p[m+1][n+1], double u_final[m][n], double v_final[m][n], double p_final[m][n]){
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			u_final[i][j] = (u[i][j] + u[i+1][j])/2.0;
			v_final[i][j] = (v[i][j] + v[i][j+1])/2.0;
			p_final[i][j] = (p[i][j] + p[i+1][j] + p[i][j+1] + p[i+1][j+1])/4.0;
		}
	}
}

void calculate_Stream_Vorticity(int m, int n, double stream_final[m][n], double vorticity_final[m][n], double u_final[m][n], double v_final[m][n], double dx, double dy){
	// Interior Points
	for(int i=1; i<m-1; i++){
		for(int j=1; j<n-1; j++){
			stream_final[i][j] = 0.0;
			vorticity_final[i][j] = (v_final[i][j+1] - v_final[i][j-1])/(2.0*dx) - (u_final[i-1][j] - u_final[i+1][j])/(2.0*dy);
		}
	}
	// Boundary
	for(int i=0; i<m; i++){
		stream_final[i][0] = 0.0;
		stream_final[i][n-1] = 0.0;
		vorticity_final[i][0] = (v_final[i][1] - v_final[i][0])/(dx);
		vorticity_final[i][n-1] = (v_final[i][n-1] - v_final[i][n-2])/(dx);
	}
	for(int j=0; j<n; j++){
		stream_final[0][j] = 0.0;
		stream_final[m-1][j] = 0.0;
		vorticity_final[0][j] = - (u_final[0][j] - u_final[1][j])/(dy);
		vorticity_final[m-1][j] = - (u_final[m-2][j] - u_final[m-1][j])/(dy);
	}
	
	
	double temp;
	double error = 0.0;
	double beta = pow(dx,2)/pow(dy,2);
	do{
	//for(int k=0; k<1000; k++){
		error = 0.0;
		for(int i=1; i<m-1; i++){
			for(int j=1; j<n-1; j++){
				temp = (vorticity_final[i][j]*dx*dx + (stream_final[i][j-1] + stream_final[i][j+1]) + (stream_final[i-1][j] + stream_final[i+1][j])*beta)/(2*(1+beta));
				error += pow(stream_final[i][j] - temp, 2.0);
				stream_final[i][j] = temp;
			}
		}
		error = sqrt(error);
	//}
	} while(error>1e-5);
	
	
}
		
		
		
		
		
		
		
		
		
		
		