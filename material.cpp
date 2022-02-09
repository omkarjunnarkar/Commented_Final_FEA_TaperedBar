/*
***********************************************************************************************************************************************************

												Finite Element Code for Tapered Bar Problem
															Material Routine

Author: Omkar Junnarkar, Semester-3 MSc. Computational Material Science
Matriculation Nr.: 66157	Email: omkar.junnarkar@student.tu-freiberg.de
IDE : Microsoft Visual Studio 2019

Objective : Handle the Material Point Computations for the Elements in the Finite Element Problem

Libraries Included :

iostream: For Input-Output of the C++ Program
Eigen/Dense: For Using Dense Matrix Interface of Linear Algebra Library Eigen
iomanip: To manipulate the number of decimal places in output
fstream: To create the stream of file and write
element.h: Element Routine header file
material.h: Material Routine header file

*/

#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include"material.h"
#include"element.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*
	To reduce the effort of specifying libraries/class for each command/data type
*/

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*
	Arguments : Stress Value, Yield Stress
	Output : Value of < Stress/Yield - 1 > , where <> are McAuley Brackets
	If in Elastic Zone, returns 0 ;	Else returns computed value
	[ Refer Report/Manual for more Details ]
*/

double Lambda(double sigma, double yield) {
	if (abs(sigma) <= yield) { return 0; }
	else return ((abs(sigma) / yield) - 1);
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*
	tuple : template from 'std' class. Ties and Passes values together ; Used to output multiple values at once

	Arguments :
			eps_p_m : Initial Plastic Strain
			epsilon : Total Strain
			eta : Viscosity
			h : Step Size for Euler Method
			t_total : time step size for force incement
			yield : yield stress
			E: Young's Modulus

	Output :
			Stress
			Plastic Strain
			Tangent Stiffness (Ct)

	Computation Strategy:

	->Divide One time step further into smaller time steps
	->Use initial plastic strain and compute stress, and next plastic strain using that stress
	->Inside Euler implicit, use the computed plastic strain to compute next stress state
	->Using sam estress update, again compute the stress, and update the plastic strain
	->Repeat untill the difference between Previous and Next Plastic Strain is less that 1e-5
	->Update Current stress and Strain
	->This represnts a Radial Return Scheme

	The Tangent Stiffness is derived (for plastic zone) from [ dStress/dStrain ] = E*Yield / (Yield + E*Viscocity*timeStep)
	[ Refer Report/Manual for more Details ]

*/
tuple<double, double, double>get_Stress_Strains(double eps_p_m, double epsilon, double eta, double h, double t_total, double yield, double E) {
	
	int s = (t_total / h);
	double Ct = E;
	double sigma_m,eps_p_next, eps_p_prev, sigma_next;
	
	// 1 step Euler Explicit

	for (int i = 0; i<s; i++) {
		sigma_m = E * (epsilon - eps_p_m);
		eps_p_next = eps_p_m + h * eta * signum(sigma_m) * Lambda(sigma_m, yield);

		// 1000 Steps Euler Implicit

		for (int j = 0; j < 1000; j++) {
			eps_p_prev = eps_p_next;
			sigma_next = E * (epsilon - eps_p_next);
			eps_p_next = eps_p_m + h * eta * signum(sigma_next) * Lambda(sigma_next, yield);

			if (abs(eps_p_next - eps_p_prev) < 1e-5) {
				break;
			}
		}
		
		sigma_m = E * (epsilon - eps_p_next);						
		eps_p_m = eps_p_next;

	}

	if (abs(sigma_m) > yield) {
		Ct = (E * yield) / (yield + E * eta * t_total);		
	}
	else Ct = E;

	return { sigma_m, eps_p_next, Ct };//,sigmas, ep_p_vals
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*
	tuple : template from 'std' class. Ties and Passes values together ; Used to output multiple values at once

	Arguments :
			elements : Number of elements
			gauss_weights : Weight for the Gauss Point
			Lengths : List of Areas
			Areas : List of Lengths
			E_mat : Stiffness(initially=Young's Modulous) of all elements
			u : Global Displacement Matrix
			time_step: Value of single time step
			eta : Viscosity
			eps_pl : Initial Plastic Strain
			yield : yield stress
			E: Young's Modulus



	Output :
			Global Tangent Stiffness (Ct)
			Global Internal Force Matrix (F_int)
			Plastic Strain of all elements (eps_pl)
			Total Strain of all elements (epsilon)
			Stress Matrix of Elements (stress_mat)

	For each element, stresses and strains are computed using get_Stress_Strains().
	F_int_temp and Kt_temp : Temperory matrices for each Element which are computed according to following for 1 Gauss Point
	F_int_element = B_transpose*A*Sigma*J*Weight & Kt_element = B_transpose*B*E*A*J*Weight where J=Length/2, and for 1 Gauss Point, Summation turns out to have weight = 2

	Conversion from Local to Global done using Assignment Matrix : Kt_Global = Summation [ A_transpose*Kt*A ]
	[ Refer Report/Manual for more Details ]

*/
tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd>material_routine(int elements, double gauss_weights, MatrixXd Lengths, MatrixXd Areas, MatrixXd  E_mat, MatrixXd u, double time_step, double eta, MatrixXd eps_pl, double yield, double E) {

	MatrixXd F_int= MatrixXd::Zero(elements + 1, 1);
	MatrixXd stress_mat = MatrixXd::Zero(elements , 1);
	MatrixXd Kt = MatrixXd::Zero(elements+1 , elements+1);
	MatrixXd epsilon = MatrixXd::Zero(elements , 1);
	MatrixXd u_temp(2, 1);

	/* Computing for each Element */
	for (int ele = 0; ele < elements; ele++) {
		
		/* u_temp : Temperory Local Displacement Matrix */
		u_temp(0, 0) = u(ele,0);
		u_temp(1, 0) = u(ele + 1, 0);
		
		/* Strain = B*u */
		MatrixXd t = B_mat(Lengths(ele,0)) * u_temp;
		epsilon(ele, 0) = t(0,0);

		double smt, ept, et;
		tie(smt, ept, et) = get_Stress_Strains(eps_pl(ele, 0), epsilon(ele, 0), eta, time_step / 100, time_step, yield, E);
		stress_mat(ele, 0) = smt;
		eps_pl(ele, 0) = ept;
		E_mat(ele, 0) = et;
		
		MatrixXd F_int_temp = B_mat(Lengths(ele, 0)).transpose() * Areas(ele, 0) * stress_mat(ele, 0) * (Lengths(ele, 0) / 2) * gauss_weights;
		MatrixXd Kt_temp = B_mat(Lengths(ele, 0)).transpose() * B_mat(Lengths(ele, 0)) * E_mat(ele, 0) * Areas(ele, 0) * (Lengths(ele, 0) / 2) * gauss_weights;

		/*Converting from Local to Global */

		F_int(ele, 0) = F_int(ele, 0) + F_int_temp(0, 0);
		F_int(ele + 1, 0) = F_int(ele + 1, 0) + F_int_temp(1, 0);

		MatrixXd A = assign(ele, elements);
		Kt = Kt + (A.transpose() * Kt_temp * A);

	}

	return { Kt, F_int, eps_pl, epsilon, stress_mat };

}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/