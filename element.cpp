/*
***********************************************************************************************************************************************************

												Finite Element Code for Tapered Bar Problem
															Element Routine

Author: Omkar Junnarkar, Semester-3 MSc. Computational Material Science
Matriculation Nr.: 66157	Email: omkar.junnarkar@student.tu-freiberg.de
IDE : Microsoft Visual Studio 2019

Objective : Handle the Element Level Operations for the Geometry considered in the Finite Element Problem

Libraries Included :

iostream: For Input-Output of the C++ Program
Eigen/Dense: For Using Dense Matrix Interface of Linear Algebra Library Eigen
iomanip: To manipulate the number of decimal places in output
fstream: To create the stream of file and write
element.h: Element Routine header file
math : For mathematical numerical operations

*/

#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include<math.h>
#include"element.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*
	To reduce the effort of specifying libraries/class for each command/data type
*/

using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*
	Arguments : Element Length, type : double
	Ouput : B-Matrix of format [ -1/L  1/L ] , type : Eigen Matrix

	B-Matrix is derived from the Partial Differentiation of Shape Function w.r.t Co-ordinate length ( dN/dX = dN/dZi * dZi/dX )
	[ Refer Report/Manual for more Details ]

*/
MatrixXd B_mat(double le) {
	MatrixXd b(1, 2);
	b(0, 0) = -1 / le;
	b(0, 1) = 1 / le;
	return b;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*
	Arguments: any number, type: double
	Ouput : 'Sign' of the number

	If input is a Non-Negetive number, returns +1
	else returns -1

	Replicates Signum Function
*/
int signum(double x) {
	if (x >= 0) return 1;
	else return -1;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*
	Arguments: Index of the Element According to C++, Total Number of Elements
	Output : The assignment matrix for each element

	Assignment Matrix is Used to map local element matrices to global geomerty .
	[ Refer Report/Manual for more Details ]

	Note : In C++, the indexing starts with '0'. Thus, the Nth element physically, would be (N-1)th element numerically.
*/
MatrixXd assign(int ele_index_cpp, int numele) {

	//ele_index is current element number, not defined here!! different for each ele !!

	MatrixXd a_element = MatrixXd::Zero(2, numele + 1);
	a_element(0, ele_index_cpp) = 1;
	a_element(1, ele_index_cpp + 1) = 1;
	return a_element;
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/


/*
	tuple : template from 'std' class. Ties and Passes values together ; Used to output multiple values at once

	Arguments :
			elements : Number of elements
			ar_1 : Area of Left End
			ar_2 : Area of Right Right
			len	 : Length of Bar

	Output :
			After Meshing -
			List of Areas, type Eigen Matrix
			List of Lengths, type Eigen Matrix

	For an Approximation of a Tapered Bar, Elements with fractional change in Area are considered.
	Approach : 
	-> Area of Left End and Right End of the bar are given -> Assuming the bar having Circular Cross Section, Radius of both Ends can be computed
	-> An approximate angle of inclination "theta" can be computed: theta = (Radius2-Radius1) / Length 
	-> Increments of radius for each progressing element can be found on the basis of previous element radius -> Area for each element is computed 

	len_per_ele : Length of each element
	lengths : List of Lengths
	r1,r2 : Radius of left & right End of Bar
	theta : Approximated Angle

*/
tuple<MatrixXd, MatrixXd>get_areas_lengths(int elements, double ar_1, double ar_2, double len) {
	MatrixXd areas(elements, 1);
	double len_per_ele = len / elements;
	MatrixXd lengths = len_per_ele*MatrixXd::Ones(elements, 1);
	double r1 = pow(ar_1 / 3.14, 0.5), r2 = pow(ar_2 / 3.14, 0.5);
	double theta = (r2 - r1) / len, y_n;

	/* Fixing the Areas of left and Right End of the Bar */
	areas(0, 0) = ar_1;
	areas(elements - 1, 0) = ar_2;
	
	/* Computing the areas of other elements */
	for (int i = 0; i < elements-2; i++) {

		/* Next Radius = First Radius + Angle*Length Progressed */

		y_n = r1 + theta * (i+1) * len_per_ele;
		areas(i+1, 0) = 3.14*pow(y_n,2);
	}

	return { areas,lengths };
}



