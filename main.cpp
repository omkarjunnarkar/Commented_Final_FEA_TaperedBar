/*
***********************************************************************************************************************************************************

												Finite Element Code for Tapered Bar Problem
																Main Program

Author: Omkar Junnarkar, Semester-3 MSc. Computational Material Science
Matriculation Nr.: 66157	Email: omkar.junnarkar@student.tu-freiberg.de
IDE : Microsoft Visual Studio 2019
													PERSONAL PROGRAMMING PROJECT

External Libraries/Resources:

Eigen : https://eigen.tuxfamily.org/dox/							Owner: TuxFamily
RapidCSV : https://github.com/d99kris/rapidcsv						Author: Kristofer Berggren
GNUPLOT-IOSTREAM : https://github.com/dstahlke/gnuplot-iostream		Author: Dan Stahlke
Microsoft VCPKG :  https://github.com/Microsoft/vcpkg				Owner: Microsoft

Objective:	 To perform Numerical simulation of a 1-D Truss problem, a Tapered Bar with Axial Force acting on the center and fixed on two ends. [ Refer Report/Manual ]
Input:		 Open 'UserInputFileTaperedBarFEA.csv' and Input the values mentioned in the file. Details/Symbols of the Input values/parameters provided below in the code.
Output:		 Plot of Stress vs. Strain & Displacement vs. Time Steps, CSV files containing Displacement of center node, Stress and Strain values of all elements,
			 Reaction Forces, Residual at each step, Areas/Lenghts after Meshing, Assignment Matrices of all Elements, B-Matrix for all ELements, Input Database File

Steps:	1. Make sure the environment is set for C++17 (Necessary for GNUPLOT Interface)
		2. Make sure the Configuration is set correctly. (Author Config.: Release x64)
		3. Install GNUPLOT from http://gnuplot.info/download.html
		4. Clone git rep https://github.com/Microsoft/vcpkg.git
		5. Open vcpkg directory from Terminal, Use Command: .\bootstrap-vcpkg.bat
		6. Install vcpkg boost from terminal, Use Command: vcpkg install boost:x64-windows boost:x86-windows
		7. Use vcpkg.exe to install & integrate, Use Commands: .\vcpkg.exe install boost:x64-windows boost:x86-windows ; .\vcpkg.exe integrate install
		8. Download Eigen and Rapid CSV from above mentioned links, Extract the Folders, Add in Properties > Additional Include Directories > Eigen,RapidCSV
		9. Edit the 'UserInputFileTaperedBar.csv' and run main.cpp.

*/
/*
*
iostream: For Input-Output of the C++ Program
Eigen/Dense: For Using Dense Matrix Interface of Linear Algebra Library Eigen
iomanip: To manipulate the number of decimal places in output
fstream: To create the stream of file and write
vector: To use vectors data type during file reading
src/rapidcsv: To read CSV files
material.h: Material Routine header file
element.h: Element Routine header file
plotme.h: Plotting Routine header file

*/

#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include<vector>
#include<src/rapidcsv.h>
#include"material.h"
#include"element.h"
#include"plotme.h"

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*
To reduce the effort of specifying libraries/class for each command/data type
*/

using namespace std;
using namespace Eigen;
using namespace rapidcsv;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*
Main Program, consists of: Newton Raphson Iterations, Boundary Condition Application, File I/O
Arguments: None

Variables :-

ForceMax : Maximum Force Applied(N)
Area_1 : csArea of left element(mm^2)
Area_2 : csArea of right element(mm^2)
Length : Length of bar(mm)
E : Youngs Modulus (MPA)
steps_num : Number of Load Steps (Will also define the strain-rate & thus, material behvaiour)
t_tot : Total time for which force is applied (s)
yield : Yield Strength (MPA)
eta : Viscosity of the Material(1/s)
elements : Number of elements for meshing

Files Created:

displacement_tapeted_10e.csv : Displacement of the Selected/Focus/Center Node at each load step
stress_10e.csv : Stress of Selected/Focus Element at each load step
strain_10e.csv : Strain of Selected/Focus Element at each load step
All_stress.csv : Stress of All Elements at each load step
All_strain.csv : Strain of All Element at each load step
Residual.csv : Residual at the last Newton Raphson Iteration at each load step
InputFileFEATaperedBar.txt : Input Database File
B-Matrix.txt : B-Matrix of all elements
Assignment-Matrix.txt : Assignment Matrix for Assembly of All Element Properties
Areas-Lengths-After-Meshing.txt : Areas and Lenghts of all elements after meshing
LoadingSteps.csv : Value of Forces at each Load Step
ReactionForces.csv : Reaction Forces at the End Points (BCs)

*/

void main() {

	cout << "Starting the Finite Element Program" << endl;
	
	/*
		ofstream : Class for the file stream from 'std' library
	*/

	ofstream mydisplacefile("displacement_tapered_10e.csv");
	ofstream mystressfile("stress_10e.csv");
	ofstream mystrainfile("strain_10e.csv");
	ofstream Allstrainfile("All_strain.csv");
	ofstream Allstressfile("All_stress.csv");
	ofstream Residualfile("Residual.csv");
	ofstream input("InputFileFEATaperedBar.txt");
	ofstream bmat("B-Matrix.txt");
	ofstream amat("Assignment-Matrix.txt");
	ofstream ar_len("Areas-Lengths-After-Meshing.txt");
	ofstream loadhist("LoadingSteps.csv");
	ofstream ReactionForces("ReactionForces.csv");

	/*
		Document : Class from 'rapdcsv' library for reading CSV Files, inp: Instance of Document class
		Reading the User Input File into a vector
	*/

	Document inpfile("UserInputFileTaperedBarFEA.csv");
	vector<double> user_input = inpfile.GetColumn<double>(0);

	/*
		Assigning Raed values to the variables
	*/
	double ForceMax = user_input[0];			
	double Area_1 = user_input[1];				
	double Area_2 = user_input[2];				
	double Length = user_input[3];				
	double E = user_input[4]*1e3;			/* Multiplied by 1e3 to convert from GPa to MPa */	
	int steps_num = 10000;					/* If needed, can be changed here */
	double steps = steps_num;
	double t_tot = user_input[5];				
	double yield = user_input[6];			
	int elements = user_input[8];				
	double force = 0;						/* Initializing Value of Force Applied */	

	/*
		gauss_weights : Weight for the Gauss Point
		Considered only 1 Gauss Point, thus the weight = 2
		In case of computing for 2 Gauss Points, change the gauss_weights value and
		manipulate the Element subroutine functions
	*/

	int gauss_weights = 2;
	double eta = user_input[7];		

	/*
		MatrixXd : Data Type from Eigen Library
		Meaning: Matrix - X - d : Matrix with 'X', i.e, expandable dimensions with 'd', i.e 'double' type elements value
		Size can be specified in Brackets, for eg.: MatrixXd Example(10,4);

		Note: In the Program, Instead of Array, Only Matrices have been used to keep avoid any errors of data type computations

		Areas	: List of Areas of Elements
		Lenghts	: List of Lengths of Elements

	*/

	MatrixXd Areas, Lengths;

	/*
		Writing the Input Database
	*/

	input << "Input Data for FEA Code: " << endl;
	input << "Max. Force (Newton) = " << ForceMax << endl;
	input << "Area of Bar 1 (mm^2) = " << Area_1 << endl;
	input << "Area of Bar 2 (mm^2) = " << Area_2 << endl;
	input << "Length of Bar (mm) = " << Length << endl;
	input << "Young's Modulus (GPa) = " << E << endl;
	input << "Yield Stress (MPa) = " << yield << endl;
	input << "Viscosity (1/s) = " << eta << endl;
	input << "Total Time (s) = " << t_tot << endl;
	input << "Time Steps = " << steps << endl;
	input << "No. of Elements = " << elements << endl;
	input << "Gauss Points, Weights = " << "1 , " << gauss_weights << endl;

	/*
		Value of time elapsed at each load step
	*/

	steps = t_tot / double(steps);

	/*
		get_areas_lengths : Function from Element Routine, Output Data Type: Tuple (class: std) : Areas, Lengths
		[ See Details in 'element.cpp' ]
		tie: function template (class: std) to recieve multple outputs from functions giving output of 'Tuple' type
	*/

	tie(Areas, Lengths) = get_areas_lengths(elements, Area_1, Area_2, Length);	

	//cout << "Out of Element Areas Lengths" << endl;

	input << "List of Areas after Seeding: \n" << Areas << endl;
	input << "List of Lengths after Seeding: \n" << Lengths << endl;

	/*
		Initialization of Quantities, All of type: Eigen Matrix:

		F_ext : External Force
		E_mat : Stiffness(initially=Young's Modulous) of all elements
		eps_pl : Platsic Strain
		epsilon : Total Strain
		u : Global Displacement Matrix
		Kt : Global Stifness Matrix
		F_int : Global Internal Force
		strain : Strain of considered element
		stress : Stress of considered element

	*/

	MatrixXd F_ext = MatrixXd::Zero(elements + 1,1);
	MatrixXd E_mat= E*MatrixXd::Ones(elements , 1);
	MatrixXd eps_pl = MatrixXd::Zero(elements , 1);
	MatrixXd epsilon = MatrixXd::Zero(elements , 1);
	MatrixXd u = MatrixXd::Zero(elements + 1, 1);
	MatrixXd Kt, F_int, strain, stress;
	MatrixXd temp_mat1, temp_mat2;
	
	int count = 0;			/* Iteration Count for Newton Raphson*/
	
	/*
		Writing B-matrix,Assignment Matrix and Areas/Lenghts to file
	*/

	for (int e = 0; e < elements; e++) {
		bmat << "B - Matrix for Element [" << e << "] of length = " << Lengths(e, 0) << " mm is : \n" << B_mat(Lengths(e, 0)) << endl;
		amat << "Assignment Matrix for Element [" << e << "] in system of " << elements << " Elements is :\n" << assign(e, elements) << endl;
		ar_len << "Area of Element [" << e << "] = " << Areas(e, 0) << " Length of Element [" << e << "] = " << Lengths(e, 0) << endl;
	}
	bmat.close();
	amat.close();
	ar_len.close();

	/*
		Computation Strategy :

		Outer Loop : Loading Steps
		->Inner Loop : Newton Raphson Iteration
		->->Innermost Loops : Through Element Properties (needed for Boundary Condition Application
	*/

	while (abs(force) <= abs(ForceMax) ) {
				
		count++;
		
		//cout << "Force = " << force << endl;

		/*
		Replace the below line by following lines for Testing Variable Boundary Conditions
		[ Refer Report/Manual for Details regarding Test Cases & Variable Boundary Conditions ]
		For ORIGINAL & CASE-1 : F_ext(int(elements / 2), 0) = force;	Force at Center
		For CASE-2 & CASE-3 :	F_ext(elements , 0) = force;			Force at Right End
		For CASE-4 :			F_ext(elements , 0) = -force;			Compressive Force at Right End
		For CASE-5 :			F_ext(int(elements/2) , 0) = -force;	Compressive force at Center
		*/

		/* Providing the value of force in current step to the Global External Force Matrix */

		F_ext(int(elements / 2), 0) = force;								
		loadhist << force << endl;

		/*
			Entering the Newton - Raphson Method
			NR : Newton raphson Iteration Number

			Newton Raphson Iteration used to solve a Non-Linear Problem with the goal of Reduction of the Residual
			Lin f(x)|xi = f(xi) + f'(xi)*(x-xi)

			Strategy:

			->Initial Guess of u -> Call Material Routine and determine Stress/Strain/Plastic Strain/ Tangent Stiffness
			->Apply BCs -> Compute Residual = F_internal-F_external -> Compute change in u -> Update u -> Use it for next NR iteration

			To check the Convergence of Newton Raphson, on the last Iteration, values of Residual are witten to the Residual.csv
			General Criterion Considered for Convergence : ||Delta u|| < 0.005 * ||u||
		*/

		for (int NR = 0; NR < 6; NR++) {

			/*
				Calling Material Routine
				Output : Global Tangent Stiffness Matrix, Global Internal Force Matrix, Plastic Strain, Total Strain & Stress Values of All Elements
				Input : No. of elements, Gauss weights, Lengths, Areas, List of Stiffness, current displacement, no. of loading steps,
				viscosity, initial plastic strain, yield stress, Young's Modulus

				[ See 'material.cpp' for Details ]
			*/

			tie(Kt, F_int, eps_pl, strain, stress) = material_routine(elements, gauss_weights, Lengths, Areas, E_mat, u, steps, eta, eps_pl, yield, E);
			
			/*
			Replace the below line by following lines for Testing Variable Boundary Conditions
			[ Refer Report/Manual for Details regarding Test Cases & Variable Boundary Conditions]

			ORIGINAL:
			MatrixXd Kt_red = Kt.block(1, 1, Kt.cols() - 2, Kt.rows() - 2);
			MatrixXd G_red(elements - 1, 1),u_red(elements-1,1);
			for (int c = 1; c < elements; c++) { G_red(c - 1, 0) = G(c, 0); }
			for (int c = 1; c < elements; c++) { u_red(c - 1, 0) = u(c, 0); }
			for (int c = 1; c < elements; c++) { u(c, 0) = u_red(c - 1, 0); }

			CASE-1 & CASE-2 & CASE-3 & CASE-4 & CASE-5 :
			MatrixXd Kt_red = Kt.block(1, 1, Kt.cols() - 1, Kt.rows() - 1);
			MatrixXd G_red(elements , 1),u_red(elements,1);
			for (int c = 1; c < elements+1; c++) { G_red(c - 1, 0) = G(c, 0); }
			for (int c = 1; c < elements+1; c++) { u_red(c - 1, 0) = u(c, 0); }
			for (int c = 1; c < elements+1; c++) { u(c, 0) = u_red(c - 1, 0); }
			*/

			/*
				Block Method : To extract Matrix of Specific size from another Matrix, Library: Eigen
				To apply Boundary Conditions, Reduced Matrix is Required.

				Kt_red : Reduced Tangent Stiffness Matrix
				inv_Kt_red : Inverse of Reduced Kt
				G: Matrix of Residuals
				G_red : Reduced G Matrix
				u_red : Reduced displacement matrix
				del_u_red : Change in displacement(delta u)

				Block Function arguements : (Start RowIndex, Start ColIndex, End RowIndex, End ColIndex)
			*/

			MatrixXd Kt_red = Kt.block(1, 1, Kt.cols() - 2, Kt.rows() - 2);	
			MatrixXd inv_Kt_red = Kt_red.inverse();
			MatrixXd G = F_int - F_ext;
			MatrixXd G_red(elements - 1, 1),u_red(elements-1,1);					
			
			for (int c = 1; c < elements; c++) { G_red(c - 1, 0) = G(c, 0); }
			for (int c = 1; c < elements; c++) { u_red(c - 1, 0) = u(c, 0); }

			MatrixXd del_u_red = inv_Kt_red * G_red;								
			u_red = u_red - del_u_red;

			for (int c = 1; c < elements; c++) { u(c, 0) = u_red(c - 1, 0); }

			if (NR == 5) {
				Residualfile << G_red.norm() << endl;
				ReactionForces << F_int(0, 0) << " , " << F_int(elements, 0) << endl;

				/*if (del_u_red.norm() < 0.005 * u_red.norm()) {
					cout << "NR Converged" << endl;
				}
				else cout << "NR NOT Converged" << endl;*/
			}


		}
		
		/*
			Appending values of Center Node / First Element to the csv files
			By changing the Row Number, Details of other elements can also be extracted
		*/

		mydisplacefile << u((elements / 2), 0) << endl;
		mystressfile << stress((elements / 2) - 1, 0) << endl;
		mystrainfile << strain((elements / 2) - 1, 0) << endl;

		temp_mat1 = strain.reshaped(1, elements);
		temp_mat2 = stress.reshaped(1, elements);
		for (int fc = 0; fc < elements; fc++) {
			Allstrainfile << temp_mat1(0, fc) << ",";
			Allstressfile << temp_mat2(0, fc) << ",";
		}
		Allstrainfile << endl;
		Allstressfile << endl;

		/* Incrementing the force value for next loading step*/

		force = force + (1 / double(steps_num)) * ForceMax;
	}

	mydisplacefile.close();
	mystressfile.close();
	mystrainfile.close();
	Allstrainfile.close();
	Allstressfile.close();
	Residualfile.close();
	input.close();
	loadhist.close();
	ReactionForces.close();

	cout << "Computation Done ! " << endl;
	cout << "Activating GNUPLOT using vcpkg and gnupuplot-iostream..\n";

	/*
	Calling the plot() function from 'plotme.h'
	Arguments: None
	[ See plotme.cpp for Details ]
	*/

	plot();
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
