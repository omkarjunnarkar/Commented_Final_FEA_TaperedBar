/*
***********************************************************************************************************************************************************

												Finite Element Code for Tapered Bar Problem
															Plotting Routine

Author: Omkar Junnarkar, Semester-3 MSc. Computational Material Science
Matriculation Nr.: 66157	Email: omkar.junnarkar@student.tu-freiberg.de
IDE : Microsoft Visual Studio 2019

Objective : Use the GNUPLOT-IOSTREAM Interface, and plot Stress-Strain & Displacement-Time Curve

Libraries Included :

iostream: For Input-Output of the C++ Program
vector: To use vectors data type during file reading
src/rapidcsv: To read CSV files
gnuplot-iostream.h : Header-Only file to Link C++ iostream and GNUPLOT interface : Author Dan Stahlke [ Link in Main Program ]
plotme.h: Plotting Routine header file

*/
#include<iostream>
#include<vector>
#include"src/rapidcsv.h"
#include"gnuplot-iostream.h"
#include"plotme.h"

/*
	To reduce the effort of specifying libraries/class for each command/data type
*/

using namespace std;
using namespace rapidcsv;

/*
	Arguments: None
	Output: Plots (Stress-Strain, Displacement-Time)

	The files generated by FEM Code are read by this function, wrote into 'vector', and then Plotted.
	[ See Documentation on github for Details(link in main.cpp) ]
*/
int plot() {

	/* 'gp' an instance of Gnuplot class, and 'dy,dstress, dstrain, etc' are instances of Document class
		To use Gnuplot, as the argument of instance, it is required to give Path of executable file of Gnuplot.
	*/

	Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");
	Document dy("displacement_tapered_10e.csv");
	Document dstress("stress_10e.csv");
	Document dstrain("strain_10e.csv");
	vector<double> getdisplacement = dy.GetColumn<double>(0);
	vector<double> getstress = dstress.GetColumn<double>(0);
	vector<double> getstrain = dstrain.GetColumn<double>(0);

	/* Creating 2D vector of doubles, according to number of computed data points */
	int data_size = getstress.size();
	vector<vector<double>>stress_strain(data_size,vector<double>(2));

	for (int i = 0; i < data_size; i++) {

		stress_strain[i][0] = getstrain[i];
		stress_strain[i][1] = getstress[i];
		//cout << getstress[i] << " , " << getstrain[i] << endl;	
		
	}

	/* Passing GNUPLOT commands through GNUPLOT syntax, in form of String to the instance 'gp' */

	gp << "set multiplot layout 2,1 rowsfirst\n";
	gp << "set xlabel 'Strain' font 'Times - Roman, 10'\n";
	gp << "set ylabel 'Stress' font 'Times - Roman, 10'\n";
	gp << "set xtics font 'Arial, 7'\n";
	gp << "set ytics font 'Arial, 7'\n";

	gp << "set title 'Stress-Strain Curve'\n";
	gp << "plot '-' with lines title 'Stress'\n";

	gp.send1d(stress_strain);

	//gp << "set size 0.36, 0.5\n";
	gp << "set title 'Displacement Curve'\n";
	gp << "set xlabel 'Time-Steps' font 'Times - Roman, 10'\n";
	gp << "set ylabel 'Displacement' font 'Times - Roman, 10'\n";
	gp << "set xtics font 'Arial, 7'\n";
	gp << "set ytics font 'Arial, 7'\n";

	gp << "plot '-' with lines title 'u'\n";
	gp.send1d(getdisplacement);

	cin.get();

};