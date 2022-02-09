/*
***********************************************************************************************************************************************************

												Finite Element Code for Tapered Bar Problem
														Header File Element Routine

Author: Omkar Junnarkar, Semester-3 MSc. Computational Material Science
Matriculation Nr.: 66157	Email: omkar.junnarkar@student.tu-freiberg.de
IDE : Microsoft Visual Studio 2019
*/

/*
	To reduce the effort of specifying libraries/class for each command/data type
*/
using namespace std;
using namespace Eigen;

MatrixXd B_mat(double le);
int signum(double x);
MatrixXd assign(int ele_index_cpp, int numele);
tuple<MatrixXd, MatrixXd>get_areas_lengths(int elements, double a_1, double a_2, double len);
