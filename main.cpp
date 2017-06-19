#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <string>
using namespace std;

int main() {
	cout << "Started... " << endl;

	ifstream readZmierzonePolozenie("measuredPosition.txt");
	ifstream readZmierzonaPredkosc("measuredVelocity.txt");

	// check for file opening errors
	if (!readZmierzonePolozenie.is_open() || !readZmierzonaPredkosc.is_open()) {
		std::cerr << "There was an error while opening the file.\n";
		return -1;
	}


	double Observations[2000][2];

	for (size_t i = 0; i < 2000; ++i) {
		readZmierzonePolozenie >> Observations[i][0];
		readZmierzonaPredkosc >> Observations[i][1];
	}

	readZmierzonePolozenie.close();
	readZmierzonaPredkosc.close();


	// GIVEN
	double V0x_velocity = { 0.00137918195926662657 }; 
	double x0_position = { 0.0002550105561488041060 };

	// INITIAL CONDITIONS
	double ax = { 0 }; // accelaration
	double delta_t = { 0.01 }; // seconds [s]  -->  100 Hz means 100 cycles per second

	// PROCESS ERRORS IN PROCESS COVARIANCE MATRIX
	double delta_Px = { 0.5 }; // meters [m]
	double delta_Pvx = { 10.0 }; // meters per second [m/s]

	// OBSERVATIONS ERRORS
	double delta_X = { 0.5 }; // meters [m]
	double Vx = { 10.0 }; // meters per second [m/s]


	//-----------------------------------------------------------||
	//-----------------------THE SOLUTION------------------------||
	//-----------------------------------------------------------||

	// LEGEND
	// for example A|x| means, that A is a normal variable and x is subscript

	// INITIAL VALUES
	double Xk1[2][1] = { // Xk1 means X|k-1|
		{ x0_position },
		{ V0x_velocity }
	};
	double Pk1[2][2] = { // Pk1 means P|k-1|
		{ pow(delta_Px, 2), 0.0 },
		{ 0.0, pow(delta_Pvx, 2) }
	};

	ofstream saveResultsVelocity("KalmanFilterVelocity.txt");
	ofstream saveResultsPosition("KalmanFilterPosition.txt");
	saveResultsVelocity.setf(ios_base::fixed, ios_base::floatfield);
	saveResultsPosition.setf(ios_base::fixed, ios_base::floatfield);

	// KALMAN FILTER LOOP
	for (size_t i = 0; i < 2000; ++i) {
		// 1) FIRST STEP IS TO CALCULATE THE PREDICTED STATE
		// GENERAL FORMULA  -->  X|xp| = A * X|k-1| + B * u|k| + w|k|

		// firstly I calculate A * X|k-1|
		double A_matrix[2][2] = {
			{ 1.0, delta_t },
			{ 0.0, 1.0 }
		};

		double result_A_matrix_times_Xk1[2][1] = {
			{ A_matrix[0][0] * Xk1[0][0] + A_matrix[0][1] * Xk1[1][0] },
			{ A_matrix[1][0] * Xk1[0][0] + A_matrix[1][1] * Xk1[1][0] },
		};

		// secondly I calculate B * u|k|
		double B_matrix[2][1] = {
			{ 0.5 * pow(delta_t, 2) },
			{ delta_t }
		};

		double uk[1][1] = { // uk means u|k|
			{ ax }
		};

		double result_B_matrix_times_uk[2][1] = {
			{ B_matrix[0][0] * uk[0][0] },
			{ B_matrix[1][0] * uk[0][0] }
		};

		// this value isn not taken into account
		double wk[1][1] = { // wk means w|k|
			{ 0.0 }
		};

		// next step is to add two results: A * X|k-1|  and  B * u|k|
		double Xkp[2][1] = { // Xkp means X|kp|
			{ result_A_matrix_times_Xk1[0][0] + result_B_matrix_times_uk[0][0] },
			{ result_A_matrix_times_Xk1[1][0] + result_B_matrix_times_uk[1][0] }
		};


		// ------------------------------------------------------------------------------


		// 2) SECOND STEP IS "THE INITIAL PROCESS COVARIANCE MATRIX"
		// THESE VALUES ARE UPDATED IN LAST STEP


		// ------------------------------------------------------------------------------


		// 3) THIRD STEP IS "THE PREDICTED PROCESS COVARIANCE MATRIX"
		// GENERAL FORMULA  -->  P|kp| = A * P|k-1| * A ^T + Q|k|

		// firstly I calculate A * P|k-1|
		// A_matrix is already initialized in step 1
		// P|k-1| is already initialized in step 2
		double result_A_matrix_times_Pk1[2][2] = { // Pk1 means in this case P|k-1|
			{ A_matrix[0][0] * Pk1[0][0] + A_matrix[0][1] * Pk1[1][0], A_matrix[0][0] * Pk1[0][1] + A_matrix[0][1] * Pk1[1][1] },
			{ A_matrix[1][0] * Pk1[0][0] + A_matrix[1][1] * Pk1[1][0], A_matrix[1][0] * Pk1[0][1] + A_matrix[1][1] * Pk1[1][1] },
		};

		// secondly I am writing how A ^T matrix look like
		double A_matrix_transposed[2][2] = {
			{ 1, 0 },
			{ delta_t, 1 }
		};

		// the last step is to multiply these two arrays
		double Pkp[2][2] = {
			{ result_A_matrix_times_Pk1[0][0] * A_matrix_transposed[0][0] + result_A_matrix_times_Pk1[0][1] * A_matrix_transposed[1][0], result_A_matrix_times_Pk1[0][0] * A_matrix_transposed[0][1] + result_A_matrix_times_Pk1[0][1] * A_matrix_transposed[1][1] },
			{ result_A_matrix_times_Pk1[1][0] * A_matrix_transposed[0][0] + result_A_matrix_times_Pk1[1][1] * A_matrix_transposed[1][0], result_A_matrix_times_Pk1[1][0] * A_matrix_transposed[0][1] + result_A_matrix_times_Pk1[1][1] * A_matrix_transposed[1][1] },
		};


		// change two values on zeros
		Pkp[0][1] = 0.0;
		Pkp[1][0] = 0.0;


		// ------------------------------------------------------------------------------


		// 4) FOURTH STEP IS TO CALCULATE THE KALMAN GAIN
		// GENERAL FORMULA  -->  K = ( P|kp| * H ^T ) / ( H * P|kp| * H ^T + R )

		// firstly I calculate P|kp| * H ^T
		// P|kp| is already in step 3
		// H ^T matrix is used to change the format P|kp| on K (Kalman gain)
		double H_matrix_transposed[2][2] = {
			{ 1, 0 },
			{ 0, 1 }
		};

		double Pkp_times_H_matrix_transposed[2][2] = {
			{ Pkp[0][0] * H_matrix_transposed[0][0] + Pkp[0][1] * H_matrix_transposed[1][0], Pkp[0][0] * H_matrix_transposed[0][1] + Pkp[0][1] * H_matrix_transposed[1][1] },
			{ Pkp[1][0] * H_matrix_transposed[0][0] + Pkp[1][1] * H_matrix_transposed[1][0], Pkp[1][0] * H_matrix_transposed[0][1] + Pkp[1][1] * H_matrix_transposed[1][1] },

		};


		// secondly I calculate H * P|kp|
		double H_matrix[2][2] = {
			{ 1, 0 },
			{ 0, 1 }
		};

		double H_matrix_times_Pkp[2][2] = {
			{ H_matrix[0][0] * Pkp[0][0] + H_matrix[0][1] * Pkp[1][0], H_matrix[0][0] * Pkp[0][1] + H_matrix[0][1] * Pkp[1][1] },
			{ H_matrix[1][0] * Pkp[0][0] + H_matrix[1][1] * Pkp[1][0], H_matrix[1][0] * Pkp[0][1] + H_matrix[1][1] * Pkp[1][1] },
		};

		// thirdly I calculate H_matrix_times_Pkp times H_matrix_transposed
		double H_matrix_times_Pkp_times_H_transposed[2][2] = {
			{ H_matrix_times_Pkp[0][0] * H_matrix_transposed[0][0] + H_matrix_times_Pkp[0][1] * H_matrix_transposed[1][0], H_matrix_times_Pkp[0][0] * H_matrix_transposed[0][1] + H_matrix_times_Pkp[0][1] * H_matrix_transposed[1][1] },
			{ H_matrix_times_Pkp[1][0] * H_matrix_transposed[0][0] + H_matrix_times_Pkp[1][1] * H_matrix_transposed[1][0], H_matrix_times_Pkp[1][0] * H_matrix_transposed[0][1] + H_matrix_times_Pkp[1][1] * H_matrix_transposed[1][1] },
		};


		// fourth step is H_matrix_times_Pkp_times_H_transposed added to R
		double R_matrix[2][2] = {
			{ pow(delta_X, 2), 0 },
			{ 0, pow(Vx, 2) }
		};

		double Kalman_gain_counter[2][2] = {
			{ H_matrix_times_Pkp_times_H_transposed[0][0] + R_matrix[0][0], H_matrix_times_Pkp_times_H_transposed[0][1] + R_matrix[0][1] },
			{ H_matrix_times_Pkp_times_H_transposed[1][0] + R_matrix[1][0], H_matrix_times_Pkp_times_H_transposed[1][1] + R_matrix[1][1] },
		};

		// the last step is to calculate Kalman gain
		double Kalman_gain[2][2] = {
			{ Pkp_times_H_matrix_transposed[0][0] / Kalman_gain_counter[0][0], Pkp_times_H_matrix_transposed[0][1] / Kalman_gain_counter[0][1] },
			{ Pkp_times_H_matrix_transposed[1][0] / Kalman_gain_counter[1][0], Pkp_times_H_matrix_transposed[1][1] / Kalman_gain_counter[1][1] },
		};

		// I have to correct values, because there might be happened that I divide 0 by 0
		if (Pkp_times_H_matrix_transposed[0][0] == 0 && Kalman_gain_counter[0][0] == 0)
			Kalman_gain[0][0] = 0;
		if (Pkp_times_H_matrix_transposed[0][1] == 0 && Kalman_gain_counter[0][1] == 0)
			Kalman_gain[0][1] = 0;
		if (Pkp_times_H_matrix_transposed[1][0] == 0 && Kalman_gain_counter[1][0] == 0)
			Kalman_gain[1][0] = 0;
		if (Pkp_times_H_matrix_transposed[1][1] == 0 && Kalman_gain_counter[1][1] == 0)
			Kalman_gain[1][1] = 0;


		// ------------------------------------------------------------------------------


		// 5) FIFTH STEP IS TO CALCULATE Y|k|  -->  THE NEW OBSERVATION
		// GENERAL FORMULA  --> Y|k| = C * Y|km| + Z|k|

		// I calculate C times Y|km|
		double C_matrix[2][2] = {
			{ 1, 0 },
			{ 0, 1 }
		};

		double Yk_matrix[2][1] = {
			{ C_matrix[0][0] * Observations[i][0] + C_matrix[0][1] * Observations[i][1] },
			{ C_matrix[1][0] * Observations[i][0] + C_matrix[1][1] * Observations[i][1] }
		};

		// this value isn not taken into account
		double Zk[1][1] = { // Zk means Z|k|
			{ 0.0 }
		};


		// ------------------------------------------------------------------------------


		// 6) SIXTH STEP IS TO CALCULATE THE CURRENT STATE
		// GENERAL FORMULA  --> X|k| = X|kp| + K * ( Y|k| - H * X|kp| )

		// firstly I calculate H * X|kp|
		double H_times_Xkp[2][1] = {
			{ H_matrix_transposed[0][0] * Xkp[0][0] + H_matrix_transposed[0][1] * Xkp[1][0] },
			{ H_matrix_transposed[1][0] * Xkp[0][0] + H_matrix_transposed[1][1] * Xkp[1][0] }
		};

		// secondly I calculate Y|k| - H_times_Xkp
		double Yk_minus_H_times_Xkp[2][1] = {
			{ Yk_matrix[0][0] - H_times_Xkp[0][0] },
			{ Yk_matrix[1][0] - H_times_Xkp[1][0] }
		};

		// thirdly I calculate K * Yk_minus_H_times_Xkp
		double K_times_Yk_minus_H_times_Xkp[2][1] = {
			{ Kalman_gain[0][0] * Yk_minus_H_times_Xkp[0][0] + Kalman_gain[0][1] * Yk_minus_H_times_Xkp[1][0] },
			{ Kalman_gain[1][0] * Yk_minus_H_times_Xkp[0][0] + Kalman_gain[1][1] * Yk_minus_H_times_Xkp[1][0] }
		};

		// the last step is to calculate the current state
		double Xk_current_state[2][1] = {
			{ Xkp[0][0] + K_times_Yk_minus_H_times_Xkp[0][0] },
			{ Xkp[0][1] + K_times_Yk_minus_H_times_Xkp[0][1] },
		};

		// Kalman gain 
		//cout << "6) The current state: " << endl;
		//cout << "| " << Xk_current_state[0][0] << " |" << endl;
		//cout << "| " << Xk_current_state[1][0] << " |" << endl;

		saveResultsVelocity << setprecision(18) << Xk_current_state[1][0] << endl;
		saveResultsPosition << setprecision(18) << Xk_current_state[0][0] << endl;

		// ------------------------------------------------------------------------------


		// 7) UPDATING THE PROCESS COVARIANCE MATRIX
		// GENERAL FORMULA  --> P|k| = ( I - K * H ) * P|kp|

		// firstly I calculate K * H
		double Kalman_gain_times_H_matrix[2][2] = {
			{ Kalman_gain[0][0] * H_matrix[0][0] + Kalman_gain[0][1] * H_matrix[1][0], Kalman_gain[0][0] * H_matrix[0][1] + Kalman_gain[0][1] * H_matrix[1][1] },
			{ Kalman_gain[1][0] * H_matrix[0][0] + Kalman_gain[1][1] * H_matrix[1][0], Kalman_gain[1][0] * H_matrix[0][1] + Kalman_gain[1][1] * H_matrix[1][1] }
		};

		// secondly I calculate I - Kalman_gain_times_H_matrix
		double I_matrix[2][2] = {
			{ 1, 0 },
			{ 0, 1 }
		};

		double I_minus_Kalman_gain_times_H_matrix[2][2] = {
			{ I_matrix[0][0] - Kalman_gain_times_H_matrix[0][0], I_matrix[0][1] - Kalman_gain_times_H_matrix[0][1] },
			{ I_matrix[1][0] - Kalman_gain_times_H_matrix[1][0], I_matrix[1][1] - Kalman_gain_times_H_matrix[1][1] }
		};

		// lastly I calculate I_minus_Kalman_gain_times_H_matrix times P|kp|
		double process_covariance_matrix[2][2] = {
			{ I_minus_Kalman_gain_times_H_matrix[0][0] * Pkp[0][0] + I_minus_Kalman_gain_times_H_matrix[0][1] * Pkp[1][0], I_minus_Kalman_gain_times_H_matrix[0][0] * Pkp[0][1] + I_minus_Kalman_gain_times_H_matrix[0][1] * Pkp[1][1] },
			{ I_minus_Kalman_gain_times_H_matrix[1][0] * Pkp[0][0] + I_minus_Kalman_gain_times_H_matrix[1][1] * Pkp[1][0], I_minus_Kalman_gain_times_H_matrix[1][0] * Pkp[0][1] + I_minus_Kalman_gain_times_H_matrix[1][1] * Pkp[1][1] }
		};


		// ------------------------------------------------------------------------------


		// 8) UPDATE VARIABLES

		x0_position = Xkp[0][0];
		V0x_velocity = Xkp[1][0];

		Xk1[0][0] = Xk_current_state[0][0];
		Xk1[1][0] = Xk_current_state[1][0];

		Pk1[0][0] = process_covariance_matrix[0][0];
		Pk1[0][1] = process_covariance_matrix[0][1];
		Pk1[1][0] = process_covariance_matrix[1][0];
		Pk1[1][1] = process_covariance_matrix[1][1];
	}

	saveResultsVelocity.close();
	saveResultsPosition.close();


	// ------------------------------------------------------------------------------

	// replace . on , (Excel requires it)
	{
		ifstream readFile("KalmanFilterPosition.txt");
		ofstream writeFile("KalmanFilterPositionCorrectness.txt");

		// check for file opening errors
		if (!readFile.is_open() || !writeFile.is_open()) {
			std::cerr << "There was an error while opening the file.\n";
			return -1;
		}

		string temp;
		for (size_t i = 0; i < 2000; ++i) {
			readFile >> temp; // read from file

			replace_if(temp.begin(), temp.end(), [](char c) { return c == '.'; }, ',');

			writeFile << temp << endl;
		}

		readFile.close();
		writeFile.close();
	}

	cout << endl << "Finished." << endl;

	cin.get();
}
