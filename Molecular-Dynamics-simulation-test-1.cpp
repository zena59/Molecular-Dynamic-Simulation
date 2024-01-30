#include<iostream>
#include<math.h>
#include<time.h>
//#include<Windows.h>
#include<fstream>
#include "Energy.h"
using namespace std;

int main()
{
	//Number of times simulation must run
	int num_times = 100;

	//time spent (do not change)
	float tim = 0;

	//delta t
	float dt = 0.01; //seconds

	//number of particle
	const int n = 100;

	//Length of the box
	const int L = 10; // sigma units (1 sigma = diameter of a molecules)

	//Volume of the box
	const float volume = float(L * L * L);

	//temperature of the system in kelvin
	float T = 30; //Kelvin 

	//boltzmen constant
	float kb = 1.380649; //H Joules ( H = 10^-23)

	//mass of the particle
	float m = 1; //unit is a mass of each molecule

	//distance between particle ( do not change)
	float r;

	//!!!DO NOT CHANGE BELOW VALUES!!!
	//force on particle i due to j and total energy
	float energy = 0;
	float energy_initial = 0;
	float forceij = 0;
	float potential = 0;
	float forcex = 0;
	float forcey = 0;
	float forcez = 0;


	//defining position of the particle
	float posx[n];
	float posy[n];
	float posz[n];

	//Distance between particle component vise
	float x_position;
	float y_position;
	float z_position;

	//for loops
	int i = 0;
	int j = 0;
	int k = 0;

	//packing fraction
	float pac_fra = (float(n) * 4.1887902 * 100) / volume;
	cout << "Packing percentage : " << pac_fra << " %" << endl;

	//density
	float density = (m * float(n)) / volume;

	//pressure
	float pressure = (float(n) * 8.3144626 * T) / volume;
	cout << "Ideal Pressure : " << pressure << endl;
	cout << "Density : " << density << endl;

	//Initializing position
	cout << "\nUpdating Initial position" << endl;

	//number of column or rows of thr grid
	int NOR = 0;
	while (true)
	{
		if ((NOR * NOR * NOR) <= n)
		{
			NOR++;
		}
		else
		{
			break;
		}
	}
	
	cout << "Dimention of atom/molecule placed : " << NOR << "X" << NOR << "X" << NOR << endl;

	// d is distance between particle
	float d = L / (float(NOR) + 2);
	
	//Now properly initializing position
	int l = 0;
	int m1 = 0;

	for (i = 0;i < n;i++)
	{
		posx[i] = (float(i % NOR) + 1) * d;

		if (j == NOR)
		{
			k++;
			j = 0;
		}

		if (m1 == (NOR * NOR))
		{
			l++;
			m1 = 0;
		}
		m1++;
		j++;
		posy[i] = (float(k % NOR) + 1) * d;
		posz[i] = (float(l % NOR) + 1) * d;
	}
	cout << "Initial position updated\n" << endl;


	//output x position for test before simulation
	ofstream test_pos("Updates/position_x.txt");
	for (i = 0;i < n;i++)
	{
		test_pos << posx[i] << endl;
	}
	test_pos.close();

	//output y position for test before simulation
	ofstream test_posy("Updates/position_y.txt");
	for (i = 0;i < n;i++)
	{
		test_posy << posy[i] << endl;
	}
	test_posy.close();

	//output z position for test before simulation
	ofstream test_posz("Updates/position_z.txt");
	for (i = 0;i < n;i++)
	{
		test_posz << posz[i] << endl;
	}
	test_posz.close();

	//defining velocity

	float velx[n];
	float vely[n];
	float velz[n];

	float ran = 0;
	float A = (kb * T) / m;

	cout << "Initializing velocity" << endl;

	srand(time(0));

	for (i = 0;i < n;i++)
	{
		ran = (rand() % 10) * 0.1;
		velx[i] = A * (ran - 0.5);
		
	}

	for (i = 0;i < n;i++)
	{
		ran = (rand() % 10) * 0.1;
		vely[i] = A * (ran - 0.5);
	}

	for (i = 0;i < n;i++)
	{
		ran = (rand() % 10) * 0.1;
		velz[i] = A * (ran - 0.5);
	}

	cout << "Initial velocity updated\n" << endl;

	//output x component of velocity for test before simulation
	ofstream test_vel("Updates/velocity_x.txt");
	for (i = 0;i < n;i++)
	{
		test_vel << velx[i] << endl;
	}
	test_vel.close();

	//output y component of velocity for test before simulation
	ofstream test_vely("Updates/velocity_y.txt");
	for (i = 0;i < n;i++)
	{
		test_vely << vely[i] << endl;
	}
	test_vely.close();


	//to calculate initial energy
	for (i = 0;i < n;i++)
	{

		for (j = 0;j < n;j++)
		{

			if (i == j)
			{
				continue;
			}
			else
			{
				if (fabsf(posx[i] - posx[j]) >= (L - 4))
				{
					x_position = L - fabsf(posx[i] - posx[j]);
				}

				else
				{
					x_position = fabsf(posx[i] - posx[j]);
				}


				if (fabsf(posy[i] - posy[j]) >= (L - 4))
				{
					y_position = L - fabsf(posy[i] - posy[j]);
				}

				else
				{
					y_position = fabsf(posy[i] - posy[j]);
				}


				if (fabsf(posz[i] - posz[j]) >= (L - 4))
				{
					z_position = L - fabsf(posz[i] - posz[j]);
				}

				else
				{
					z_position = fabsf(posz[i] - posz[j]);
				}


				r = sqrt(pow(x_position, 2) + pow(y_position, 2) + pow(z_position, 2));
				potential = potential + psi(r);

			}
		
		}

		energy_initial = energy_initial + (pow(velx[i], 2) + pow(vely[i], 2) + pow(velz[i], 2)) + potential;

	}
	//SIMULATION
	for (k = 0;k <= num_times;k++)
	{

		energy = 0;
		cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
		cout << "Progress :" << (float(k) / float(num_times)) * 100 << "%"; //progress bar

		for (i = 0;i < n;i++)
		{

			for (j = 0;j < n;j++)
			{

				if (i == j)
				{
					continue;
				}
				else
				{
					/*if (fabsf(posx[i] - posx[j]) >= (L - 3))
					{
						x_position = L - fabsf(posx[i] - posx[j]);
					}
					else
					{
						x_position = fabsf(posx[i] - posx[j]);
					}
					if (fabsf(posy[i] - posy[j]) >= (L - 3))
					{
						y_position = L - fabsf(posy[i] - posy[j]);
					}
					else
					{
						y_position = fabsf(posy[i] - posy[j]);
					}
					if (fabsf(posz[i] - posz[j]) >= (L - 3))
					{
						z_position = L - fabsf(posz[i] - posz[j]);
					}
					else
					{
						z_position = fabsf(posz[i] - posz[j]);
					}
					*/

					//reassigning values to zero before going to next particle
					x_position = 0;
					y_position = 0;
					z_position = 0;
					forcex = 0;
					forcey = 0;
					forcez = 0;
					forceij = 0;
					potential = 0;

					x_position = fabs(posx[i] - posx[j]);
					y_position = fabs(posy[i] - posy[j]);
					z_position = fabs(posz[i] - posz[j]);

					if (x_position >= (L - 3.2) && y_position >= (L - 3.2) && z_position >= (L - 3.2))
					{
						x_position = L - x_position;
						y_position = L - y_position;
						z_position = L - z_position;
					}

					//distance between ith and jth particle
					r = sqrt(pow(x_position, 2) + pow(y_position, 2) + pow(z_position, 2));

					if (r == 0) //check for failed simulation
					{
						cout << "\nParticles have overlapped!!\nSimulation failed!" << endl;
						exit(0);
					}

					potential = potential + psi(r);
					forceij = force(r) / r;

					//component of the forces
					forcex = forcex + (forceij * x_position);
					forcey = forcey + (forceij * y_position);
					forcez = forcez + (forceij * z_position);
				}
			}

			//updating new position
			posx[i] = posx[i] + (velx[i] * dt);
			posy[i] = posy[i] + (vely[i] * dt);
			posz[i] = posz[i] + (velz[i] * dt);

			//updating new velocity
			velx[i] = velx[i] + ((forcex * dt) / m);
			vely[i] = vely[i] + ((forcey * dt) / m);
			velz[i] = velz[i] + ((forcez * dt) / m);

			//constraining particles
			if (posx[i] < 0 && velx[i] < 0) //x-position changing. Bringing back to box
			{
				posx[i] = L;
			}

			if (posx[i] > L  && velx[i] > 0)
			{
				posx[i] = 0;
			}


			if (posy[i] < 0 && vely[i] < 0) //y-position changing
			{
				posy[i] = L;
			}

			if (posy[i] > L && vely[i] > 0)
			{
				posy[i] = 0;
			}


			if (posz[i] < 0 && velz[i] < 0) //z-position changing
			{
				posz[i] = L;
			}

			if (posz[i] > L && velz[i] > 0)
			{
				posz[i] = 0;
			}
			

			//energy calculation
			energy = energy + (pow(velx[i], 2) + pow(vely[i], 2) + pow(velz[i], 2)) + potential;

		}
		/*if ((energy_initial - energy) >= 1 && (energy_initial - energy) <= -1)
		{
			cout << "Total Energy of system has changed significantly" << endl;
			cout << "Simulation Failed" << endl;
			break;
		}
		*/
		tim = tim + dt;
	}

	cout << "\nInitital Energy : " << energy_initial << endl;
	cout << "Final Energy : " << energy << endl;
	cout << "Total time simulated : " << tim << endl;

	/*calculationg pressure
	float pre1 = 0;
	float pre2 = 0;
	pressure = 0;
	cout << "\nCalculating pressure" << endl;
	for (i = 0;i < n;i++)
	{
		pre1 = (velx[i] * velx[i]) + (vely[i] * vely[i]) + (velz[i] * velz[i]);
		pre2 = pre2 + pre1;
	}
	pressure = (pre2 * m * density) / 3;
	cout << "Pressure : " << pressure << endl;
	*/

	//output x position for test after simulation
	ofstream test_pos_final("Updates/position_x_final.txt");
	for (i = 0;i < n;i++)
	{
		test_pos_final << posx[i] << endl;
	}
	test_pos_final.close();

	//output y position for test after simulation
	ofstream test_posy_final("Updates/position_y_final.txt");
	for (i = 0;i < n;i++)
	{
		test_posy_final << posy[i] << endl;
	}
	test_posy_final.close();

	//output z position for test after simulation
	ofstream test_posz_final("Updates/position_z_final.txt");
	for (i = 0;i < n;i++)
	{
		test_posz_final << posz[i] << endl;
	}
	test_posz_final.close();

	//output x component of velocity for test after simulation
	ofstream test_vel_final("Updates/velocity_x_final.txt");
	for (i = 0;i < n;i++)
	{
		test_vel_final << velx[i] << endl;
	}
	test_vel_final.close();

	//output y component of velocity for test after simulation
	ofstream test_vely_final("Updates/velocity_y_final.txt");
	for (i = 0;i < n;i++)
	{
		test_vely_final << vely[i] << endl;
	}
	test_vely_final.close();

	cout << "Uploaded final position and velocity!" << endl;

	//system("pause");

	return 0;
}
