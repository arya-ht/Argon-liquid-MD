#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#pragma warning(disable:4996)


/*Joel C
 *Malena Fassnach
 *Arya HajiTaheri
 *Nico Marioni
 */

// Definitions
#define N 1000
#define dt 0.005
#define endT 10  // we go from time=0 to time=1
#define cutoff 4


//Global variables
double frameSize = 10.0;
double temp = 1.0;
double PEnergy = 0.0;
double KEnergy = 0.0;
double TEnergy = 0.0;

struct atom {
	double
		Fx, Fy, Fz,
		x, y, z,
		xo, yo, zo,
		vx, vy, vz;
};
struct atom a[N];


void initialize();
void force();
void integrate();

int main(void)
{
	// Files
	FILE *fpos = fopen("MDPosition.xyz", "w");


	double i;
	int j;

	initialize();

	printf("**********************************************************\n");
	printf("\t\t\tNumber of particles: %d\n\n", N);
	printf("\t\t\tL Box:  %.1f\n\n", frameSize);
	printf("(Reduced Units:positions, velocities)\n");
	printf("**********************************************************\n\n\n");

	for (i = 0; i <= endT; i += dt)
	{
		PEnergy = 0, KEnergy = 0;
		force();
		integrate();
		TEnergy = PEnergy + KEnergy;
	/*	printf("Total Energy = %lf\t\tTemp = %lf\n", TEnergy, temp);*/
		fprintf(fpos, "1000\nfilename: position.xyz\n");
		for (j = 0; j<N; j++)
		{
			fprintf(fpos, "ATOM%04d%12.5f%12.5f%12.5f\n", j, a[j].x, a[j].y, a[j].z);
		}

	}

	fclose(fpos);
	return 0;
}

void initialize()
{
	srand(time(NULL));
	int counter = 0, i, j, k;
	double sumvx = 0.0, sumvy = 0.0, sumvz = 0.0;
	double sumvx2 = 0.0, sumvy2 = 0.0, sumvz2 = 0.0;


	for (i = -(int)frameSize / 2; i<(int)frameSize / 2; i++)
	{
		for (j = -(int)frameSize / 2; j<(int)frameSize / 2; j++)
		{
			for (k = -(int)frameSize / 2; k<(int)frameSize / 2; k++)
			{
				//printf("%d %d %d\n",i, j, k);
				// Positions
				a[counter].z = k;
				a[counter].y = j;
				a[counter].x = i;
				// Velocities between  -1, 1.
				a[counter].vx = (double)rand() / RAND_MAX * 2.0 - 1.0;
				a[counter].vy = (double)rand() / RAND_MAX * 2.0 - 1.0;
				a[counter].vz = (double)rand() / RAND_MAX * 2.0 - 1.0;
				// Finding the c.g. Motion
				sumvx += a[counter].vx;
				sumvy += a[counter].vy;
				sumvz += a[counter].vz;
				// Finding v2 of motion
				sumvx2 += pow(a[counter].vx, 2);
				sumvy2 += pow(a[counter].vy, 2);
				sumvz2 += pow(a[counter].vz, 2);

				counter++;
			}
		}

	} // end of the initial settings. Time to start correcting for a specific temperature.



	  // Findin’ velocity c.g. And mean-squared velocity
	sumvx = sumvx / N; sumvy = sumvy / N; sumvz = sumvz / N;
	sumvx2 = sumvx2 / N; sumvy2 = sumvy2 / N; sumvz2 = sumvz2 / N;
	double correctVx = sqrt(3 * 1 / sumvx2);
	double correctVy = sqrt(3 * 1 / sumvy2);
	double correctVz = sqrt(3 * 1 / sumvz2);

	// Correcting the individual velocities to kill net movement
	for (i = 0; i<N; i++)
	{
		a[i].vx = (a[i].vx - sumvx)*correctVx;
		a[i].vy = (a[i].vy - sumvy)*correctVy;
		a[i].vz = (a[i].vz - sumvz)*correctVz;
		// Determining the previous position (Used in Motion fxn)
		a[i].xo = a[i].x - a[i].vx * dt;
		a[i].yo = a[i].y - a[i].vy * dt;
		a[i].zo = a[i].z - a[i].vz * dt;
	}

}

void force()
{
	int i, j;
	double Rx, Ry, Rz, R2, F;

	for (i = 0; i<N; i++)
	{
		a[i].Fx = 0;
		a[i].Fy = 0;
		a[i].Fz = 0;
	}

	for (i = 0; i<N - 1; i++)
	{
		for (j = i + 1; j<N; j++)
		{
			Rx = a[i].x - a[j].x;
			Ry = a[i].y - a[j].y;
			Rz = a[i].z - a[j].z;

			//update the boundary condition
			Rx = Rx - frameSize * (round(Rx / frameSize));
			Ry = Ry - frameSize * (round(Ry / frameSize));
			Rz = Rz - frameSize * (round(Rz / frameSize));

			R2 = Rx * Rx + Ry * Ry + Rz * Rz;
			double ecut = 4.0 / (cutoff*cutoff*cutoff*cutoff*cutoff*cutoff) * ((1 / (cutoff*cutoff*cutoff*cutoff*cutoff*cutoff)) - 1);
			PEnergy += 4.0 / (R2*R2*R2) * ((1 / (R2*R2*R2)) - 1) - ecut;
			if (R2 > cutoff*cutoff)
				continue;
			else
			{
				F = 48.0 / (R2*R2*R2*R2) * ((1 / (R2*R2*R2)) - 0.5);

				a[i].Fx += F * Rx;
				a[j].Fx -= F * Rx;

				a[i].Fy += F * Ry;
				a[j].Fy -= F * Ry;

				a[i].Fz += F * Rz;
				a[j].Fz -= F * Rz;
			}
		}
	}
}

void integrate()
{
	int i;
	double xn, yn, zn, V2 = 0;

	//repeat for x,y,z
	for (i = 0; i<N; i++)
	{
		xn = 2 * a[i].x - a[i].xo + a[i].Fx*dt*dt;
		a[i].vx = (double)(xn - a[i].xo) / (2 * dt);

		if (xn>0.5*frameSize)
		{
			a[i].xo = a[i].x - frameSize;
			a[i].x = xn - frameSize;
		}
		else if (xn<-0.5*frameSize)
		{
			a[i].xo = a[i].x + frameSize;
			a[i].x = xn + frameSize;
		}
		else
		{
			a[i].xo = a[i].x;
			a[i].x = xn;
		}

		yn = 2 * a[i].y - a[i].yo + a[i].Fy*dt*dt;
		a[i].vy = (double)(yn - a[i].yo) / (2 * dt);
		if (yn>0.5*frameSize)
		{
			a[i].yo = a[i].y - frameSize;
			a[i].y = yn - frameSize;
		}
		else if (yn<-0.5*frameSize)
		{
			a[i].yo = a[i].y + frameSize;
			a[i].y = yn + frameSize;
		}
		else
		{
			a[i].yo = a[i].y;
			a[i].y = yn;
		}

		zn = 2 * a[i].z - a[i].zo + a[i].Fz*dt*dt;
		a[i].vz = (double)(zn - a[i].zo) / (2 * dt);
		if (zn>0.5*frameSize)
		{
			a[i].zo = a[i].z - frameSize;
			a[i].z = zn - frameSize;
		}
		else if (zn<-0.5*frameSize)
		{
			a[i].zo = a[i].z + frameSize;
			a[i].z = zn + frameSize;
		}
		else
		{
			a[i].zo = a[i].z;
			a[i].z = zn;
		}

		V2 += (a[i].vx*a[i].vx + a[i].vy*a[i].vy + a[i].vz*a[i].vz);
	}
	temp = V2 / 2 / (3 * N);
	KEnergy = .5*V2;
}
