#define WALL false
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <armadillo>
#define CHUNKSIZE 1

using namespace std;

typedef vector< vector<int> > vec2Dint;
typedef vector< vector<double> > vec2Ddob;

void read_system(string fname, vec2Dint &triangles, vec2Ddob &points, vec2Ddob &velocities,
			 vec2Ddob &spline_points, vec2Ddob &spline_velocities, vector<double> &spline_width,
			 vector<double> &X,vector<double> &Y,vector<double> &Z) {
	ifstream meshfile(fname.c_str());
	string line;
	int n;

	getline(meshfile, line); // Triangles
	n = atoi(line.c_str()); 
	for (int i=0;i<n;i++) {
		getline(meshfile, line);
		istringstream buffer(line);
		vector<int> v(3);
		buffer >> v[0] >> v[1] >> v[2];
		v[0]--;v[1]--;v[2]--;
		triangles.push_back(v);
	}
	int triangle_n = n;
	
	getline(meshfile, line); // Mesh points
	n = atoi(line.c_str());
	for (int i=0;i<n;i++) {
		getline(meshfile, line);
		istringstream buffer(line);
		vector<double> v(3);
		buffer >> v[0] >> v[1] >> v[2];
		points.push_back(v);
	}

	getline(meshfile, line); // Mesh velocities
	n = atoi(line.c_str());
	for (int i=0;i<n;i++) {
		getline(meshfile, line);
		istringstream buffer(line);
		vector<double> v(3);
		buffer >> v[0] >> v[1] >> v[2];
		velocities.push_back(v);
	}

	getline(meshfile, line); // Spline points
	n = atoi(line.c_str());
	for (int i=0;i<n;i++) {
		getline(meshfile, line);
		istringstream buffer(line);
		vector<double> v(3);
		buffer >> v[0] >> v[1] >> v[2];
		spline_points.push_back(v);
	}
	int spline_n = n;

	getline(meshfile, line); // Spline velocities
	n = atoi(line.c_str());
	for (int i=0;i<n;i++) {
		getline(meshfile, line);
		istringstream buffer(line);
		vector<double> v(3);
		buffer >> v[0] >> v[1] >> v[2];
		spline_velocities.push_back(v);
	}

	getline(meshfile, line);
	n = atoi(line.c_str()); 
	for (int i=0;i<n;i++) {
		getline(meshfile, line);
		istringstream buffer(line);
		double w;
		buffer >> w;
		spline_width.push_back(w);
	}

	getline(meshfile, line);
	n = atoi(line.c_str()); 
	for (int i=0;i<n;i++) {
		getline(meshfile, line);
		istringstream buffer(line);
		vector<double> v(3);
		buffer >> v[0] >> v[1] >> v[2];
		X.push_back(v[0]);
		Y.push_back(v[1]);
		Z.push_back(v[2]);
	}
	int evaluate_n = n;

	cout << "System loaded (" << triangle_n << " triangles and " << spline_n << " spline points)" << endl;
	cout << "Evaluation of velocity field at " << evaluate_n << " points." << endl;
	return;
}

inline vector<double> diff(const vector<double> &a, const vector<double> &b) {
	vector<double> c(a.size(),0);
	for (int i=0;i<a.size();i++) {
		c[i] = a[i] - b[i];
	}
	return c;
}

double normcross(const vector<double> &u, const vector<double> &v) {
	double x = u[1]*v[0] - u[2]*v[1];
	double y = u[2]*v[1] - u[0]*v[2];
	double z = u[0]*v[2] - u[1]*v[0];
	return sqrt(x*x + y*y + z*z);
}

inline double sqr(double x) {
	return x*x;
}

inline double oseen(int i, int j, const vector<double> &r, double eps) {
	double I0 = 0;
	double n2 = sqr(r[0])+sqr(r[1])+sqr(r[2]);
	if (i==j) {
		I0 += n2+2*sqr(eps);
	}
	I0 += (r[i])*(r[j]);
	I0 /= pow(n2 + sqr(eps),1.5);
	return I0;
}

inline double oseen_wall_tensor(int i, int k, double h, const vector<double> &R, double eps) {
	double dx = eps*1e-5;
	double n2, plus, minus;
	// Plus:
	vector<double> R_p(R);
	R_p[k] += dx;
	n2 = sqr(R_p[0])+sqr(R_p[1])+sqr(R_p[2]);
	plus = h*R_p[i]/pow(n2,1.5) - oseen(i,2,R_p,eps);
	// Minus:
	vector<double> R_m(R);
	R_m[k] -= dx;
	n2 = sqr(R_m[0])+sqr(R_m[1])+sqr(R_m[2]);
	minus = h*R_m[i]/pow(n2,1.5) - oseen(i,2,R_m,eps);
	return (plus-minus)/(2*dx);
}

double oseen_wall(int i, int j, const vector<double> &y, const vector<double> &x, double eps) {
	vector<double> r = diff(y,x);
	vector<double> R = r;
	R[2] += 2*x[2];
	double R5 = pow(sqr(R[0])+sqr(R[1])+sqr(R[2])+sqr(eps),2.5);

	double I = 0;
	I += oseen(i,j,r,eps);
	I -= oseen(i,j,R,eps);
	if (j<=1) {
		I += 2*x[2]*oseen_wall_tensor(i,j,x[2],R,eps);
		if (i==j) {
			I -= 6*sqr(x[2]*eps)/R5;
		}
	} else {
		I -= 2*x[2]*oseen_wall_tensor(i,j,x[2],R,eps);
		if (i==j) {
			I += 6*sqr(x[2]*eps)/R5;
		}
	}
	double fac = 6*x[2]*sqr(eps)/R5;
	if (i==2) {
		I -= fac*R[j];
	}
	if (i==j) {
		I += fac*R[2];
	}
	return I;
}

inline double oseen_xy(int i, int j, const vector<double> &y, const vector<double> &x, double eps) {
	if (WALL) {
		return oseen_wall(i,j,y,x,eps);
	} else {
		return oseen(i,j,diff(y,x),eps);
	}
}

double int_triangle_oseen(int i, int j, const vector<double> &y, const vector<double> &xa,
							const vector<double> &xb, const vector<double> &xc,double eps) {
	const int N = 6;
	
	double weights[N] = {0.22338158967801, 0.22338158967801, 0.22338158967801, 
						 0.10995174365532, 0.10995174365532, 0.10995174365532};
	double pos1[N] = {0.44594849091597, 0.44594849091597, 0.10810301816807, 
					  0.09157621350977, 0.09157621350977, 0.81684757298046};
    double pos2[N] = {0.44594849091597, 0.10810301816807, 0.44594849091597,
    			      0.09157621350977, 0.81684757298046, 0.09157621350977};

	double area = normcross(diff(xb,xa),diff(xc,xa));
	double I = 0;
	vector<double> x(3);
	double n2,I0;
	for (int k=0;k<N;k++) {
		x[0] = xa[0]*(1-pos1[k]-pos2[k]) + xb[0]*pos1[k] + xc[0]*pos2[k];
		x[1] = xa[1]*(1-pos1[k]-pos2[k]) + xb[1]*pos1[k] + xc[1]*pos2[k];
		x[2] = xa[2]*(1-pos1[k]-pos2[k]) + xb[2]*pos1[k] + xc[2]*pos2[k];
		I0 = oseen_xy(i,j,y,x,eps);
		I += I0*weights[k];
	}
	return I;
}

double spline_oseen(int i, int j, const vector<double> &y, const vector<double> &x, double eps) {
	return oseen_xy(i,j,y,x,eps);
}

vec2Ddob assemble_matrix(const vec2Dint &t, const vec2Ddob &p,const vec2Ddob &sx,
				const vector<double> &sw, double eps) {
	int tn = t.size() + sx.size();

	vec2Ddob M(3*tn);
	for (int i=0;i<3*tn;i++) {
		M[i] = vector<double>(3*tn);
	}
	double c = 75;
	double print_percent = 0.5;
	cout << "Assembling: ";
	#pragma omp parallel for
	for (int i=0;i<tn;i++) {
		bool on_mesh = (i<t.size());
		int si = i-t.size();

		if (i/(float)tn*100 > c) {
			#pragma omp critical
			{
				cout << "." << flush;
				c += print_percent;
			}
		}
		// Find midpoint or spline point
		vector<double> tx(3);
		if (on_mesh) {
			vector<int> t1 = t[i];
			for (int j=0;j<3;j++) {
				tx[j] = (p[t1[0]][j]+p[t1[1]][j]+p[t1[2]][j])/3.0;
			}
		} else { // spline
			for (int j=0;j<3;j++) {
				tx[j] = sx[si][j];
			}
		}
		// Assemble i-th row of matrix
		for (int j=0;j<tn;j++) {
			int sj = j-t.size();
			bool on_mesh2 = sj<0;

			if (on_mesh2){
				vector<int> t2 = t[j];
				for (int ii=0;ii<3;ii++) {
					for (int jj=0;jj<3;jj++) {
						if (WALL) {
							M[tn*ii+i][tn*jj+j] = int_triangle_oseen(ii,jj,tx,p[t2[0]],p[t2[1]],p[t2[2]],eps);
						} else {
							if (jj<=ii) {
								M[tn*ii+i][tn*jj+j] = int_triangle_oseen(ii,jj,tx,p[t2[0]],p[t2[1]],p[t2[2]],eps);
								if (jj!=ii) {
									M[tn*jj+i][tn*ii+j] = M[tn*ii+i][tn*jj+j];
								}
							}
						}
					}
				}
			} else {
				for (int ii=0;ii<3;ii++) {
					for (int jj=0;jj<3;jj++) {
						if (WALL) {
							M[tn*ii+i][tn*jj+j] = spline_oseen(ii,jj,tx,sx[sj],sw[sj]);
						} else {
							if (jj<=ii) {
								M[tn*ii+i][tn*jj+j] = spline_oseen(ii,jj,tx,sx[sj],sw[sj]);
								if (jj!=ii) {
									M[tn*jj+i][tn*ii+j] = M[tn*ii+i][tn*jj+j];
								}
							}
						}
					}
				}
			}
		}
	}
	cout << endl;
	return M;
}

vector<double> solve_for_forces(const vec2Ddob &M, const vector<double> &u) {
	vector<double> M_vec(M.size()*M.size());
	for (int i=0;i<M.size();i++) {
		for (int j=0;j<M.size();j++) {
			M_vec[i+M.size()*j] = M[i][j]; // M is almost symmetric, but not quite, so order is important
		}
	}
	arma::mat A(M_vec);
	A.reshape(M.size(),M.size());
	arma::vec b(u);
	cout << "Solving for forces." << endl;
	arma::vec x = arma::solve(A,b);
	vector<double> F(x.size(),0);
	for (int i=0;i<x.size();i++) {
		F[i] = x(i);
	}
	return F;
}

void find_velocity_field(const vector<double> &X, const vector<double> &Y, const vector<double> &Z,
						 vector<double> &UX, vector<double> &UY, vector<double> &UZ,
					     const vector<double> &F, vec2Dint &t, vec2Ddob &p, const vec2Ddob &sx,
					     const vector<double> &sw,double eps) {
	// Assumes UX,UY is zeros
	
	int tn = t.size();
	int sxn = sx.size();
	int n = tn+sxn;

	double c = 75;
	double print_percent = 0.5;
	cout << "Evaluating: ";
	#pragma omp parallel for
	for (int i=0;i<X.size();i++) {
		if (i/(float)X.size()*100 > c) {
			#pragma omp critical
			{
				cout << "." << flush;
				c += print_percent;
			}
		}
		vector<double> x(3);
		x[0] = X[i];
		x[1] = Y[i];
		x[2] = Z[i];
		double s01,s02,s12;
		for (int j=0;j<tn;j++) {
			vector<int> t2 = t[j];
			if (WALL) {
				UX[i] += int_triangle_oseen(0,0,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j] 
					   + int_triangle_oseen(0,1,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j+n]
					   + int_triangle_oseen(0,2,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j+2*n];
				UY[i] += int_triangle_oseen(1,0,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j] 
					   + int_triangle_oseen(1,1,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j+n]
					   + int_triangle_oseen(1,2,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j+2*n];
				UZ[i] += int_triangle_oseen(2,0,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j] 
					   + int_triangle_oseen(2,1,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j+n]
					   + int_triangle_oseen(2,2,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j+2*n];

			} else {
				s01 = int_triangle_oseen(0,1,x,p[t2[0]],p[t2[1]],p[t2[2]],eps);
				s12 = int_triangle_oseen(1,2,x,p[t2[0]],p[t2[1]],p[t2[2]],eps);
				s02 = int_triangle_oseen(0,2,x,p[t2[0]],p[t2[1]],p[t2[2]],eps);
				UX[i] += int_triangle_oseen(0,0,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j] + s01*F[j+n] + s02*F[j+2*n];
				UY[i] += s01*F[j] + int_triangle_oseen(1,1,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j+n] + s12*F[j+2*n];
				UZ[i] += s02*F[j] + s12*F[j+n] + int_triangle_oseen(2,2,x,p[t2[0]],p[t2[1]],p[t2[2]],eps)*F[j+2*n];
			}
		}
		for (int j=0;j<sxn;j++) {
			int Fj = j+tn;
			if (WALL) {
				UX[i] += spline_oseen(0,0,x,sx[j],sw[j])*F[Fj]
					   + spline_oseen(0,1,x,sx[j],sw[j])*F[Fj+n]
					   + spline_oseen(0,2,x,sx[j],sw[j])*F[Fj+2*n];
				UY[i] += spline_oseen(1,0,x,sx[j],sw[j])*F[Fj] 
					   + spline_oseen(1,1,x,sx[j],sw[j])*F[Fj+n]
					   + spline_oseen(1,2,x,sx[j],sw[j])*F[Fj+2*n];
				UZ[i] += spline_oseen(2,0,x,sx[j],sw[j])*F[Fj] 
					   + spline_oseen(2,1,x,sx[j],sw[j])*F[Fj+n]
					   + spline_oseen(2,2,x,sx[j],sw[j])*F[Fj+2*n];
			} else {
				s01 = spline_oseen(0,1,x,sx[j],sw[j]);
				s12 = spline_oseen(1,2,x,sx[j],sw[j]);
				s02 = spline_oseen(0,2,x,sx[j],sw[j]);
				
				UX[i] += spline_oseen(0,0,x,sx[j],sw[j])*F[Fj] + s01*F[Fj+n] + s02*F[Fj+2*n];
				UY[i] += s01*F[Fj] + spline_oseen(1,1,x,sx[j],sw[j])*F[Fj+n] + s12*F[Fj+2*n];
				UZ[i] += s02*F[Fj] + s12*F[Fj+n] + spline_oseen(2,2,x,sx[j],sw[j])*F[Fj+2*n];	
			}
		}
	}
	cout << endl;
}

int main() {
	// Mesh
	vec2Dint t;
	vec2Ddob p;
	vec2Ddob v;
	vec2Ddob sx;
	vec2Ddob sv;
	vector<double> sw,X,Y,Z;
	read_system("system.txt",t,p,v,sx,sv,sw,X,Y,Z);

	// Assemble
	double mesh_eps = 0.1;
	vec2Ddob M = assemble_matrix(t,p,sx,sw,mesh_eps);

	// Boundary conditions
	vector<double> u(M.size(),0);
	int nM = t.size() + sx.size();
	for (int i=0;i<nM;i++) {
		int si = i-t.size();
		if (i<t.size()) {
			u[i] = v[i][0];
			u[i+nM] = v[i][1];
			u[i+2*nM] = v[i][2];
		} else {
			u[i] = sv[si][0];
			u[i+nM] = sv[si][1];
			u[i+2*nM] = sv[si][2];
		}
	}
	// Solve
	vector<double> F = solve_for_forces(M,u);

	// Find velocity field
	vector<double> UX(X.size(),0);
	vector<double> UY(X.size(),0);
	vector<double> UZ(X.size(),0);
	find_velocity_field(X,Y,Z,UX,UY,UZ,F,t,p,sx,sw,mesh_eps);

	ofstream velfile("velocity.txt");
	for (int i=0;i<X.size();i++) {
		velfile << X[i] << " " << Y[i] << " " << Z[i] << " " << UX[i] << " " << UY[i] << " " << UZ[i] << "\n";
	}
	velfile.close();

	// Find velocity at BC'sfield
	if (false) {
		vector<double> Xbc;
		vector<double> Ybc;
		vector<double> Zbc;

		int tn = t.size() + sx.size();
		for (int i=0;i<tn;i++) {
			bool on_mesh = (i<t.size());
			int si = i-t.size();
			if (on_mesh) {
				vector<int> t1 = t[i];
				Xbc.push_back((p[t1[0]][0]+p[t1[1]][0]+p[t1[2]][0])/3.0);
				Ybc.push_back((p[t1[0]][1]+p[t1[1]][1]+p[t1[2]][1])/3.0);
				Zbc.push_back((p[t1[0]][2]+p[t1[1]][2]+p[t1[2]][2])/3.0);
			} else { // spline
				Xbc.push_back(sx[si][0]);
				Ybc.push_back(sx[si][1]);
				Zbc.push_back(sx[si][2]);
			}
		}
		vector<double> UXbc(Xbc.size(),0);
		vector<double> UYbc(Xbc.size(),0);
		vector<double> UZbc(Xbc.size(),0);
		find_velocity_field(Xbc,Ybc,Zbc,UXbc,UYbc,UZbc,F,t,p,sx,sw,mesh_eps);
		// Save to txt
		ofstream velbcfile("velocity_bc.txt");
		for (int i=0;i<Xbc.size();i++) {
			velbcfile << Xbc[i] << " " << Ybc[i] << " " << Zbc[i] << " " << UXbc[i] << " " << UYbc[i] << " " << UZbc[i] << "\n";
		}
		velbcfile.close();
	}

	return 0;
}