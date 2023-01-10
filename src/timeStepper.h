#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include "elasticRod.h"

extern "C" void dgbsv_( int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info );

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);


class timeStepper
{
public:
	timeStepper(elasticRod &m_rod);
	~timeStepper();
	double* getForce();
	double* getJacobian();
	void setZero();
	void addForce(int ind, double p);
    void addForce1(int ind, double p);

	void addJacobian(int ind1, int ind2, double p);
	void addJacobian1(int ind1, int ind2, double p);

	void integrator();

	void pardisoSolver();

	void monitor();

	VectorXd force;

	void update(int num);
	int freeDOF;

	MatrixXd Jacobian;
	VectorXd Force;

	double computeD();
	VectorXd DX;

	void getdirection();

	VectorXd F_ei;

    void setFe();

    void setFi();

    void addForce_c(int ind, double p);

    void addJacobian_c(int ind1, int ind2, double p);

    void updateforeigen(int num);
    void integrator_eigen();

    // double scalesolver;



private:
	elasticRod *rod;
	int kl, ku;
	
	double *totalForce;
	double *jacobian;

	
	// utility variables
	int mappedInd, mappedInd1, mappedInd2;
	int row, col, offset;
	int NUMROWS;
	int jacobianLen;
	int nrhs;
    int *ipiv;
    int info;
    int ldb;
};

#endif
