#include "timeStepper.h"

timeStepper::timeStepper(elasticRod &m_rod)
{
	rod = &m_rod;
	kl = 10; // lower diagonals
	ku = 10; // upper diagonals
	freeDOF = rod->uncons;
	ldb = freeDOF;
	NUMROWS = 2 * kl + ku + 1;
	totalForce = new double[freeDOF];
	nrhs = 1;
    ipiv = new int[freeDOF];
    info = 0;

    //jacobian
    jacobianLen = (2 * kl + ku + 1) * freeDOF;
    jacobian = new double [jacobianLen];

    // scalesolver = rod->massArray[0]*(rod->rodLength/rod->ne) /(pow(rod->dt,2));
}

timeStepper::~timeStepper()
{
	;
}

double* timeStepper::getForce()
{
	return totalForce;
}

double* timeStepper::getJacobian()
{
	return jacobian;
}

// void timeStepper::addForce(int ind, double p)
// {
//     if (ind < rod->ndof)
//     {
//         if (rod->getIfConstrained(ind) == 0) // free dof
//         {
//             mappedInd = rod->fullToUnconsMap[ind];
//             // totalForce[mappedInd] = totalForce[mappedInd] + p; // subtracting elastic force
//             Force[mappedInd] = Force[mappedInd] + p;
//         }
//         force(ind) = force(ind) + p;

//     }
//     // else
//     // {
//     //     ind = ind - rod->ndof + rod->uncons;
//     //     // totalForce[ind] = totalForce[ind] + p; // subtracting elastic force
//     //     Force[ind] = Force[ind] + p;
//     // }
//     // force(ind) = force(ind) + p;
// }

void timeStepper::addForce(int ind, double p)
{
    if (rod->getIfConstrained(ind) == 0) // free dof
    {
        mappedInd = rod->fullToUnconsMap[ind];
        totalForce[mappedInd] = totalForce[mappedInd] + p; // subtracting elastic force
        Force[mappedInd] = Force[mappedInd] + p;
    }
    force(ind) = force(ind) + p;
}

void timeStepper::addJacobian(int ind1, int ind2, double p)
{
    mappedInd1 = rod->fullToUnconsMap[ind1];
    mappedInd2 = rod->fullToUnconsMap[ind2];
    if (rod->getIfConstrained(ind1) == 0 && rod->getIfConstrained(ind2) == 0) // both are free
    {
        row = kl + ku + mappedInd2 - mappedInd1;
        col = mappedInd1;
        offset = row + col * NUMROWS;
        jacobian[offset] = jacobian[offset] + p;
        Jacobian(mappedInd2, mappedInd1) = Jacobian(mappedInd2, mappedInd1) + p;
    }
}


// void timeStepper::addJacobian(int ind1, int ind2, double p)
// {

//     mappedInd1 = rod->fullToUnconsMap[ind1];
//     mappedInd2 = rod->fullToUnconsMap[ind2];
//     if (rod->getIfConstrained(ind1) == 0 && rod->getIfConstrained(ind2) == 0) // both are free
//     {
//         // row = kl + ku + mappedInd2 - mappedInd1;
//         // col = mappedInd1;
//         // offset = row + col * NUMROWS;
//         // jacobian[offset] = jacobian[offset] + p;
//         Jacobian(mappedInd2, mappedInd1) = Jacobian(mappedInd2, mappedInd1) + p;
//     }

//     // bool flag = true;
//     // if (ind1 < rod->ndof)
//     // {
//     //     if (rod->getIfConstrained(ind1) == 0)
//     //     {
//     //         mappedInd1 = rod->fullToUnconsMap[ind1];
//     //     }
//     //     else
//     //     {
//     //         flag = false;
//     //     }
//     // }
//     // else
//     // {
//     //     mappedInd1 = ind1 - rod->ndof + rod->uncons;
//     // }

//     // if (ind2 < rod->ndof)
//     // {
//     //     if (rod->getIfConstrained(ind2) == 0)
//     //     {
//     //         mappedInd2 = rod->fullToUnconsMap[ind2];
//     //     }
//     //     else
//     //     {
//     //         flag = false;
//     //     }
//     // }
//     // else
//     // {
//     //     mappedInd2 = ind2 - rod->ndof + rod->uncons;
//     // }

//     // if (flag)
//     // {
//     //     Jacobian(mappedInd1, mappedInd2) = Jacobian(mappedInd1, mappedInd2) + p;
//     // }
// }

void timeStepper::addForce_c(int ind, double p)
{
    Force(ind) = Force(ind) + p;
}

void timeStepper::addJacobian_c(int ind1, int ind2, double p)
{
    Jacobian(ind1, ind2) = Jacobian(ind1, ind2) + p;
}



void timeStepper::setZero()
{
	// for (int i=0; i < freeDOF; i++)
	// 	totalForce[i] = 0;
	Force = VectorXd::Zero(freeDOF);
	Jacobian = MatrixXd::Zero(freeDOF,freeDOF);
	force = VectorXd::Zero(rod->ndof);

	for (int i=0; i < freeDOF; i++)
	{
		totalForce[i] = 0;
	}
  for (int i=0; i < jacobianLen; i++)
	{
		jacobian[i] = 0;
	}
	force = VectorXd::Zero(rod->ndof);
}

void timeStepper::update(int num)
{
    freeDOF = num;
    delete []totalForce;
    totalForce = NULL;
    totalForce = new double[freeDOF];
    DX = VectorXd::Zero(freeDOF);
    // cout<<"freeDOF_new: "<<freeDOF<<endl;
}


void timeStepper::updateforeigen(int num)
{
    freeDOF = rod->uncons;
    ldb = freeDOF;
    delete [] totalForce;
    delete [] jacobian;
    delete [] ipiv;
    totalForce = new double[freeDOF];
    jacobianLen = (2 * kl + ku + 1) * freeDOF;
    jacobian = new double [jacobianLen];
    nrhs = 1;
    ipiv = new int[freeDOF];
    info = 0;
}

//void timeStepper::pardisoSolver()
//{
//
//    int    n = freeDOF;
//
//    int ia[n+1];
//    ia[0] = 0;
//
//    int temp = 0;
//
//    for (int i =0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            if (Jacobian(i,j) != 0)
//            {
//                temp = temp + 1;
//            }
//        }
//        ia[i+1] = temp;
//    }
//
//    int ja[ia[n]];
//    double a[ia[n]];
//
//    temp = 0;
//
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            if (Jacobian(i,j) != 0)
//            {
//                ja[temp] = j;
//                a[temp] = Jacobian(i,j);
//
//                temp = temp + 1;
//            }
//        }
//    }
//
//    int      nnz = ia[n];
//    int      mtype = 11;        /* Real symmetric matrix */
//
//    /* RHS and solution vectors. */
//    double   b[n], x[n];
//    int      nrhs = 1;          /* Number of right hand sides. */
//
//    /* Internal solver memory pointer pt,                  */
//    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
//    /* or void *pt[64] should be OK on both architectures  */
//    void    *pt[64];
//
//    /* Pardiso control parameters. */
//    int      iparm[64];
//    double   dparm[64];
//    int      maxfct, mnum, phase, error, msglvl, solver;
//
//    /* Number of processors. */
//    int      num_procs;
//
//    /* Auxiliary variables. */
//    char    *var;
//    int      i;
//
//    double   ddum;              /* Double dummy */
//    int      idum;              /* Integer dummy. */
//
//
///* -------------------------------------------------------------------- */
///* ..  Setup Pardiso control parameters.                                */
///* -------------------------------------------------------------------- */
//
//    error = 0;
//    solver = 0; /* use sparse direct solver */
//    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);
//
//    /* Numbers of processors, value of OMP_NUM_THREADS */
//    var = getenv("OMP_NUM_THREADS");
//    if(var != NULL)
//        sscanf( var, "%d", &num_procs );
//    else {
//        printf("Set environment OMP_NUM_THREADS to 1");
//        exit(1);
//    }
//    iparm[2]  = num_procs;
//
//    maxfct = 1;     /* Maximum number of numerical factorizations.  */
//    mnum   = 1;         /* Which factorization to use. */
//
//    msglvl = 0;         /* Print statistical information  */
//    error  = 0;         /* Initialize error flag */
//
///* -------------------------------------------------------------------- */
///* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
///*     notation.                                                        */
///* -------------------------------------------------------------------- */
//    for (i = 0; i < n+1; i++) {
//        ia[i] += 1;
//    }
//    for (i = 0; i < nnz; i++) {
//        ja[i] += 1;
//    }
//
//    /* Set right hand side to one. */
//    for (i = 0; i < n; i++) {
//        b[i] = Force[i];
//    }
//
///* -------------------------------------------------------------------- */
///*  .. pardiso_chk_matrix(...)                                          */
///*     Checks the consistency of the given matrix.                      */
///*     Use this functionality only for debugging purposes               */
///* -------------------------------------------------------------------- */
//
//    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
//    if (error != 0) {
//        printf("\nERROR in consistency of matrix: %d", error);
//        exit(1);
//    }
//
///* -------------------------------------------------------------------- */
///* ..  pardiso_chkvec(...)                                              */
///*     Checks the given vectors for infinite and NaN values             */
///*     Input parameters (see PARDISO user manual for a description):    */
///*     Use this functionality only for debugging purposes               */
///* -------------------------------------------------------------------- */
//
//    pardiso_chkvec (&n, &nrhs, b, &error);
//    if (error != 0) {
//        printf("\nERROR  in right hand side: %d", error);
//        exit(1);
//    }
//
//
///* -------------------------------------------------------------------- */
///* .. pardiso_printstats(...)                                           */
///*    prints information on the matrix to STDOUT.                       */
///*    Use this functionality only for debugging purposes                */
///* -------------------------------------------------------------------- */
//
//    // pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
//    // if (error != 0) {
//    //     printf("\nERROR right hand side: %d", error);
//    //     exit(1);
//    // }
//
///* -------------------------------------------------------------------- */
///* ..  Reordering and Symbolic Factorization.  This step also allocates */
///*     all memory that is necessary for the factorization.              */
///* -------------------------------------------------------------------- */
//    phase = 11;
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//         &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error, dparm);
//
//    if (error != 0) {
//        printf("\nERROR during symbolic factorization: %d", error);
//        exit(1);
//    }
//    // printf("\nReordering completed ... ");
//    // printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
//    // printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
//
///* -------------------------------------------------------------------- */
///* ..  Numerical factorization.                                         */
///* -------------------------------------------------------------------- */
//    phase = 22;
//    iparm[32] = 1; /* compute determinant */
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//
//    if (error != 0) {
//        printf("\nERROR during numerical factorization: %d", error);
//        exit(2);
//    }
//    // printf("\nFactorization completed ...\n ");
//
///* -------------------------------------------------------------------- */
///* ..  Back substitution and iterative refinement.                      */
///* -------------------------------------------------------------------- */
//    phase = 33;
//
//    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, b, x, &error,  dparm);
//
//    if (error != 0) {
//        printf("\nERROR during solution: %d", error);
//        exit(3);
//    }
//    for (i = 0; i < n; i++)
//    {
//        totalForce[i] = x[i];
//        DX(i) = x[i];
//    }
//
//
///* -------------------------------------------------------------------- */
///* ..  Convert matrix back to 0-based C-notation.                       */
///* -------------------------------------------------------------------- */
//    for (i = 0; i < n+1; i++) {
//        ia[i] -= 1;
//    }
//    for (i = 0; i < nnz; i++) {
//        ja[i] -= 1;
//    }
//
///* -------------------------------------------------------------------- */
///* ..  Termination and release of memory.                               */
///* -------------------------------------------------------------------- */
//    phase = -1;                 /* Release internal memory. */
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, &ddum, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//
//
//}

void timeStepper::pardisoSolver()
{
    int n = freeDOF;
    int ia[n+1];
    ia[0] = 1;

    int temp = 0;
    for (int i =0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (Jacobian(i,j) != 0)
            {
                temp = temp + 1;
            }
        }
        ia[i+1] = temp+1;
    }

    int ja[ia[n]];
    double a[ia[n]];
    temp = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (Jacobian(i,j) != 0)
            {
                ja[temp] = j + 1;
                a[temp] = Jacobian(i,j);
                temp = temp + 1;
            }
        }
    }
    MKL_INT mtype = 11;       /* Real unsymmetric matrix */
    // Descriptor of main sparse matrix properties
    double b[n], x[n], bs[n], res, res0;
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i, j;
    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 0;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */


    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information  */
    error = 0;            /* Initialize error flag */\
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: " IFORMAT, error);
        exit (1);
    }
    // printf ("\nReordering completed ... ");
    // printf ("\nNumber of nonzeros in factors = " IFORMAT, iparm[17]);
    // printf ("\nNumber of factorization MFLOPS = " IFORMAT, iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during numerical factorization: " IFORMAT, error);
        exit (2);
    }
    // printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
    phase = 33;

// descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
// descrA.mode = SPARSE_FILL_MODE_UPPER;
// descrA.diag = SPARSE_DIAG_NON_UNIT;
// mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ONE, n, n, ia, ia+1, ja, a );

    /* Set right hand side to one. */
    for ( i = 0; i < n; i++ )
    {
        b[i] = Force[i];
    }
//  Loop over 3 solving steps: Ax=b, AHx=b and ATx=b
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during solution: " IFORMAT, error);
        exit (3);
    }

    /* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);

    for (int i = 0; i < n; i++)
    {
        totalForce[i] = x[i];
        DX(i) = x[i];
    }

}



void timeStepper::integrator()
{
    // cout<<"jacobian: "<<Jacobian.size()<<endl;
    // MatrixXd temp;
    // temp = Jacobian - Jacobian.transpose();
    // cout<<temp.norm()<<endl;
	// dgbsv_(&freeDOF, &kl, &ku, &nrhs, jacobian, &NUMROWS, ipiv, totalForce, &ldb, &info);
    pardisoSolver();
}

void timeStepper::integrator_eigen()
{
    // dgbsv_(&freeDOF, &kl, &ku, &nrhs, jacobian, &NUMROWS, ipiv, totalForce, &ldb, &info);
    dgbsv_(&freeDOF, &kl, &ku, &nrhs, jacobian, &NUMROWS, ipiv, totalForce, &ldb, &info);
}


void timeStepper::monitor()
{
    double f = 0;
    for (int i = 0; i< freeDOF; i++)
    {
        f = f + totalForce[i] * totalForce[i];

    }
    cout<<"totalForce: "<<f<<endl;;
}



void timeStepper::setFe()
{

    for (int i = 0; i< rod->ndof; i++)
    {
        if (rod->isConstrained[i]!= 0)
        {
            F_ei(i) = force(i);
        }
    }
}

void timeStepper::setFi()
{
    F_ei = VectorXd::Zero(rod->ndof);
    for (int i = 0; i< rod->ndof; i++)
    {
        if (rod->isConstrained[i] ==0)
        {
            F_ei(i) = -force(i);
        }

    }
}
