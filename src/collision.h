#ifndef COLLISION_H
#define COLLISION_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"


struct collisionP{
    int node1;
    int node2;
    Vector3d P; //P is the distance vector between two segments
    double distance;

    double sc; //ratio in segment i -> i+1
    double tc; //ratio in segment j -> j+1
};

struct candidateP{
    int i;
    int j;

    double sc;
    double tc;

    Vector3d n;
    Vector3d t1;
    Vector3d t2;

    double distance;
};


class collision
{
public:
    collision(elasticRod &m_rod, timeStepper &m_stepper, double friction, double m_r);
    ~collision();
    elasticRod *rod;
    timeStepper *stepper;

    void setZero();
    bool candidateTest();

    void computeFc();
    void computeJc();

    void collisionTest();

    void initialize();
    void updateForce(VectorXd &dx, double a);

    void initialguess();
    double updaterho(double rho);
    double evaluation(double rho_k);
    void projector(VectorXd &z);
    double fixbound(double num);
    void prepareForSolve();
    void prepareForIteration();
    void computeFac();
    void computeFacandJac();
    void possibleCandidate();
    void updatedn();

    int p_pairs;
    VectorXd lamda;
    VectorXd lamda_old;
    VectorXd Fc;
    MatrixXd jacobian;



private:
    // parameter for collision check
    double SMALL_NUMBER;
    Vector3d p0;
    Vector3d p1;
    Vector3d q0;
    Vector3d q1;
    Vector3d u;
    Vector3d v;
    Vector3d w;
    double tn;
    double td;
    double sn;
    double sd;
    double tc;
    double sc;
    double a;
    double b;
    double c;
    double d;
    double e;
    double D;
    collisionP detectionP;
    int interval;

    double friction;


    //mininum distance between rod segment
    double d0;
    //ref length between two nodes;
    double delta_l;

    //vector for collision candidate
    candidateP c_p;
    vector <candidateP> possibleC;

    VectorXd Fn;

    MatrixXd J11;
    MatrixXd J12;
    MatrixXd J21;
    MatrixXd J22;

    //Vector for calculate distance
    VectorXd dn;
    double r; //penalty coefficient

    MatrixXd q;

    MatrixXd H;

    VectorXd F_k;
    VectorXd F_ki;

    VectorXd z_k;
    VectorXd z_ki;

    //direction
    Vector3d n;
    Vector3d t1;
    Vector3d t2;

    Vector3d l_n;
    Vector3d l_t;

    MatrixXd Id3;

    Vector3d sigma_n;
    Vector3d sigma_t;

    Vector3d delta_t;

    double f_n;

    MatrixXd W;
    VectorXd b1;
    // MatrixXd Jacobian;

    VectorXd u_k;
    VectorXd b_k;

    MatrixXd P_n;
    MatrixXd P_t;

    MatrixXd rm;
    MatrixXd um;
    double rn;
    double un;

    Vector2d rt;
    Vector2d ut;

    Vector3d lam;

    MatrixXd jac_u;
    MatrixXd jac_r;

    Vector3d fac;

    MatrixXd jacobian_u;
    MatrixXd jacobian_r;

    double delta;
    Vector2d y;

    Matrix2d Id2;
    MatrixXd Pn;
    MatrixXd Pt;


    void projectFk(VectorXd &F);

    double computeAC(VectorXd &F, VectorXd &Z);
    void computeDistance();



};
#endif
