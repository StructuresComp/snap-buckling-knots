#include "collision.h"


collision::collision(elasticRod &m_rod, timeStepper &m_stepper, double m_friction, double m_r)
{
	rod = &m_rod;
	stepper = &m_stepper;
	SMALL_NUMBER = 1e-4;
	d0 = 2*rod->rodRadius;
	delta_l = rod->rodLength/(rod->nv - 1);
  interval = 1.5*d0/delta_l + 1;
  Fc = VectorXd::Zero(rod->uncons);
  Id3 = MatrixXd::Identity(3, 3);
  Id2 = MatrixXd::Identity(2, 2);
	r = m_r;
  friction = m_friction;
}

collision::~collision()
{
	;
}

double collision::fixbound(double num)
{
    if (num < 0)
        num = 0;
    else if (num > 1)
        num = 1;

    return num;
}


void collision::computeDistance()
{
    u = p1 - p0; //d1
    v = q1 - q0; //d2
    w = q0 - p0; //d12

    a = u.dot(u); //D1
    b = u.dot(v); //R
    c = v.dot(v); //D2
    d = u.dot(w); //S1
    e = v.dot(w); //S2
    D = a*c - b*b; //den

    double uf;
    if (D == 0)
    {
        sc = 0;
        tc = -e/c;
        uf = fixbound(tc);

        if (uf != tc)
        {
            sc = (uf * b + d)/a;
            sc = fixbound(sc);
            tc = uf;
        }
    }
    else
    {
        sc = (d*c - e *b)/D;
        sc = fixbound(sc);
        tc = (sc * b - e)/c;
        uf = fixbound(tc);

        if (uf != tc)
        {
            sc = (uf * b + d)/a;
            sc = fixbound(sc);
            tc = uf;
        }
    }

    detectionP.P = u*sc-tc*v-w;
    detectionP.distance = detectionP.P.norm();

    detectionP.P = detectionP.P/detectionP.distance;
    detectionP.sc = sc;
    detectionP.tc = tc;
}

void collision::setZero()
{
    possibleC.clear();
    p_pairs = 0;
}

bool collision::candidateTest()
{
    // constraint candidates selection

    double min = 1;
    for (int i = 0; i < rod->nv-2; i++)
    {
        p0 = rod->getVertex(i);
        p1 = rod->getVertex(i+1);
        for (int j = i+2; j < rod->nv-1; j++)
        {
            q0 = rod->getVertex(j);
            q1 = rod->getVertex(j+1);
            computeDistance();
            if (detectionP.distance <1.0*d0  &&  abs(i - j) > interval)
            {
                c_p.i = i;
                c_p.j = j;
                possibleC.push_back(c_p);
            }
        }
    }

    if (possibleC.size()!=0)
    {
    	p_pairs = possibleC.size();
      return true;
    }
    else
        return false;
}

void collision::possibleCandidate()
{
	int c = 0;
	vector<candidateP>::iterator it;
	dn = VectorXd::Zero(p_pairs);
	for (it = possibleC.begin(); it != possibleC.end(); it++)
	{
		int i = it->i;
    int j = it->j;

    p0 = rod->getVertex(i);
    p1 = rod->getVertex(i+1);
    q0 = rod->getVertex(j);
    q1 = rod->getVertex(j+1);

    computeDistance();

    it->sc = detectionP.sc;
    it->tc = detectionP.tc;
    n = detectionP.P;
    Vector3d temp(1, 0, 0);
    t1 = n.cross(temp);
    if (t1.norm()==0)
    {
			temp = Vector3d(1,1,0);
      t1 = n.cross(temp);
		}
		t1 = t1/t1.norm();
    t2 = n.cross(t1);
    t2 = t2/t2.norm();
    it->n = n;
    it->t1 = t1;
    it->t2 = t2;
		dn(c) = (d0 - detectionP.distance)/rod->dt;
		c = c + 1;
  }
}

void collision::updatedn()
{
    int c = 0;
    vector<candidateP>::iterator it;
    dn = VectorXd::Zero(p_pairs);
    for (it = possibleC.begin(); it != possibleC.end(); it++)
    {
        int i = it->i;
        int j = it->j;

        p0 = rod->getVertex(i);
        p1 = rod->getVertex(i+1);
        q0 = rod->getVertex(j);
        q1 = rod->getVertex(j+1);

        u = p1 - p0; //d1
        v = q1 - q0; //d2
        w = q0 - p0; //d12

        detectionP.P = u*it->sc-it->tc*v-w;
        detectionP.distance = detectionP.P.norm();
        dn(c) = (d0 - detectionP.distance)/rod->dt;

        c = c + 1;
    }

}

void collision::prepareForSolve()
{
    H = MatrixXd::Zero(rod->ndof, 3*p_pairs);
    Pn = MatrixXd::Zero(3*p_pairs, p_pairs);
    Pt = MatrixXd::Zero(3*p_pairs, 2*p_pairs);

    int c = 0;
    vector<candidateP>::iterator it;
    for (it = possibleC.begin(); it != possibleC.end(); it++)
    {
        int i, j;
        double sc, tc;
        i = it->i;
        j = it->j;
        sc = it->sc;
        tc = it->tc;
        n = it->n;
        t1 = it->t1;
        t2 = it->t2;

        H.block<3,3>(4*i, 3*c) = Id3 * (1-sc);
        H.block<3,3>(4*(i+1), 3*c) = Id3 * sc;
        H.block<3,3>(4*j, 3*c) = -Id3 *(1-tc);
        H.block<3,3>(4*(j+1), 3*c) = -Id3 * tc;
        Pn.block<3,1>(3*c, c) = n;
        Pt.block<3,1>(3*c, 2*c) = t1;
        Pt.block<3,1>(3*c, 2*c+1) = t2;

        c = c + 1;
    }

    W = rod->dt * H.transpose() * rod->inv_M * H; //3m x 3m
    b1 = rod->dt* H.transpose() * rod->inv_M * stepper->F_ei +H.transpose() * rod->u;
    b1 = b1 - Pn * dn;
}

void collision::prepareForIteration()
{
    u_k = W * lamda + b1;
}


void collision::computeFacandJac()
{
    int c = 0;
    vector<candidateP>::iterator it;
    Fc = VectorXd::Zero(3*p_pairs);
    jacobian = MatrixXd::Zero(3*p_pairs, 3*p_pairs);
    jacobian_u = MatrixXd::Zero(3*p_pairs, 3*p_pairs);
    jacobian_r = MatrixXd::Zero(3*p_pairs, 3*p_pairs);
    // cout<<"tdz_b1"<<endl;
    for (it = possibleC.begin(); it != possibleC.end(); it++)
    {
        // cout<<"c: "<<c<<endl;
        VectorXd temp;
        n = it->n;
        t1= it->t1;
        t2 = it->t2;
        P_n = MatrixXd::Zero(1,3);
        P_t = MatrixXd::Zero(2,3);

        P_n.block<1,3>(0,0) = n;
        P_t.block<1,3>(0,0) = t1;
        P_t.block<1,3>(1,0) = t2;

        lam = lamda.segment<3>(3*c);
        u = u_k.segment<3>(3*c);

        jac_u = MatrixXd::Zero(3,3);
        jac_r = MatrixXd::Zero(3,3);

        rm =  P_n * lam;
        rn = rm(0);
        um = P_n * u ;
        un = um(0);


        if (rn - r * un <= 0) //gap
        {
            temp = -P_n *lam;
            fac(0) = temp(0);

            jac_r.block<1,3>(0,0) = - P_n;
        }
        else //contact
        {
            temp = - r * P_n * u;
            fac(0) = temp(0);
            jac_u.block<1,3>(0,0) = -r * P_n;
        }

        delta = friction *rn;
        rt = P_t * lam;
        ut = P_t * u ;

        y = rt - r * ut;

        if (y.norm()<= delta) //stick
        {
            fac.segment<2>(1) = -r * P_t * u ; //2x1
            jac_u.block<2,3>(1,0) = -r * P_t;
        }

        if (y.norm() > delta && delta > 0) //slip
        {
            fac.segment<2>(1) = delta * y/y.norm() - rt; //2x1
            jac_u.block<2,3>(1,0) = -r * delta/y.norm() * (Id2 - y * y.transpose()/pow(y.norm(), 2))
                                    * P_t;
            jac_r.block<2,3>(1,0) = friction * y/y.norm() * P_n  +  delta/y.norm() * (Id2 -
                                   y*y.transpose()/pow(y.norm(), 2))*P_t - P_t;
        }

        if (y.norm()>=0 && delta<=0)//gap
        {
            fac.segment<2>(1) = -P_t*lam;
            jac_r.block<2,3>(1,0) = -P_t;
        }

        Fc.segment<3>(3*c) = fac;
        jacobian_u.block<3,3>(3*c, 3*c) = jac_u;
        jacobian_r.block<3,3>(3*c, 3*c) = jac_r;

        c = c+1;
    }
    jacobian = jacobian_u * W +  jacobian_r;

    MatrixXd temp = jacobian_u * W;


    for (int i = 0; i<3*p_pairs; i++)
    {
        stepper->addForce_c(i, Fc(i));
    }

    for (int i = 0; i<3*p_pairs; i++)
    {
        for (int j = 0; j<3*p_pairs; j++)
        {
            stepper->addJacobian_c(i,j, jacobian(i, j));
        }
    }
}



void collision::computeFac()
{
    int c = 0;
    vector<candidateP>::iterator it;
    Fc = VectorXd::Zero(3*p_pairs);
    for (it = possibleC.begin(); it != possibleC.end(); it++)
    {
        VectorXd temp;
        n = it->n;
        t1= it->t1;
        t2 = it->t2;
        P_n = MatrixXd::Zero(1,3);
        P_t = MatrixXd::Zero(2,3);

        P_n.block<1,3>(0,0) = n;
        P_t.block<1,3>(0,0) = t1;
        P_t.block<1,3>(1,0) = t2;

        lam = lamda.segment<3>(3*c);
        u = u_k.segment<3>(3*c);

        rm =  P_n * lam;
        rn = rm(0);
        um = P_n * u ;
        un = um(0);
        //calculate normal part
        if (rn - r * un <= 0) //gap
        {
            // cout<<"gap"<<endl;
            temp = -P_n *lam;
            fac(0) = temp(0);
        }
        else //contact
        {
            temp = - r * P_n * u;
            fac(0) = temp(0);
        }

        //calculate friction part;
        delta = friction *rn;
        rt = P_t * lam;
        ut = P_t * u ;

        y = rt - r * ut;

        if (y.norm()<= delta) //stick
        {
            fac.segment<2>(1) = -r * P_t * u ; //2x1
        }

        if (y.norm() > delta && delta > 0) //slip
        {
            fac.segment<2>(1) = delta * y/y.norm() - rt; //2x1
        }

        if (y.norm()>=0 && delta<=0)//gap
        {
            fac.segment<2>(1) = -P_t*lam;
        }

        Fc.segment<3>(3*c) = fac;
        c = c+1;
    }
    MatrixXd temp = jacobian_u * W;


    for (int i = 0; i<3*p_pairs; i++)
    {
        stepper->addForce_c(i, Fc(i));
    }
}

void collision::initialize()
{
	lamda = VectorXd::Zero(3*p_pairs);

}


void collision::collisionTest()
{
    ;
}


void collision::computeFc()
{
	Fn = VectorXd::Zero(rod->ndof);
  Fn = H * lamda;

	for (int i = 0; i< rod->ndof; i++)
	{
		stepper->addForce(i, -Fn(i));
	}
}


void collision::computeJc()
{
	for (int i = 0; i<rod->ndof; i++)
	{
		for (int j = 0; j<rod->ndof; j++)
		{
			stepper->addJacobian(i, j, J11(i,j));
		}
	}

	for (int i = 0; i<rod->ndof; i++)
	{
		for (int j = 0; j<3*p_pairs; j++)
		{
			stepper->addJacobian(i, j+rod->ndof, J12(i,j));
			stepper->addJacobian(j+rod->ndof, i, J21(j,i));
		}
	}

	for (int i = 0; i<3*p_pairs; i++)
	{
		stepper->addJacobian(i+rod->ndof, i+rod->ndof, J22(i, i));
	}
}


void collision::updateForce(VectorXd &dx, double a)
{

	for (int i = 0; i< 3 * p_pairs; i++)
	{
		lamda(i) = lamda(i) - a * dx[i];
	}
}


void collision::initialguess()
{
    z_k = lamda;
    F_k = W * z_k + b1; //relative velocity

    double tol = 1e-9;
    int iter_max = 100;
    int iter = 0;

    double error;
    error = computeAC(F_k, z_k);
    // cout<<"error: "<<error<<endl;
    double rho_k;
    rho_k = 500;

    while (error > tol && iter<iter_max)
    {
        //update rho_k
        rho_k = updaterho(rho_k);
        z_ki = z_k - rho_k * F_k;
        projector(z_ki);
        F_ki = W * z_ki + b1;
        z_k = z_k - rho_k *F_ki;
        projector(z_k);
        F_k = W * z_k + b1;
        error = computeAC(F_k, z_k);
        iter = iter + 1;

        // cout<<"iter: "<<iter<<" error: "<<error<<" rho_k: "<<rho_k<<endl;
    }
    // exit(0);
    r = rho_k;
    lamda = z_k;
}


double collision::computeAC(VectorXd &F, VectorXd &Z)
{
	int c = 0;

	VectorXd Fac;
	Fac = VectorXd::Zero(3*p_pairs);
    vector<candidateP>::iterator it;
    for (it = possibleC.begin(); it != possibleC.end(); it++)
    {
    	n = it->n;
        t1= it->t1;
        t2 = it->t2;

        P_n = MatrixXd::Zero(1,3);
        P_t = MatrixXd::Zero(2,3);

        P_n.block<1,3>(0,0) = n;
        P_t.block<1,3>(0,0) = t1;
        P_t.block<1,3>(1,0) = t2;

        fac = Vector3d(0,0,0);

        lam = Z.segment<3>(3*c);

        u = F.segment<3>(3*c);

        rm =  P_n * lam;
        rn = rm(0);
        um = P_n * u;
        un = um(0);
        //calculate normal part
        if (rn - r * un <= 0) //gap
        {
            // temp = -P_n *lam;
            fac(0) = -rn;
        }
        else //contact
        {
            // temp = - r * P_n * u;
            fac(0) = -r * un;
        }

        //calculate friction part;
        delta = friction * rn;
        rt = P_t * lam;
        ut = P_t * u ;

        y = rt - r * ut;

        if (y.norm()<= delta) //stick
        {
            fac.segment<2>(1) = -r *ut; //2x1
        }

        if (y.norm() > delta && delta > 0) //slip
        {
            fac.segment<2>(1) = delta * y/y.norm() - rt; //2x1
        }

        if (y.norm()>=0 && delta<=0)
        {
            fac.segment<2>(1) = -rt;
        }
        Fac.segment<3>(3*c) = fac;

        c = c+1;
    }

    double error = Fac.norm();
    return error;

}

void collision::projectFk(VectorXd &F)
{
    int c = 0;
    vector<candidateP>::iterator it;
    for (it = possibleC.begin(); it != possibleC.end(); it++)
    {
        n = it->n;
        t1= it->t1;
        t2 = it->t2;

        u = F.segment<3>(3*c);

        P_t = MatrixXd::Zero(2,3);
        P_t.block<1,3>(0,0) = t1;
        P_t.block<1,3>(1,0) = t2;
        ut = P_t * u;



        c = c +1;
    }
}


double collision::updaterho(double rho)
{
    double L_min = 0.2;
    double L = 0.8;
    double v = 0.5;

    double rho_k = rho;

    z_ki = z_k - rho_k * F_k;
    projector(z_ki);
    double rk;

    F_ki = W * z_ki + b1;

    rk = evaluation(rho_k);
    // cout<<rk<<endl;
    // exit(0);

    // int iter = 0;
    while (rk > L)
    {
        rho_k = v * rho_k;
        F_k = W * z_k + b1;
        z_ki = z_k - rho_k * F_k;
        projector(z_ki);
        F_ki = W * z_ki + b1;
        rk = evaluation(rho_k);
    }
    if (rk < L_min)
        rho_k = 1.0/v * rho_k;

    return rho_k;
}

double collision::evaluation(double rho_k)
{
    double r;
    VectorXd num, den;
    num = F_k - F_ki;
    den = z_k - z_ki;

    r = rho_k * num.norm()/den.norm();

    return r;
}


void collision::projector(VectorXd &z)
{
    int c = 0;
    vector<candidateP>::iterator it;
    for (it = possibleC.begin(); it != possibleC.end(); it++)
    {
        n = it->n;
        t1= it->t1;
        t2 = it->t2;

        lam = z.segment<3>(3*c); //for one contact pair

        P_n = MatrixXd::Zero(1,3);
        P_t = MatrixXd::Zero(2,3);

        P_n.block<1,3>(0,0) = n;
        P_t.block<1,3>(0,0) = t1;
        P_t.block<1,3>(1,0) = t2;

        rm =  P_n * lam;
        rn = rm(0); //normal direction force

        rt = P_t * lam; //tangent direction


        if (rn > 0 && rt.norm() <= friction * rn )
        {
        	;
        }
        else if (-rn > 0 && rt.norm() <= 1.0/friction * (-rn))
        {
        	rn = 0;
        	rt = Vector2d(0, 0);
        }
        else
        {
        	rn = 1.0/(1+friction*friction)*(rn + friction * rt.norm());
        	rt = 1.0/(1+friction*friction)*(rn + friction * rt.norm())*friction*rt/rt.norm();
        }

        z.segment<3>(3*c) = n * rn  + P_t.transpose() * rt;
        c = c +1;
    }
}
