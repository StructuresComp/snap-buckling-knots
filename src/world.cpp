#include "world.h"
#include <sstream>
#include <iomanip>

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean

	// Physical parameters
	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter
  gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
  maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
	rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter
	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s
	pulltime = m_inputData.GetScalarOpt("pulltime");  //get time of pulling
  mu = m_inputData.GetScalarOpt("friction"); // get coefficient of friction
  r = m_inputData.GetScalarOpt("r");
  statictime = m_inputData.GetScalarOpt("statictime");
  v_constant = m_inputData.GetScalarOpt("v_constant");
  knot_config = m_inputData.GetStringOpt("filename");

	data_resolution = m_inputData.GetScalarOpt("data_resolution");

	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus

	// Viscous drag coefficients using Resistive Force Theory
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile, string file_type) {
    int systemRet = system("mkdir datafiles"); //make the directory
    if (systemRet == -1) {
        cout << "Error in creating directory\n";
    }

    // Open an input file named after the current time
    ostringstream file_name;
    file_name.precision(6);
    file_name << "datafiles/" << file_type;
    file_name << "_" << knot_config;
    file_name << "_dt_" << deltaTime;
    file_name << "_totalTime_" << totalTime;
    file_name << "_mu_" << mu;
    file_name << ".txt";
    outfile.open(file_name.str().c_str());
    outfile.precision(10);
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false)
		return;

	outfile.close();
}

void world::CoutDataC(ofstream &outfile)
{
	if (saveData==false || currentTime < statictime)
		return;

	if (currentTime - statictime > recordTime)
	{
		recordTime += data_resolution;
	}
	else{
		return;
	}

	for (int i = 0; i < rod->nv; i++)
	{
		if (i<rod->ne)
		{
			outfile  << rod->x(4*i) << " " <<rod->
		x(4*i+1) << " " << rod->x(4*i+2) <<" "<< rod->x(4*i+3)<<endl;
		}
		else
		{
			outfile<< rod->x(4*i) << " " <<rod->
		x(4*i+1) << " " << rod->x(4*i+2) <<" "<< 0<<endl;
		}
	}
}


void world::CoutData(ofstream &outfile)
{
	if (saveData==false || currentTime < statictime)
		return;

	if (currentTime - statictime > recordTime)
	{
		recordTime += data_resolution;
	}
	else{
		return;
	}
	double f, f1;
	Vector3d temp;
	temp[0] =stepper->force[0] + stepper->force[4];
	temp[1] = stepper->force[1] + stepper->force[5];
	temp[2] = stepper->force[2] + stepper->force[6];
	f = temp.norm();

	temp[0] = stepper->force[rod->ndof-3] + stepper->force[rod->ndof-7];
	temp[1] = stepper->force[rod->ndof-2] + stepper->force[rod->ndof-6];
	temp[2] = stepper->force[rod->ndof-1] + stepper->force[rod->ndof-5];
	f1 = temp1.norm();
	temp = rod->getVertex(numVertices-1)-rod->getVertex(0);
	outfile<<currentTime<<" "<< f<<" "<<f1<<" "<<temp.norm()<<endl;

}


void world::setRodStepper()
{
	// Set up geometry
	rodGeometry();

	// Create the rod
	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength, theta);

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;

	// Set up boundary condition
	rodBoundaryCondition();

	// setup the rod so that all the relevant variables are populated
	rod->setup();
	// End of rod setup

	// set up the time stepper
	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);
	m_bendingForce = new elasticBendingForce(*rod, *stepper);
	m_twistingForce = new elasticTwistingForce(*rod, *stepper);
	m_inertialForce = new inertialForce(*rod, *stepper);
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity);

  check = new collision(*rod, *stepper, mu, r);

	// Allocate every thing to prepare for the first iteration
	rod->updateTimeStep();

	currentTime = 0.0;
	recordTime = 0;
}

// Setup geometry
void world::rodGeometry()
{
	vertices = MatrixXd(numVertices, 3);

	ifstream myfile(("knot_configurations/" + knot_config).c_str());

	int row1 =  numVertices;
	MatrixXd data = MatrixXd(row1, 4);
  double a ;
	if (myfile.is_open())
	{
		for (int i = 0; i<row1* 4; i++)
		{
			myfile>>a;
			// cout<<a<<endl;
			if (i%4 == 0)
				data(i/4, 0) = a;
			else if (i%4 == 1)
				data(i/4, 1) = a;
			else if (i%4 == 2)
				data(i/4, 2) = a;
			else if (i%4 == 3)
				data(i/4, 3) = a;
		}
	}
  theta = VectorXd::Zero(numVertices - 1);
	for (int i = 0; i< numVertices; i++)
	{
		vertices(i, 0) = data(i,0);
		vertices(i, 1) = data(i,1);
		vertices(i, 2) = data(i,2);
	}


	theta = VectorXd::Zero(numVertices - 1);

}


void world::rodBoundaryCondition()
{
	rod->setVertexBoundaryCondition(rod->getVertex(0),0);
	rod->setVertexBoundaryCondition(rod->getVertex(1),1);
	rod->setThetaBoundaryCondition(rod->getTheta(0), 0);


	rod->setVertexBoundaryCondition(rod->getVertex(numVertices-1),numVertices-1);
  rod->setVertexBoundaryCondition(rod->getVertex(numVertices-2),numVertices-2);
  rod->setThetaBoundaryCondition(rod->getTheta(numVertices-2), numVertices-2);
}


void world::updateBoundary()
{
	Vector3d u;
	u(0) = 0;
	u(1)  = v_constant;
	u(2) = 0;

	if (currentTime>statictime && currentTime<=statictime + pulltime)
	{
		rod->setVertexBoundaryCondition(rod->getVertex(0)-u*deltaTime,0);
		rod->setVertexBoundaryCondition(rod->getVertex(1)-u*deltaTime,1);

		rod->setVertexBoundaryCondition(rod->getVertex(numVertices-1)+u*deltaTime,numVertices-1);
		rod->setVertexBoundaryCondition(rod->getVertex(numVertices-2)+u*deltaTime,numVertices-2);
	}

	if (currentTime>statictime+pulltime)
	{
		currentTime = totalTime;
	}

}

void world::updatesolver(int num, bool flag)
{
	if (flag)
		stepper->update(num);
	else
		stepper->updateforeigen(num);
    totalForce = stepper->getForce();
}



void world::updateTimeStep()
{
	bool solved = false;
	iter = 0;

	// Start with a trial solution for our solution x
	rod->updateGuess(); // x = x0 + u * dt
	updateBoundary();
	updatesolver(rod->uncons, false);
	newtonMethod(solved, false);
	check->setZero();
	bool flag = check->candidateTest();

	if (flag)
	{
		rod->initialguess();
		intializeR();
	  check->possibleCandidate();//define some parameters
		check->prepareForSolve();
		check->initialize();
		updatesolver(3*check->p_pairs, true);
		//solve for lamda
	  check->initialguess();
	  solved = false;
	  newtonMethodFac(solved);
	  updatesolver(rod->uncons, false);
	  solved = false;
	  rod->updateGuess();
	  updateBoundary();
	  newtonMethod(solved, true);
	}
	rod->updateTimeStep();

	if (render) cout << "time: " << currentTime << " iter=" << iter <<" converged force: "<<force_record<<endl;

	currentTime += deltaTime;

	if (solved == false)
	{
		currentTime = totalTime; // we are exiting
	}
}

void world::intializeR()
{
	// compute Force
	rod->prepareForIteration();

	stepper->setZero();

	m_stretchForce->computeFs();
	m_bendingForce->computeFb();
	m_twistingForce->computeFt();
	m_gravityForce->computeFg();
	m_dampingForce->computeFd();
	stepper->setFi();

	stepper->setZero();
  m_inertialForce->computeFi();
	m_stretchForce->computeFs();
	m_bendingForce->computeFb();
	m_twistingForce->computeFt();
	m_gravityForce->computeFg();
	m_dampingForce->computeFd();

	stepper->setFe();
}




void world::newtonMethodFac(bool &solved)
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	double a = 10;
	iter = 0;

	double tol_c = check->lamda.norm()*1e-4;
	while (solved == false)
	{
		check->prepareForIteration();
		stepper->setZero();
		// Compute the forces and the jacobians
		check->computeFacandJac();

		// Compute norm of the force equations.
		normf = 0;
		for (int i=0; i < stepper->freeDOF; i++)
		{
			normf += stepper->Force[i] * stepper->Force[i];
		}
		normf = sqrt(normf);

		if (iter == 0)
		{
			normf0 = normf;
		}
		if (normf <= 1e-12 )
		{
			solved = true;
		}
		else if(iter > 0 && normf <= 1e-12 )
		{
			solved = true;
		}
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			if (iter < 10)
			{
				a = 1;
			}
			else
			{
				a = linesearch();
			}
			check->updateForce(stepper->DX, a);
			iter++;
		}

		if (iter > 100 && normf<=1e-5)
		{
			solved = true;
		}
		if (iter > maxIter)
		{
			cout << "Error. Fac could not converge. Exiting.\n";
			break;
		}
	}

	force_record = normf;
}

void world::newtonMethod(bool &solved, bool contact_flag)
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	iter = 0;
	while (solved == false)
	{
		rod->prepareForIteration();

		stepper->setZero();

		// Compute the forces and the jacobians
		m_inertialForce->computeFi();
		m_inertialForce->computeJi();

		m_stretchForce->computeFs();
		m_stretchForce->computeJs();

		m_bendingForce->computeFb();
		m_bendingForce->computeJb();

		m_twistingForce->computeFt();
		m_twistingForce->computeJt();

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();

		m_dampingForce->computeFd();
		m_dampingForce->computeJd();

		if (contact_flag)
		{
			check->computeFc();
		}


		// Compute norm of the force equations.
		normf = 0;
		for (int i=0; i < stepper->freeDOF; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}
		normf = sqrt(normf);
		if (iter == 0)
		{
			normf0 = normf;
		}

		if (normf <= forceTol )
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		if (solved == false)
		{
			stepper->integrator_eigen(); // Solve equations of motion
			rod->updateNewtonX(totalForce, 1); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}
}


int world::simulationRunning()
{
	if (currentTime<totalTime)
		return 1;
	else
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i)
{
	return rod->x[i] / (0.6*RodLength);
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}


double world::linesearch()
{
	check->lamda_old = check->lamda;

	double amax = 2;
  double amin = 1e-2;
  double al = 0;
  double au = 1;

  double a = 1;

	double q0 = 0.5 * pow(stepper->Force.norm(), 2);
  if (q0 == 0) return 1;

	double dq0 = -(stepper->Force.transpose() * stepper->Jacobian * stepper->DX)(0);

	bool success = false;
	double m2 = 0.9;
  double m1 = 0.1;

	while (!success)
  {
		stepper->setZero();
		check->lamda = check->lamda_old;
		check->updateForce(stepper->DX, a);
		check->computeFac();
		double q = 0.5 * pow(stepper->Force.norm(), 2);

    double slope = (q - q0)/a;

		if (slope >= m2 * dq0 && slope <= m1 * dq0)
    {
      success = true;
    }
    else
    {
      if (slope < m2 * dq0)
      {
        al = a;
      }
      else
      {
        au = a;
      }

      if (au < amax)
      {
        a = 0.5 * (al + au);
      }
      else
      {
        a = 10 * a;
      }
    }
    if (a > amax || a < amin)
    {
      break;
    }
  }

	check->lamda = check->lamda_old;

	return a;
}
