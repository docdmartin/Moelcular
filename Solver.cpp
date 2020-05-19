#include "Solver.h"

using namespace std;

Solver::Solver( size_t n, vector<Matrix3x3>& hessian, vector<double> frequency, double* electric_potential, Kernel& kernel ):mHessian(hessian), mKernel(kernel)
{
  mLength            = n;
  mFrequencyVec      = frequency;
  mFrequency         = mFrequencyVec[0];
  mFreqSq            = mFrequency * mFrequency;
  mElectricPotential = electric_potential;

  r_hat = new double[mLength]();
  r 	  = new double[mLength]();
  v 	  = new double[mLength]();
  p 	  = new double[mLength]();
  h 	  = new double[mLength]();
  s 	  = new double[mLength]();
  t 	  = new double[mLength]();

  intermediate_vec = new double[mLength]();
  real_response    = new double[mLength]();
  imag_response    = new double[mLength]();
  mb               = new double[mLength]();
  mAb              = new double[mLength]();

  r_hat_end = &r_hat[mLength - 1];
  t_end     = &t    [mLength - 1];
  v_end     = &v    [mLength - 1];
}

Solver::~Solver()
{
  delete[] r_hat;
  delete[] r;
  delete[] v;
  delete[] p;
  delete[] h;
  delete[] s;
  delete[] t;
  delete[] intermediate_vec;
  delete[] real_response;
  delete[] imag_response;
  delete[] mb;
  delete[] mAb;
}

void Solver::Print(){
  /*
  for( Matrix3x3& hess : mHessian ){
    hess.Print();
  }
  */

  //double* test_mult = new double[mLength]();
  //for( Matrix3x3& hess : mHessian ){
  //  hess.Multiply( test_mult, mElectricPotential );
  //}
  for( size_t cnt = 0; cnt < mLength; ++cnt ){
      cout << " " << mElectricPotential[cnt];
      if( cnt%3 == 2 )
      cout << endl;
  }

  //delete[] test_mult;
}


vector<pair<double, double>> Solver::CalculateResponse( double weight )
{
  solve_response = real_response;
  for( size_t cnt = 0; cnt < mLength; ++cnt )
    mAb[cnt] = 0.0;
  for( Matrix3x3& hess : mHessian ){
    hess.Multiply( mAb, mElectricPotential );
  }

  mConvergenceThreshold = 0.0;
  for( size_t cnt = 0; cnt < mLength; ++cnt ){
    mb[cnt] = mAb[cnt];
    mConvergenceThreshold += mb[cnt] * mb[cnt];
  }
  mConvergenceThreshold *= 1e-12;

  vector<pair<double, double>> response;

  double r_response, i_response;
  for( double frq : mFrequencyVec ){
/*
    for( size_t cnt = 0; cnt < mLength; ++cnt ){
      r_hat[cnt] = 0.0;
      r 	 [cnt] = 0.0;
      v 	 [cnt] = 0.0;
      p 	 [cnt] = 0.0;
      h 	 [cnt] = 0.0;
      s 	 [cnt] = 0.0;
      t 	 [cnt] = 0.0;

      intermediate_vec[cnt] = 0.0;
      imag_response   [cnt] = 0.0;
    }
*/

    mFrequency = frq        * weight;
    mFreqSq    = mFrequency * mFrequency;

    Solve();

    if( mFrequency >= 1.0 ){
      for( size_t cnt = 0; cnt < mLength; ++cnt )
        imag_response[cnt] = 0.0;
      for( Matrix3x3& hess : mHessian ){
        hess.Multiply( imag_response, real_response );
      }
      for( size_t cnt = 0; cnt < mLength; ++cnt )
        imag_response[cnt] -= mElectricPotential[cnt];
    }
    else{
      solve_response = imag_response;

      double tmp_threshold = mConvergenceThreshold;
      mConvergenceThreshold = 0.0;
      for( size_t cnt = 0; cnt < mLength; ++cnt ){
        mb[cnt] = -mFrequency * mElectricPotential[cnt];
        mConvergenceThreshold += mb[cnt] * mb[cnt];
      }
      mConvergenceThreshold *= 1e-6;

      Solve();

      for( size_t cnt = 0; cnt < mLength; ++cnt ){
        mb[cnt] = mAb[cnt];
      }
      solve_response        = real_response;
      mConvergenceThreshold = tmp_threshold;
    }

    r_response = 0.0;
    i_response = 0.0;
    for( size_t cnt = 0; cnt < mLength; ++cnt ){
      r_response += mElectricPotential[cnt] * real_response[cnt];
      i_response += mElectricPotential[cnt] * imag_response[cnt];
    }
    if( mFrequency >= 1.0 )
      i_response /= mFrequency;

    response.push_back( pair<double, double>(r_response, i_response) );
  }

  return response;
}

void Solver::Reset(){
  for( size_t cnt = 0; cnt < mLength; ++cnt ){
    r_hat[cnt] = 0.0;
    r 	 [cnt] = 0.0;
    v 	 [cnt] = 0.0;
    p 	 [cnt] = 0.0;
    h 	 [cnt] = 0.0;
    s 	 [cnt] = 0.0;
    t 	 [cnt] = 0.0;

    intermediate_vec[cnt] = 0.0;
    real_response   [cnt] = 0.0;
    imag_response   [cnt] = 0.0;
    mb              [cnt] = 0.0;
  }
}

bool Solver::Solve(){
  (void) mKernel;

  RealMatrixMultiply( solve_response, v );

	threshold = 0.0;
	for (size_t cnt = 0; cnt < mLength; ++cnt) {
			r    [cnt] = mb[cnt] - v[cnt];
			r_hat[cnt] = r[cnt];
			v    [cnt] = 0.;
			p    [cnt] = 0.;
      h 	 [cnt] = 0.0;
      s 	 [cnt] = 0.0;
      t 	 [cnt] = 0.0;

  //    imag_response[cnt] = 0.0;

			threshold += r[cnt] * r[cnt];
	}

	if(mConvergenceThreshold == 0.0){
		for( size_t cnt = 0; cnt < mLength; ++cnt )
			solve_response[cnt] = 0.0;

		return true;
	}

	alpha = 1.0;
	omega = 1.0;
	rho   = 1.0;

	for (size_t curr_iter = 0; curr_iter < mLength+1000; ++curr_iter) {
		iteration();

//    mKernel.RemoveProjection( r );
//    mKernel.RemoveProjection( solve_response );

    if( curr_iter < 10 )
      continue;

		// Check to see if this has completed
		if (threshold < mConvergenceThreshold){
//			cout << "Converged to " << threshold << " (< " << mConvergenceThreshold << ") in " << curr_iter+1 << " iterations" << endl;

			return true;
		}
	}

	cout << "Convergence fail " << threshold << " not less than " << mConvergenceThreshold << " after " <<  mLength << " iterations"<< endl;

	return false;
}


void Solver::RealMatrixMultiply( double* input, double* output ){

  for( size_t cnt = 0; cnt < mLength; ++cnt ){
    intermediate_vec[cnt] = 0.0;
    output[cnt] = mFreqSq * input[cnt];
  }
  for( Matrix3x3& hess : mHessian ){
    hess.Multiply( intermediate_vec, input );
  }
  for( Matrix3x3& hess : mHessian ){
    hess.Multiply( output, intermediate_vec );
  }
}



void Solver::iteration(){

	beta = omega * rho;
	if( beta == 0.0 ){ // protect against division by zero
		return;//throw string("Beta fail");
	}
	a   = &r[0];
	b   = &r_hat[0];
	rho = (*a) * (*b);
	do{
			rho += *(++a) * *(++b);
	}while(b != r_hat_end);
	beta = (alpha * rho) / beta;

	a  = &r[0];
	b  = &p[0];
	c  = &v[0];
	*b = *a + beta * (*b - omega * *c);
	do{
			++b;
			*b = *(++a) + beta * (*b - omega * *(++c));
	}while(c != v_end);

  RealMatrixMultiply( p, v );

	a     = &v[0];
	b     = &r_hat[0];
	alpha = (*a) * (*b);
	do{
			alpha += *(++a) * *(++b);
	}while(b != r_hat_end);
	if( alpha == 0.0 ){ // protect against division by zero
		return;//throw string("Alpha fail");
	}
	alpha = rho / alpha;

	a     = &h[0];
	b     = &p[0];
	c     = &v[0];
	d     = &solve_response[0];
	g     = &s[0];
	f     = &r[0];
	*a = *d + alpha * *b;
	*g = *f - alpha * *c;
	do{
		*(++a) = *(++d) + alpha * *(++b);  // a and d are solutions
		*(++g) = *(++f) - alpha * *(++c);  // f and e are residuals
	}while(c != v_end);

	RealMatrixMultiply( s, t );

	a  = &s[0];
	b  = &t[0];
	ts = (*a) * (*b);
	tt = (*b) * (*b);
	do{
				ts += *(++a) * *(++b);
				tt +=   (*b) *   (*b);
	}while(b != t_end);
	if( tt == 0.0 ){ // protect against division by zero
		return;//throw string("Omega fail");
	}
	omega = ts / tt;

	a  = &h[0];
	b  = &t[0];
	c  = &r[0];
	d  = &solve_response[0];
	g  = &s[0];
	*d = *a + omega * *g;
	*c = *g - omega * *b;
	threshold = (*c) * (*c);
	do{
		*(++d) = *(++a) + omega * *(++g);
		*(++c) = *g     - omega * *(++b);
		threshold += (*c) * (*c);
	}while(b != t_end);
}
