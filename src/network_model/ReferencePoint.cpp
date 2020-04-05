#include <cmath>
#include <fstream>

#include "network_model/ReferencePoint.h"

ReferencePoint::ReferencePoint(double x, double y, double z, vector<Node> &nodes, HessianMatrix &hessianMatrix) {
    mPosition.push_back(x);
    mPosition.push_back(y);
    mPosition.push_back(z);

    mPotentialL2    = 0.0;
    mPotentialSigma = 0.0;

    mTxxL2    = 0.0;
    mTxxSigma = 0.0;

    mTyyL2    = 0.0;
    mTyySigma = 0.0;

    mTzzL2    = 0.0;
    mTzzSigma = 0.0;

    preCompute( nodes, hessianMatrix );
}

ReferencePoint::ReferencePoint(const ReferencePoint& rp) {
    mPosition        = rp.GetPosition();

    mPotentialL2     = rp.GetPotentialL2();
    mPotentialSigma  = rp.GetPotentialSigma();
    mPotentialLambda = rp.GetPotentialLambda();

    mTxxL2     = rp.GetTxxL2();
    mTxxSigma  = rp.GetTxxSigma();
    mTxxLambda = rp.GetTxxLambda();

    mTyyL2     = rp.GetTyyL2();
    mTyySigma  = rp.GetTyySigma();
    mTyyLambda = rp.GetTyyLambda();

    mTzzL2     = rp.GetTzzL2();
    mTzzSigma  = rp.GetTzzSigma();
    mTzzLambda = rp.GetTzzLambda();
}

ReferencePoint::~ReferencePoint()
{
}


void ReferencePoint::preCompute(vector<Node> &nodes, HessianMatrix &hessianMatrix) {
    vector<double> electricPotential;
    vector<double> dipoleTensorX;
    vector<double> dipoleTensorY;
    vector<double> dipoleTensorZ;

    electricPotential.reserve( 3 * nodes.size() );
    dipoleTensorX    .reserve( 3 * nodes.size() );
    dipoleTensorY    .reserve( 3 * nodes.size() );
    dipoleTensorZ    .reserve( 3 * nodes.size() );

    for(auto node : nodes) {
        vector<double> ep = node.GetElectricPotential( mPosition );
        electricPotential.push_back( ep[0] );
        electricPotential.push_back( ep[1] );
        electricPotential.push_back( ep[2] );

        //double q = node.GetQ();

        double dx = mPosition[0] - node.GetX();
        double dy = mPosition[1] - node.GetY();
        double dz = mPosition[2] - node.GetZ();

        if(dx == 0.0 && dy == 0.0 && dz == 0.0){
          cout << node.GetID() << " is the probe point" << endl;
          // electricPotential.push_back( 0.0 );
          // electricPotential.push_back( 0.0 );
          // electricPotential.push_back( 0.0 );

          dipoleTensorX.push_back( 0.0 );
          dipoleTensorX.push_back( 0.0 );
          dipoleTensorX.push_back( 0.0 );

          dipoleTensorY.push_back( 0.0 );
          dipoleTensorY.push_back( 0.0 );
          dipoleTensorY.push_back( 0.0 );

          dipoleTensorZ.push_back( 0.0 );
          dipoleTensorZ.push_back( 0.0 );
          dipoleTensorZ.push_back( 0.0 );

          continue;
        }

        double dx_sq = dx*dx;
        double dy_sq = dy*dy;
        double dz_sq = dz*dz;
        double r_sq  = dx_sq + dy_sq + dz_sq;
        double r     = sqrt( r_sq );

        double dipole_constant    = r * r_sq;
    //    double potential_constant = q   / dipole_constant;
        dipole_constant           = dipole_constant * r_sq;
        double diag_constant      = r_sq / dipole_constant;
        dipole_constant           = -3.0 / dipole_constant;
//cout << "Parts of potential: " << dx * potential_constant << " " << dy * potential_constant << " " << dz * potential_constant << endl;
        // electricPotential.push_back( dx * potential_constant );
        // electricPotential.push_back( dy * potential_constant );
        // electricPotential.push_back( dz * potential_constant );

        dipoleTensorX.push_back( diag_constant + dx_sq * dipole_constant );
        double yy              = diag_constant + dy_sq * dipole_constant;
        double zz              = diag_constant + dz_sq * dipole_constant;

        double tmp_xy = dx * dy * dipole_constant;
        dipole_constant *= dz;
        double tmp_xz = dx * dipole_constant;
        double tmp_yz = dy * dipole_constant;

        dipoleTensorX.push_back( tmp_xy );
        dipoleTensorX.push_back( tmp_xz );

        dipoleTensorY.push_back( tmp_xy );
        dipoleTensorY.push_back( yy );
        dipoleTensorY.push_back( tmp_yz );

        dipoleTensorZ.push_back( tmp_xz );
        dipoleTensorZ.push_back( tmp_yz );
        dipoleTensorZ.push_back( zz );
    }

    fstream fs;
    fs.open ("ElectricPotential.txt", std::fstream::out );

    for(int cnt = 0; cnt < static_cast<int>(electricPotential.size()); cnt+=3)
      fs << electricPotential[cnt] << " " << electricPotential[cnt+1] << " " << electricPotential[cnt+2] << endl;

    fs.close();

    hessianMatrix.RemoveProjection( electricPotential );

    calculateFrequencyTerms( electricPotential, mPotentialL2, mPotentialSigma, mPotentialLambda, nodes, hessianMatrix );
//    calculateFrequencyTerms( dipoleTensorX    , mTxxL2      , mTxxSigma      , mTxxLambda      , nodes, hessianMatrix );
//    calculateFrequencyTerms( dipoleTensorY    , mTyyL2      , mTyySigma      , mTyyLambda      , nodes, hessianMatrix );
//    calculateFrequencyTerms( dipoleTensorZ    , mTzzL2      , mTzzSigma      , mTzzLambda      , nodes, hessianMatrix );
}


void ReferencePoint::calculateFrequencyTerms(
  vector<double>& b_vector,
  double        & bT_b, 
  double        & sigma,
  vector<double>& lambda,
  vector<Node>  & nodes,
  HessianMatrix & hessianMatrix ) {

  vector<double> input_vec;
  vector<double> output_vec;
  input_vec .resize( b_vector.size() );
  output_vec.resize( b_vector.size() );

  bT_b = innerProduct(b_vector, b_vector);
  double mu = sqrt( bT_b );
  vector<double> u_hat;
  u_hat.resize( b_vector.size() );
  for(int cnt = 0; cnt < static_cast<int>(b_vector.size()); ++cnt){
    u_hat[cnt] = b_vector[cnt] / mu;

    input_vec [cnt] = u_hat[cnt];
    output_vec[cnt] = 0.0;
  }


  for(int cnt = 0; cnt < 100; ++cnt){
    hessianMatrix.MultiplyMatrix(input_vec, output_vec);

    double tmp_lambda = 0.0;
    for(int cnt = 0; cnt < static_cast<int>(b_vector.size()); ++cnt){
      tmp_lambda += output_vec[cnt] * u_hat[cnt];
    }
    double tmp_sigma = 0.0;
    for(int cnt = 0; cnt < static_cast<int>(b_vector.size()); ++cnt){
      input_vec[cnt] = output_vec[cnt] - tmp_lambda * u_hat[cnt];
      output_vec[cnt] = 0.0;
      tmp_sigma += input_vec[cnt] * input_vec[cnt];
    }
    tmp_sigma = sqrt( tmp_sigma );
    double unit_test = 0.0;
    for(int cnt = 0; cnt < static_cast<int>(b_vector.size()); ++cnt){
      input_vec[cnt] /= tmp_sigma;
      unit_test += input_vec[cnt] * input_vec[cnt];
    }

    cout << "Lambda = " << tmp_lambda << ", sigma = " << tmp_sigma << ", unit test: " << unit_test << endl;;
  }







  vector<double> result;
  result.resize( 3 * nodes.size() );

  vector<double> result2;
  result2.resize( 3 * nodes.size() );

  hessianMatrix.MultiplyMatrix(b_vector, result);
  lambda.push_back( innerProduct(b_vector, result) / bT_b );

  double last_lambda[3] = {0.0, 0.0, 0.0};
  double scalar_value;

  int cnt      = 0;
  int term_cnt = 0;
  while( term_cnt < 2         && cnt < 20 ) {
      last_lambda[2] = *(lambda.rbegin());
//cout << result[0] << " -> ";
      // Two steps per loop to address ping-pong of vectors
      grahamSchmit( result, b_vector, *(lambda.rbegin()) );
//  cout << result[0] << " -> ";
      scalar_value = rescale( result );
//  cout << result[0] << " -> ";
      hessianMatrix.MultiplyMatrix(result, result2);
//  cout << result[0] << " -> " << result2[0] << endl;
//  cout << cnt << " : " << scalar_value << ", " << innerProduct(b_vector, result2) << endl;
      lambda.push_back( scalar_value * innerProduct(b_vector, result2) / bT_b );
      last_lambda[1] = *(lambda.rbegin());
//cout << result2[0] << " -> ";
      grahamSchmit( result2, b_vector, *(lambda.rbegin()) );
//  cout << result2[0] << " -> ";
      scalar_value = rescale( result2 );
//  cout << result2[0] << " -> ";
      hessianMatrix.MultiplyMatrix(result2, result);
//  cout << result2[0] << " -> " << result[0] << endl;
//  cout << cnt+1 << " : " << scalar_value << ", " << innerProduct(b_vector, result) << endl;
      lambda.push_back( scalar_value * innerProduct(b_vector, result) / bT_b );
      last_lambda[0] = *(lambda.rbegin());

//      cout << "\tLambdas: " << last_lambda[0] << " " << last_lambda[1] << " " << last_lambda[2] << endl;

      cnt += 2;
      if(cnt < 20)
          continue;

      if(fabs( (last_lambda[1]*last_lambda[1])/(last_lambda[0]*last_lambda[2]) - 1.0) < 1e-6) {
          term_cnt += 1;
      }
      else{
          term_cnt = 0;
      }
  }
//  cout << "Iterative convergence was achieved in " << cnt << " iterations" << endl;
//  cout << "\tLambdas: " << last_lambda[0] << " " << last_lambda[1] << " " << last_lambda[2] << endl;

  sigma = last_lambda[0] / last_lambda[1];

}


double ReferencePoint::innerProduct( vector<double>& a, vector<double>& b ) {
    double ip = 0.0;
    for(int cnt = 0; cnt < static_cast<int>(a.size()); ++cnt) {
        ip += a[cnt] * b[cnt];
    }

    return ip;
}

void ReferencePoint::grahamSchmit( vector<double>& a, vector<double>& b, double x ) {

    for(int cnt = 0; cnt < static_cast<int>(a.size()); ++cnt) {
        a[cnt] -= x * b[cnt];
    }
}

double ReferencePoint::rescale( vector<double>& a ) {
    double max_value = 0.0;

    for(int cnt = 0; cnt < static_cast<int>(a.size()); ++cnt) {
        if(max_value < fabs(a[cnt]))
            max_value = fabs(a[cnt]);
    }

    for(int cnt = 0; cnt < static_cast<int>(a.size()); ++cnt) {
        a[cnt] /= max_value;
    }

    return max_value;
}

vector< pair<double, double> > ReferencePoint::SingleModeFrequencyResponse(double omega) {
/*
    cout << "Frequency: " << omega << endl;
    cout << "Electric potential: L2 = " << mPotentialL2 << ", sigma = " << mPotentialSigma << ", num of lambdas = " << mPotentialLambda.size() << endl;
    cout << "Txx      potential: L2 = " << mTxxL2       << ", sigma = " << mTxxSigma       << ", num of lambdas = " << mTxxLambda.size() << endl;
    cout << "Tyy      potential: L2 = " << mTyyL2       << ", sigma = " << mTyySigma       << ", num of lambdas = " << mTyyLambda.size() << endl;
    cout << "Tzz      potential: L2 = " << mTzzL2       << ", sigma = " << mTzzSigma       << ", num of lambdas = " << mTzzLambda.size() << endl;
*/

    pair<double, double> electric_potential = singleModeFrequencyResponse( omega, mPotentialL2, mPotentialSigma, mPotentialLambda );

    pair<double, double> Txx_potential = singleModeFrequencyResponse( omega, mTxxL2, mTxxSigma, mTxxLambda );
    pair<double, double> Tyy_potential = singleModeFrequencyResponse( omega, mTyyL2, mTyySigma, mTyyLambda );
    pair<double, double> Tzz_potential = singleModeFrequencyResponse( omega, mTzzL2, mTzzSigma, mTzzLambda );

    Txx_potential.first  += Tyy_potential.first  + Tzz_potential.first ;
    Txx_potential.second += Tyy_potential.second + Tzz_potential.second;

    vector< pair<double, double> > complex_potentials;
    complex_potentials.push_back( electric_potential );
    complex_potentials.push_back( Txx_potential      );

    return complex_potentials;
}

vector< pair<double, double> > ReferencePoint::DualModeFrequencyResponse(double omega1, double omega2, double mix_ratio) {

  pair<double, double> electric_potential = dualModeFrequencyResponse( omega1, omega2, mix_ratio, mPotentialL2, mPotentialSigma, mPotentialLambda );

  pair<double, double> Txx_potential = dualModeFrequencyResponse( omega1, omega2, mix_ratio, mTxxL2, mTxxSigma, mTxxLambda );
  pair<double, double> Tyy_potential = dualModeFrequencyResponse( omega1, omega2, mix_ratio, mTyyL2, mTyySigma, mTyyLambda );
  pair<double, double> Tzz_potential = dualModeFrequencyResponse( omega1, omega2, mix_ratio, mTzzL2, mTzzSigma, mTzzLambda );

  Txx_potential.first  += Tyy_potential.first  + Tzz_potential.first ;
  Txx_potential.second += Tyy_potential.second + Tzz_potential.second;

  vector< pair<double, double> > complex_potentials;
  complex_potentials.push_back( electric_potential );
  complex_potentials.push_back( Txx_potential      );

  return complex_potentials;
}





pair<double, double> ReferencePoint::singleModeFrequencyResponse(double omega, double l2, double sigma, vector<double>& lambda) {
    pair<double, double> complex_response(0.0, 0.0);

    vector<double> four_by_four(16, 0.0);
    /* 0 4  8 12
       1 5  9 13
       2 6 10 14
       3 7 11 15 */

    vector<double> solution_vector(4, 0.0);
    double  param_1_value = 1.0;
    double  param_2_value = 1.0;
    double *r_ptr = &param_1_value;

    // initialize
    four_by_four[ 2] =  sigma;             four_by_four[ 3] = omega; four_by_four[ 6] = -four_by_four[ 3]; four_by_four[ 7] = four_by_four[ 2];
    four_by_four[ 8] =  lambda[0];         four_by_four[ 9] = omega; four_by_four[12] = -four_by_four[ 9]; four_by_four[13] = four_by_four[ 8];
    four_by_four[10] = *(lambda.rbegin());                                                                 four_by_four[15] = four_by_four[10];

    solution_vector[0] = 1.0;

    // iterate
    double tmp_value;
    for(int cnt = 0; cnt < static_cast<int>(lambda.size() - 1); ++cnt) {
        tmp_value = lambda[cnt+1] / omega;
        if( r_ptr == &param_1_value ) {
            four_by_four[12] -= tmp_value * param_1_value;
            four_by_four[ 9] += tmp_value * param_2_value;

            param_1_value /= -omega;
            param_2_value /=  omega;

            r_ptr = &param_2_value;
        }
        else {
          four_by_four[ 8] += tmp_value * param_1_value;
          four_by_four[13] -= tmp_value * param_2_value;

          param_1_value /=  omega;
          param_2_value /= -omega;

          r_ptr = &param_1_value;
        }
    }
    if( r_ptr == &param_1_value ) {
        four_by_four[0] = param_1_value;
        four_by_four[5] = param_2_value;
    }
    else {
      four_by_four[4] = param_1_value;
      four_by_four[1] = param_2_value;
    }

    // solve for alpha, beta
    double alpha = 0.0;
    double beta  = 0.0;
    partialSolver(four_by_four, solution_vector, l2, &alpha, &beta);

    complex_response.first  = alpha;
    complex_response.second =  beta;
    return complex_response;
}
pair<double, double> ReferencePoint::dualModeFrequencyResponse(double omega1, double omega2, double mix_ratio, double l2, double sigma, vector<double>& lambda) {
    pair<double, double> complex_response(0.0, 0.0);

    double omega_plus = omega1 + omega2;
    double omega_prod = omega1 * omega2;
    double omega_bar  = (1 - mix_ratio) * omega1 + mix_ratio * omega2;

    vector<double> six_by_six(36, 0.0);
    /* 0  6 12 18 24 30
       1  7 13 19 25 31
       2  8 14 20 26 32
       3  9 15 21 27 33
       4 10 16 22 28 34
       5 11 17 23 29 35 */

    vector<double> solution_vector(6, 0.0);
    vector<double> r_0(4, 0.0);
    vector<double> r_1(4, 0.0);
    vector<double> r_2(4, 0.0);
    vector<double> s_0(4, 0.0);
    vector<double> s_1(4, 0.0);
    vector<double> s_2(4, 0.0);

    int lowest_index;

    // initialize
    six_by_six[ 4] = sigma*sigma - omega_prod;                 six_by_six[ 5] = omega_plus*sigma;              six_by_six[10] = -six_by_six[ 5]; six_by_six[11] = six_by_six[ 4];
    six_by_six[12] = lambda[0];                                six_by_six[13] = omega_plus;                    six_by_six[18] = -six_by_six[13]; six_by_six[19] = six_by_six[12];
    six_by_six[14] = lambda[1] - omega_prod;                                                                                                     six_by_six[21] = six_by_six[14];
    six_by_six[16] = *(lambda.rbegin());                                                                                                         six_by_six[23] = six_by_six[16];
    six_by_six[24] = lambda[0]*lambda[0]+lambda[1]-omega_prod; six_by_six[25] = omega_plus*lambda[0];          six_by_six[30] = -six_by_six[25]; six_by_six[31] = six_by_six[24];
    six_by_six[26] = lambda[0]*lambda[1]+lambda[2];            six_by_six[27] = omega_plus*lambda[1];          six_by_six[32] = -six_by_six[27]; six_by_six[33] = six_by_six[26];
    six_by_six[28] = *(lambda.rbegin())*(lambda[0]+sigma);     six_by_six[29] = *(lambda.rbegin())*omega_plus; six_by_six[34] = -six_by_six[29]; six_by_six[35] = six_by_six[28];

    solution_vector[0] = lambda[0];
    solution_vector[1] = omega_bar;
    solution_vector[2] = lambda[1];
    solution_vector[4] = *(lambda.rbegin());

    r_1[0] = 1.0; r_1[3] =  omega_plus;
    s_1[1] = 1.0; s_1[2] = -omega_plus;
    r_2[2] = 1.0;
    s_2[3] = 1.0;

    lowest_index = 1;

    // iterate
    double bb1, bb2, v0, v1, v2;
    v2  = 1.0 / omega_prod;
    v1  = omega_plus / omega_prod;
    for(int cnt = 1; cnt < static_cast<int>(lambda.size() - 1); ++cnt) {
        if(cnt == static_cast<int>(lambda.size() - 2)) {
            bb1 = (lambda[0] + sigma) * lambda[cnt+1] * v2;
        }
        else {
            bb1 = (lambda[0]*lambda[cnt+1] + lambda[cnt+2]) * v2;
        }
        bb2 = lambda[cnt+1] * v1;
        v0  = lambda[cnt+1] * v2;

        if(lowest_index == 0) {
            six_by_six[24] += r_0[0] * bb1; six_by_six[25] += r_0[1] * bb1; six_by_six[26] += r_0[2] * bb1; six_by_six[27] += r_0[3] * bb1;
            six_by_six[30] -= r_0[0] * bb2; six_by_six[31] -= r_0[1] * bb2; six_by_six[32] -= r_0[2] * bb2; six_by_six[33] -= r_0[3] * bb2;
            six_by_six[12] += r_0[0] * v0 ; six_by_six[13] += r_0[1] * v0 ; six_by_six[14] += r_0[2] * v0 ; six_by_six[15] += r_0[3] * v0 ;

            s_1[0] -= r_0[0] * v1; s_1[1] -= r_0[1] * v1; s_1[2] -= r_0[2] * v1; s_1[3] -= r_0[3] * v1;
            r_2[0]  = r_0[0] * v2; r_2[1]  = r_0[1] * v2; r_2[2]  = r_0[2] * v2; r_2[3]  = r_0[3] * v2;

            solution_vector[0] -= r_0[0] * v0 ; solution_vector[1] -= r_0[1] * v0 ; solution_vector[2] -= r_0[2] * v0 ; solution_vector[3] -= r_0[3] * v0 ;

            six_by_six[24] += s_0[0] * bb2; six_by_six[25] += s_0[1] * bb2; six_by_six[26] += s_0[2] * bb2; six_by_six[27] += s_0[3] * bb2;
            six_by_six[30] += s_0[0] * bb1; six_by_six[31] += s_0[1] * bb1; six_by_six[32] += s_0[2] * bb1; six_by_six[33] += s_0[3] * bb1;
            six_by_six[18] += s_0[0] * v0 ; six_by_six[19] += s_0[1] * v0 ; six_by_six[20] += s_0[2] * v0 ; six_by_six[21] += s_0[3] * v0 ;

            r_1[0] += s_0[0] * v1; r_1[1] += s_0[1] * v1; r_1[2] += s_0[2] * v1; r_1[3] += s_0[3] * v1;
            s_2[0]  = s_0[0] * v2; s_2[1]  = s_0[1] * v2; s_2[2]  = s_0[2] * v2; s_2[3]  = s_0[3] * v2;

            lowest_index = 1;
        }
        else if(lowest_index == 1) {
            six_by_six[24] += r_1[0] * bb1; six_by_six[25] += r_1[1] * bb1; six_by_six[26] += r_1[2] * bb1; six_by_six[27] += r_1[3] * bb1;
            six_by_six[30] -= r_1[0] * bb2; six_by_six[31] -= r_1[1] * bb2; six_by_six[32] -= r_1[2] * bb2; six_by_six[33] -= r_1[3] * bb2;
            six_by_six[12] += r_1[0] * v0 ; six_by_six[13] += r_1[1] * v0 ; six_by_six[14] += r_1[2] * v0 ; six_by_six[15] += r_1[3] * v0 ;

            s_2[0] -= r_1[0] * v1; s_2[1] -= r_1[1] * v1; s_2[2] -= r_1[2] * v1; s_2[3] -= r_1[3] * v1;
            r_0[0]  = r_1[0] * v2; r_0[1]  = r_1[1] * v2; r_0[2]  = r_1[2] * v2; r_0[3]  = r_1[3] * v2;

            solution_vector[0] -= r_1[0] * v0 ; solution_vector[1] -= r_1[1] * v0 ; solution_vector[2] -= r_1[2] * v0 ; solution_vector[3] -= r_1[3] * v0 ;

            six_by_six[24] += s_1[0] * bb2; six_by_six[25] += s_1[1] * bb2; six_by_six[26] += s_1[2] * bb2; six_by_six[27] += s_1[3] * bb2;
            six_by_six[30] += s_1[0] * bb1; six_by_six[31] += s_1[1] * bb1; six_by_six[32] += s_1[2] * bb1; six_by_six[33] += s_1[3] * bb1;
            six_by_six[18] += s_1[0] * v0 ; six_by_six[19] += s_1[1] * v0 ; six_by_six[20] += s_1[2] * v0 ; six_by_six[21] += s_1[3] * v0 ;

            r_2[0] += s_1[0] * v1; r_2[1] += s_1[1] * v1; r_2[2] += s_1[2] * v1; r_2[3] += s_1[3] * v1;
            s_0[0]  = s_1[0] * v2; s_0[1]  = s_1[1] * v2; s_0[2]  = s_1[2] * v2; s_0[3]  = s_1[3] * v2;

            lowest_index = 2;
        }
        else {
            six_by_six[24] += r_2[0] * bb1; six_by_six[25] += r_2[1] * bb1; six_by_six[26] += r_2[2] * bb1; six_by_six[27] += r_2[3] * bb1;
            six_by_six[30] -= r_2[0] * bb2; six_by_six[31] -= r_2[1] * bb2; six_by_six[32] -= r_2[2] * bb2; six_by_six[33] -= r_2[3] * bb2;
            six_by_six[12] += r_2[0] * v0 ; six_by_six[13] += r_2[1] * v0 ; six_by_six[14] += r_2[2] * v0 ; six_by_six[15] += r_2[3] * v0 ;

            s_0[0] -= r_2[0] * v1; s_0[1] -= r_2[1] * v1; s_0[2] -= r_2[2] * v1; s_0[3] -= r_2[3] * v1;
            r_1[0]  = r_2[0] * v2; r_1[1]  = r_2[1] * v2; r_1[2]  = r_2[2] * v2; r_1[3]  = r_2[3] * v2;

            solution_vector[0] -= r_2[0] * v0 ; solution_vector[1] -= r_2[1] * v0 ; solution_vector[2] -= r_2[2] * v0 ; solution_vector[3] -= r_2[3] * v0 ;

            six_by_six[24] += s_2[0] * bb2; six_by_six[25] += s_2[1] * bb2; six_by_six[26] += s_2[2] * bb2; six_by_six[27] += s_2[3] * bb2;
            six_by_six[30] += s_2[0] * bb1; six_by_six[31] += s_2[1] * bb1; six_by_six[32] += s_2[2] * bb1; six_by_six[33] += s_2[3] * bb1;
            six_by_six[18] += s_2[0] * v0 ; six_by_six[19] += s_2[1] * v0 ; six_by_six[20] += s_2[2] * v0 ; six_by_six[21] += s_2[3] * v0 ;

            r_0[0] += s_2[0] * v1; r_0[1] += s_2[1] * v1; r_0[2] += s_2[2] * v1; r_0[3] += s_2[3] * v1;
            s_1[0]  = s_2[0] * v2; s_1[1]  = s_2[1] * v2; s_1[2]  = s_2[2] * v2; s_1[3]  = s_2[3] * v2;

            lowest_index = 0;
        }
    }

    if(lowest_index == 0) {
        six_by_six[0] += r_0[0] + sigma * r_1[0]; six_by_six[1] += r_0[1] + sigma * r_1[1]; six_by_six[2] += r_0[2] + sigma * r_1[2]; six_by_six[3] += r_0[3] + sigma * r_1[3];
        six_by_six[6] += s_0[0] + sigma * s_1[0]; six_by_six[7] += s_0[1] + sigma * s_1[1]; six_by_six[8] += s_0[2] + sigma * s_1[2]; six_by_six[9] += s_0[3] + sigma * s_1[3];
    }
    else if(lowest_index == 1) {
        six_by_six[0] += r_1[0] + sigma * r_2[0]; six_by_six[1] += r_1[1] + sigma * r_2[1]; six_by_six[2] += r_1[2] + sigma * r_2[2]; six_by_six[3] += r_1[3] + sigma * r_2[3];
        six_by_six[6] += s_1[0] + sigma * s_2[0]; six_by_six[7] += s_1[1] + sigma * s_2[1]; six_by_six[8] += s_1[2] + sigma * s_2[2]; six_by_six[9] += s_1[3] + sigma * s_2[3];
    }
    else {
        six_by_six[0] += r_2[0] + sigma * r_0[0]; six_by_six[1] += r_2[1] + sigma * r_0[1]; six_by_six[2] += r_2[2] + sigma * r_0[2]; six_by_six[3] += r_2[3] + sigma * r_0[3];
        six_by_six[6] += s_2[0] + sigma * s_0[0]; six_by_six[7] += s_2[1] + sigma * s_0[1]; six_by_six[8] += s_2[2] + sigma * s_0[2]; six_by_six[9] += s_2[3] + sigma * s_0[3];
    }

    // solve for alpha, beta
    double alpha = 0.0;
    double beta  = 0.0;
    partialSolver(six_by_six, solution_vector, l2, &alpha, &beta);

    complex_response.first  = alpha;
    complex_response.second =  beta;
    return complex_response;
}

void ReferencePoint::partialSolver(vector<double>& A, vector<double>& b, double l2, double *alpha, double *beta) {
    int length = static_cast<int>(b.size());

    // Forward elimination
    for(int cnt = 0; cnt < length-2; ++cnt) {
        int column = cnt * length;

        int max_index = cnt;
        int max_value = fabs( A[column + cnt] );
        for(int row = cnt + 1; row < length; ++row) {
            if(max_value < fabs( A[column + row] )) {
                max_value = fabs( A[column + row] );
                max_index = row;
            }
        }

        if(max_index != cnt) {
            // swap rows
            double tmp;
            for(int i = cnt; i < length; ++i) {
                double cc = i * length;

                tmp = A[cc + cnt];
                A[cc + cnt] = A[cc + max_index];
                A[cc + max_index] = tmp;
            }
            tmp = b[cnt];
            b[cnt] = b[max_index];
            b[max_index] = tmp;
        }

        // elimination
        for(int row = cnt + 1; row < length; ++row) {
            double multiplier = -A[column + row] / A[column + cnt];

            A[column + row] = 0.0;
            for(int col = cnt + 1; col < length; ++col) {
                A[col * length + row] += multiplier * A[col * length + cnt];
            }

            b[row] += multiplier * b[cnt];
        }
    }

    // Solve 2 x 2
    int index1 = length-2;
    int index2 = length-1;

    double scalar = l2 / (A[index1 * length + index1] * A[index2 * length + index2] -
                          A[index1 * length + index2] * A[index2 * length + index1]);
    *alpha = (  A[index2 * length + index2] * b[index1] - A[index2 * length + index1] * b[index2] ) * scalar;
    *beta  = ( -A[index1 * length + index2] * b[index1] + A[index1 * length + index1] * b[index2] ) * scalar;
}


void ReferencePoint::Print() {
  cout << "Reference node located at ( " << mPosition[0] << ", " << mPosition[1] << ", " << mPosition[2] << " )" << endl;
}
