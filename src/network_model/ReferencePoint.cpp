#include <cmath>

#include "network_model/ReferencePoint.h"

ReferencePoint::ReferencePoint(double x, double y, double z) {
    mPosition.push_back(x);
    mPosition.push_back(y);
    mPosition.push_back(z);

    mBTb   = 0.0;
    mSigma = 0.0;
}
ReferencePoint::~ReferencePoint()
{
}


void ReferencePoint::BuildPotential(vector<Node> &nodes, HessianMatrix &hessianMatrix) {
    vector<double> electricPotential;
    electricPotential.reserve( 3 * nodes.size() );

    for(auto node : nodes) {
        double q = node.GetQ();

        double dx = mPosition[0] - node.GetX();
        double dy = mPosition[1] - node.GetY();
        double dz = mPosition[2] - node.GetZ();

        double constant_value = q * pow(dx*dx + dy*dy + dz*dz, -1.5);
        electricPotential.push_back( dx * constant_value );
        electricPotential.push_back( dy * constant_value );
        electricPotential.push_back( dz * constant_value );
    }

    mBTb = innerProduct(electricPotential, electricPotential);
//cout << "Potential magnitude = " << mBTb << endl;

    vector<double> result;
    result.resize( 3 * nodes.size() );

    vector<double> result2;
    result2.resize( 3 * nodes.size() );

    hessianMatrix.MultiplyMatrix(electricPotential, result);
    mLambda.push_back( innerProduct(electricPotential, result) / mBTb );

    double last_lambda[3] = {0.0, 0.0, 0.0};
    double scalar_value;

    int cnt      = 0;
    int term_cnt = 0;
    while( term_cnt < 2 ) {
        last_lambda[2] = *(mLambda.rbegin());

        // Two steps per loop to address ping-pong of vectors
        grahamSchmit( result, electricPotential, *(mLambda.rbegin()) );
        scalar_value = rescale( result );
        hessianMatrix.MultiplyMatrix(result, result2);
        mLambda.push_back( scalar_value * innerProduct(electricPotential, result2) / mBTb );
        last_lambda[1] = *(mLambda.rbegin());

        grahamSchmit( result2, electricPotential, *(mLambda.rbegin()) );
        scalar_value = rescale( result2 );
        hessianMatrix.MultiplyMatrix(result2, result);
        mLambda.push_back( scalar_value * innerProduct(electricPotential, result) / mBTb );
        last_lambda[0] = *(mLambda.rbegin());

//cout << cnt+1 << ") " << last_lambda[1] << " ratio = " << last_lambda[1] / last_lambda[2] << endl;
//cout << cnt+2 << ") " << last_lambda[0] << " ratio = " << last_lambda[0] / last_lambda[1] << endl;
//cout << "Convergence criteria = " << fabs( (last_lambda[1]*last_lambda[1])/(last_lambda[0]*last_lambda[2]) - 1.0) << endl;
        cnt += 2;
        if(cnt < 10)
            continue;

        if(fabs( (last_lambda[1]*last_lambda[1])/(last_lambda[0]*last_lambda[2]) - 1.0) < 1e-6) {
            term_cnt += 1;
        }
        else{
            term_cnt = 0;
        }


    }

    mSigma = last_lambda[0] / last_lambda[1];

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


pair<double, double> ReferencePoint::SingleModeFrequencyResponse(double omega) {
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
    four_by_four[ 2] =  mSigma;             four_by_four[ 3] = omega; four_by_four[ 6] = -four_by_four[ 3]; four_by_four[ 7] = four_by_four[ 2];
    four_by_four[ 8] =  mLambda[0];         four_by_four[ 9] = omega; four_by_four[12] = -four_by_four[ 9]; four_by_four[13] = four_by_four[ 8];
    four_by_four[10] = *(mLambda.rbegin());                                                                 four_by_four[15] = four_by_four[10];

    solution_vector[0] = 1.0;

    // iterate
    double tmp_value;
    for(int cnt = 0; cnt < static_cast<int>(mLambda.size() - 1); ++cnt) {
        tmp_value = mLambda[cnt+1] / omega;
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
    partialSolver(four_by_four, solution_vector, &alpha, &beta);

    complex_response.first  = alpha;
    complex_response.second =  beta;
    return complex_response;
}
pair<double, double> ReferencePoint::DualModeFrequencyResponse(double omega1, double omega2, double mix_ratio) {
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
    six_by_six[ 4] = mSigma*mSigma - omega_prod;                  six_by_six[ 5] = omega_plus*mSigma;              six_by_six[10] = -six_by_six[ 5]; six_by_six[11] = six_by_six[ 4];
    six_by_six[12] = mLambda[0];                                  six_by_six[13] = omega_plus;                     six_by_six[18] = -six_by_six[13]; six_by_six[19] = six_by_six[12];
    six_by_six[14] = mLambda[1] - omega_prod;                                                                                                        six_by_six[21] = six_by_six[14];
    six_by_six[16] = *(mLambda.rbegin());                                                                                                            six_by_six[23] = six_by_six[16];
    six_by_six[24] = mLambda[0]*mLambda[0]+mLambda[1]-omega_prod; six_by_six[25] = omega_plus*mLambda[0];          six_by_six[30] = -six_by_six[25]; six_by_six[31] = six_by_six[24];
    six_by_six[26] = mLambda[0]*mLambda[1]+mLambda[2];            six_by_six[27] = omega_plus*mLambda[1];          six_by_six[32] = -six_by_six[27]; six_by_six[33] = six_by_six[26];
    six_by_six[28] = *(mLambda.rbegin())*(mLambda[0]+mSigma);     six_by_six[29] = *(mLambda.rbegin())*omega_plus; six_by_six[34] = -six_by_six[29]; six_by_six[35] = six_by_six[28];

    solution_vector[0] = mLambda[0];
    solution_vector[1] = omega_bar;
    solution_vector[2] = mLambda[1];
    solution_vector[4] = *(mLambda.rbegin());

    r_1[0] = 1.0; r_1[3] =  omega_plus;
    s_1[1] = 1.0; s_1[2] = -omega_plus;
    r_2[2] = 1.0;
    s_2[3] = 1.0;

    lowest_index = 1;

    // iterate
    double bb1, bb2, v0, v1, v2;
    v2  = 1.0 / omega_prod;
    v1  = omega_plus / omega_prod;
    for(int cnt = 1; cnt < static_cast<int>(mLambda.size() - 1); ++cnt) {
        if(cnt == static_cast<int>(mLambda.size() - 2)) {
            bb1 = (mLambda[0] + mSigma) * mLambda[cnt+1] * v2;
        }
        else {
            bb1 = (mLambda[0]*mLambda[cnt+1] + mLambda[cnt+2]) * v2;
        }
        bb2 = mLambda[cnt+1] * v1;
        v0  = mLambda[cnt+1] * v2;

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
        six_by_six[0] += r_0[0] + mSigma * r_1[0]; six_by_six[1] += r_0[1] + mSigma * r_1[1]; six_by_six[2] += r_0[2] + mSigma * r_1[2]; six_by_six[3] += r_0[3] + mSigma * r_1[3];
        six_by_six[6] += s_0[0] + mSigma * s_1[0]; six_by_six[7] += s_0[1] + mSigma * s_1[1]; six_by_six[8] += s_0[2] + mSigma * s_1[2]; six_by_six[9] += s_0[3] + mSigma * s_1[3];
    }
    else if(lowest_index == 1) {
        six_by_six[0] += r_1[0] + mSigma * r_2[0]; six_by_six[1] += r_1[1] + mSigma * r_2[1]; six_by_six[2] += r_1[2] + mSigma * r_2[2]; six_by_six[3] += r_1[3] + mSigma * r_2[3];
        six_by_six[6] += s_1[0] + mSigma * s_2[0]; six_by_six[7] += s_1[1] + mSigma * s_2[1]; six_by_six[8] += s_1[2] + mSigma * s_2[2]; six_by_six[9] += s_1[3] + mSigma * s_2[3];
    }
    else {
        six_by_six[0] += r_2[0] + mSigma * r_0[0]; six_by_six[1] += r_2[1] + mSigma * r_0[1]; six_by_six[2] += r_2[2] + mSigma * r_0[2]; six_by_six[3] += r_2[3] + mSigma * r_0[3];
        six_by_six[6] += s_2[0] + mSigma * s_0[0]; six_by_six[7] += s_2[1] + mSigma * s_0[1]; six_by_six[8] += s_2[2] + mSigma * s_0[2]; six_by_six[9] += s_2[3] + mSigma * s_0[3];
    }

    // solve for alpha, beta
    double alpha = 0.0;
    double beta  = 0.0;
    partialSolver(six_by_six, solution_vector, &alpha, &beta);

    complex_response.first  = alpha;
    complex_response.second =  beta;
    return complex_response;
}

void ReferencePoint::partialSolver(vector<double>& A, vector<double>& b, double *alpha, double *beta) {
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
                double col = i * length;

                tmp = A[col + cnt];
                A[col + cnt] = A[col + max_index];
                A[col + max_index] = tmp;
            }
            tmp = b[cnt];
            b[cnt] = b[max_index];
            b[max_index] = tmp;
        }

        // elimination
        for(int row = cnt + 1; row < length; ++row) {
            double multiplier = -A[column + row] / A[column + cnt];

            A[column + row] = 0.0;
            for(int col = cnt + 1; cnt < length; ++col) {
                A[col * length + row] += multiplier * A[col * length + cnt];
            }

            b[row] += multiplier * b[cnt];
        }
    }

    // Solve 2 x 2
    int index1 = length-2;
    int index2 = length-1;
    double scalar = mBTb / (A[index1 * length + index1] * A[index2 * length + index2] - A[index1 * length + index2] * A[index2 * length + index1]);
    *alpha = (  A[index2 * length + index2] * b[index1] - A[index2 * length + index1] * b[index2] ) * scalar;
    *beta  = ( -A[index1 * length + index2] * b[index1] + A[index1 * length + index1] * b[index2] ) * scalar;
}


void ReferencePoint::Print() {
    cout << "Reference position: (" << mPosition[0]
         << ", " << mPosition[1]
         << ", " << mPosition[2] << ")" << endl;
    cout << "Potential magnitude (bTb) = " << mBTb << endl;
    for(int cnt = 0; cnt < static_cast<int>(mLambda.size()); ++cnt) {
        cout << "Lambda[" << cnt << "] = " << mLambda[cnt] << endl;
    }
    cout << "Sigma = " << mSigma << endl;
}
