#ifndef ____REFPNT__
#define ____REFPNT__

#include <iostream>
#include <vector>

#include "network_model/Node.h"
#include "network_model/HessianMatrix.h"

using namespace std;

class ReferencePoint{
public:
    ReferencePoint(double x, double y, double z, vector<Node> &, HessianMatrix &);
    ReferencePoint(const ReferencePoint& rp);
    ~ReferencePoint();

    vector< pair<double, double> > SingleModeFrequencyResponse(double omega);
    vector< pair<double, double> > DualModeFrequencyResponse  (double omega1, double omega2, double a);

    vector<double> GetPosition       () const { return mPosition       ; }

    double         GetPotentialL2    () const { return mPotentialL2    ; }
    double         GetPotentialSigma () const { return mPotentialSigma ; }
    vector<double> GetPotentialLambda() const { return mPotentialLambda; }

    double         GetTxxL2    () const { return mTxxL2    ; }
    double         GetTxxSigma () const { return mTxxSigma ; }
    vector<double> GetTxxLambda() const { return mTxxLambda; }

    double         GetTyyL2    () const { return mTyyL2    ; }
    double         GetTyySigma () const { return mTyySigma ; }
    vector<double> GetTyyLambda() const { return mTyyLambda; }

    double         GetTzzL2    () const { return mTzzL2    ; }
    double         GetTzzSigma () const { return mTzzSigma ; }
    vector<double> GetTzzLambda() const { return mTzzLambda; }

private:
    void   preCompute  ( vector<Node> &, HessianMatrix & );

    double innerProduct( vector<double>&, vector<double>& );
    void   grahamSchmit( vector<double>&, vector<double>&, double );
    double rescale     ( vector<double>& );
    void   partialSolver(vector<double>&, vector<double>&, double, double*, double*);

    void calculateFrequencyTerms(vector<double>& b_vector, double& bT_b, double& sigma, vector<double>& lambda, vector<Node> &nodes, HessianMatrix &hessianMatrix );

    pair<double, double> singleModeFrequencyResponse(double omega                          , double l2, double sigma, vector<double>& lambda);
    pair<double, double> dualModeFrequencyResponse  (double omega1, double omega2, double a, double l2, double sigma, vector<double>& lambda);

    vector<double> mPosition;

    // Potential
    double         mPotentialL2;
    double         mPotentialSigma;
    vector<double> mPotentialLambda;

    // Txx
    double         mTxxL2;
    double         mTxxSigma;
    vector<double> mTxxLambda;

    // Tyy
    double         mTyyL2;
    double         mTyySigma;
    vector<double> mTyyLambda;

    // Tzz
    double         mTzzL2;
    double         mTzzSigma;
    vector<double> mTzzLambda;
};

#endif
