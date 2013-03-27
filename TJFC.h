#ifndef TJFC_DEFINE_THINGY
#define TJFC_DEFINE_THINGY

// This is TJFC.h!  Enjoy :)

// TJFC.h / TJFC.cc
// Originally created August, 2009 by Josh Albert

///////////////////////////////////////////////////////////////////////////////
//
// class to calculate FC confidence intervals
// differs from TFeldmanCousins in 
// profile likelihoods to deal with systematic uncertainties in the background
//
///////////////////////////////////////////////////////////////////////////////
// Last updated Josh Albert September 25, 2009

#include "TObject.h"
#include "TNamed.h"
#include "TMath.h"
#include "TRandom.h"

#include <vector>
#include <map>
#include <iostream>
#include <cmath>

class TJFC : public TNamed
{
private:
  Double_t fCL;              // Confidence level as a fraction (0.9 for 90%)
  Double_t fUpperLimit;      // Calculated Upper Limit
  Double_t fLowerLimit;      // Calculated Lower limit
  Int_t fMaxBins;            // Maximum number of bins we're willing to look at
  Int_t fErrorCode;          // What's wrong this time?
  bool fProfile;             // Do we want to use profile likelihood?
  bool fSmear;               // Do we want to use the uncertainty smear?
  bool fBkgSmear;            // Do we want to smear the background uncertianty?
  bool fSigSmear;            // Do we want to smear the signal uncertainty?
  Int_t fTopBin;             // The highest bin we expect to need for mu & b
  Double_t fPrecision;       // How precisely must we know the intervals?
  Int_t fMaxIterations;      // How many steps to take to find interval?
  bool fVerbose;             // Do we output extra text?
  Double_t fBkgSysSmear;     // Background Systematic error (only for Smear)
  Double_t fSigSysSmear;     // Signal Systematic error (only for Smear)
  bool fFullMinimize;        // Do we minimize the ratio denomenator fully, or
                             // are we restricted to b_best from the numerator?
  Int_t fMaxCheckIter;    // Number of steps we take above initial found UL
  bool fDoubleSmear;         // Smear the likelihood ratio?

  // Parameters used in the distribution calculations
  Double_t f_mu;             // Expected signal
  Double_t f_b;              // Expected background

  // Upper and lower band limits
  Int_t fBandUp;
  Int_t fBandDn;

  // Parameters used for smearing and the separate absolute minimum
  bool fSeparateSmear;       // Will we be doing the full sig & bkg smearing?
  bool fUseAbsMin;           // Do we have an absolute minimum expected evt 
  bool fNAbsMin;             // Absolute minimum value for N_exp (for N_best)

  // It is well worth saying here and everywhere that limits are inclusive
  void CalcTopBin();         // Set the highest bin we will need for mu & b
  Double_t GaussianProb(Double_t, Double_t, Double_t);
  Double_t MB(Double_t, Double_t, Int_t);  // Calculate mu_best
  void CalculateBandNoUncertainty(Double_t, Double_t);
  Double_t CalculateLowerLimit(Double_t, Double_t);
  Double_t CalculateUpperLimit(Double_t, Double_t);
  void CalcTopBinGenerous(); // Set the highest bin...a bit higher
  Double_t CheckUpperLimit(Double_t, Double_t, Double_t); // check upwards
  Double_t PoissonProbSepSmear(Int_t, Double_t, Double_t);

public:
  /* Constructor */
  TJFC(Double_t CL = 0.9, Option_t *option = "");

  /* Destructor */
  virtual ~TJFC();

  /* Get and set the Confidence Level */
  Double_t GetCL() const     {  // Get the CL
    return fCL;
  }
  
  Double_t PoissonProb(Int_t i , Double_t sigmean, Double_t bkgmean); // Get Appropriate Poisson Probability
  Double_t PoissonProbSmear(Int_t i, Double_t sigmean, Double_t bkgmean); // Get Smeared Poisson Probability

  void SetCL(Double_t CL)    { // Change the CL
    if (CL <= 0 || CL >= 1) {
      std::cout << "What are you doing!?\n";
      std::cout << "Setting CL to 90%, you moron!\n";
      fCL = 0.9;
    }
    fCL = CL;
  }
  Double_t GetUpperLimit()   {  // Return last calculated UL
    return fUpperLimit;
  }
  Double_t GetLowerLimit()   {  // Return last calculated LL
    return fLowerLimit;
  }
  Int_t GetBandUp()          {  // Return last calculated band top
    return fBandUp;
  }
  Int_t GetBandDn()          {  // Return last calculated band bottom
    return fBandDn;
  }
  Double_t GetPrecision()    {  // Return the precision we are using
    return fPrecision;
  }
  void SetPrecision(Double_t prec)        {  // Set Precision
    if(prec>0.0000001 && prec<0.5) fPrecision = prec;
    else std::cout << "Invalid Precision!\n";
  }

  Int_t GetTopBin() const    {  // Top bin to be used in calculations
    return fTopBin;
  } 

  void CalculateBand(Double_t mu, Double_t b)  {
    if (fCL <= 0 || fCL >= 1) {
      std::cout << "What are you doing!?\n";
      std::cout << "Setting CL to 90%, you moron!\n";
      fCL = 0.9;
    }
    CalculateBandNoUncertainty(mu, b);
  }

  Int_t CalculateBandDn(Double_t mu, Double_t b)  {
    CalculateBand(mu, b);
    return fBandDn;
  }
  Int_t CalculateBandUp(Double_t mu, Double_t b)  {
    CalculateBand(mu, b);
    return fBandUp;
  }

  Double_t GetLowerLimit(Double_t n, Double_t b, Double_t CL){
    if (CL <= 0 || CL >= 1) {
      std::cout << "What are you doing!?\n";
      std::cout << "Setting CL to 90%, you moron!\n";
      fCL = 0.9;
    }
    else fCL = CL;
    return CalculateLowerLimit(n, b);
  }

  Double_t GetLowerLimit(Double_t n, Double_t b){
    return CalculateLowerLimit(n, b);
  }

  Double_t GetUpperLimit(Double_t n, Double_t b, Double_t CL){
    if (CL <= 0 || CL >= 1) {
      std::cout << "What are you doing!?\n";
      std::cout << "Setting CL to 90%, you moron!\n";
      fCL = 0.9;
    }
    else fCL = CL;
    return CalculateUpperLimit(n, b);
  }

  Double_t GetUpperLimit(Double_t n, Double_t b){
    return CalculateUpperLimit(n, b);
  }

  void SetSmearOn(Double_t sys) {
    if (sys <= 0.0) {
      std::cout << "Invalid Background Systematic, You Suck at This Game!\n";
      return;
    }
    fBkgSysSmear = sys;
    fSmear = true;
  }

  void SetDoubleSmear(bool doublesmear) {
    // Do we want to smear both the prob distribution and the likelihood ratio?
    fDoubleSmear = doublesmear;
  }

  bool GetDoubleSmear() {
    // Are we smearing the likelihood ratio like we're smearing the prob dist?
    return fDoubleSmear;
  }

  void SetSmearOff() {
    fSmear = false;
  }

  void SetBkgSysSmear(Double_t sys) {
    if (sys <= 0.0) {
      std::cout << "Invalid Background Systematic, You Suck at This Game!\n";
      return;
    }
    fBkgSysSmear = sys;
  }

  Double_t GetBkgSysSmear() {
    return fBkgSysSmear;
  }

  void SetSigSysSmear(Double_t sys) {
    if (sys <= 0.0) {
      std::cout << "Invalid Signal Systematic, You Suck at This Game!\n";
      return;
    }
    fSigSysSmear = sys;
  }

  Double_t GetSigSysSmear() {
    return fSigSysSmear;
  }

  void SetAbsoluteMin( Double_t absmin) {
    if (absmin < 0.0) {
      std::cout << "Invalid absolute minimum, stop sucking so much!\n";
      return;
    }
    fNAbsMin = absmin;
    fUseAbsMin = true;
  }

  void SetAbsoluteMinOff() {
    fUseAbsMin = false;
  }

  Double_t GetAbsoluteMin() {
    if (!fUseAbsMin) return -1.0;
    return fNAbsMin;
  }

  void SetSeparateSmearOn() {
    fSeparateSmear = true;
  }

  void SetSeparateSmearOff() {
    fSeparateSmear = false;
  }

  bool GetSeparateSmear() {
    return fSeparateSmear;
  }

  void SetVerbose(bool verb) {
    fVerbose = verb;
  }

  Int_t GetMaxCheckIter() {
    return fMaxCheckIter;
  }

  void SetMaxCheckIter( Int_t maxit ) { // Number of steps to take upwards
      // to make sure the interval finder didn't get stuck on a fluctuation.
      // Setting higher slows things but increases resolution.
      // Honestly, you probably never need to touch this.
    fMaxCheckIter = maxit;
  }

  ClassDef(TJFC, 1)
};

// Improved Calculation of FC confidence intervals
//
#endif
