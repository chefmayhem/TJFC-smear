#include "TJFC.h"
#include <iostream>
using namespace std;

ClassImp(TJFC)

//**********  CONSTRUCTOR  **********************
// Make the gorram class constructor and destructor!
TJFC::TJFC(Double_t CL, Option_t * /*option*/)
:  fCL(CL),
   fUpperLimit(0.0),
   fLowerLimit(0.0)
{
  // Upon construction, initialize things
  f_mu = 0;
  f_b = 0;
  fTopBin = 20;
  fPrecision = 0.0005;
  fMaxIterations = 100;
  fVerbose = false;
  fSmear = false;
  fBkgSysSmear = 0.0;
  fSigSysSmear = 0.0;
  fSeparateSmear = false;
  fMaxCheckIter = 50;
  fDoubleSmear = false;
  fUseAbsMin = false;
  fNAbsMin = 0.0;
}

//**********  DESTRUCTOR  ***********************
TJFC::~TJFC()
{
}

//***********************************************
// What is the highest bin we should look at for given s & b
void TJFC::CalcTopBin() {
  // Function is now modified to be a function of CL
  // Obviously, for higher CL, we need to check more.
  fTopBin = (int)(2 * floor(f_mu + f_b) + 4);
  if (fCL > 0.9) {
    fTopBin += (int)(floor((f_mu + f_b + 1) * 0.2 / (1.00-fCL)) - 2 * floor(f_mu + f_b));
  }
  if (fTopBin > 100) {
    cout << "Warning: Checking over 100 bins!  You sure this is right?\n";
    cout << "Top bin is: " << fTopBin << "\n";
  }
}

//***********************************************
void TJFC::CalcTopBinGenerous() {
  // For profile use when the normalization requirement
  // demands knowing all non-zero probabilities.
  fTopBin = (int)(3 * floor(f_mu + 2 * f_b) + 8);
}

// Some calculation internals

//*********  FIND MU BEST  **********************
Double_t TJFC::MB(Double_t mu, Double_t b, Int_t n) {
  // If we have the abs min, we need to be ready for the 
  // case of the abs min being less than b, yielding a 
  // negative mu_best.  I think this is ok, but not positive.
  if (fUseAbsMin) {
    if (n > fNAbsMin) { // business as usual
      return n-b;
    }
    else if (n <= fNAbsMin) { // the absolute minimum comes into play!
      return fNAbsMin - b; // so that S + B = N_absmin
    }
  }
  // And then the normal case...
  if (n > b) return n-b; // if best signal is positive
  return 0; // because signal must be non-negative
}

//********* POISSON PROB MC *********************
Double_t TJFC::PoissonProbSmear(Int_t i, Double_t sigmean, Double_t bkgmean) {
  if (bkgmean <= 0) { cout << "This could be bad!!!!\n";}
  Int_t interpnmin = 50; Int_t interpnmax = 5000;
  Int_t interpn = TMath::Max(interpnmin, //...
    (Int_t)(interpnmax/(Double_t)((i+1)*(i+1)*(i+1))));
  Double_t gaussmean = bkgmean;
  Double_t gausswidth = fBkgSysSmear * bkgmean;
  Double_t interpmin = TMath::Max(0.0, gaussmean - 4.0 * gausswidth);
  Double_t interpmax = gaussmean + 4.0 * gausswidth;
  Double_t interpstep = (interpmax - interpmin)/((Double_t)interpn);
//  cout << interpmin << " " << interpmax << " " << interpstep << "\n";
  Double_t gausscorr = 0.5 * (1 + TMath::Erf(gaussmean/sqrt(2)/gausswidth));
  Double_t probtot = 0;
  for (int i_interp = 0; i_interp != interpn; i_interp++) {
    Double_t interpmean = interpmin + (i_interp + 0.5) * interpstep;
    probtot += TMath::Poisson(i, interpmean + sigmean) //...
      * TMath::Gaus(interpmean, gaussmean, gausswidth, kTRUE) //...
      * interpstep / gausscorr;
  } // end interpolation over gaussian for loop
  return probtot;
}  // End TJFC::PoissonProbSmear

//********* POISSON PROB MC *********************
Double_t TJFC::PoissonProbSepSmear(Int_t i, Double_t sigmean, Double_t bkgmean) {
  // Combine sig & bkg sys errors to make smear in both sig & bkg
  if (bkgmean + sigmean <= 0) { 
    cout << "This could be bad!!!! " << bkgmean << " " << sigmean << "\n";
    }
  if (fDoubleSmear) {cout << "Not compatible with double smear!\n";}
  Int_t interpnmin = 50; Int_t interpnmax = 5000;
  Int_t interpn = TMath::Max(interpnmin, //...
    (Int_t)(interpnmax/(Double_t)((i+1)*(i+1)*(i+1))));
  Double_t gaussmean = bkgmean+sigmean;
  Double_t gausswidth = pow(fBkgSysSmear * bkgmean, 2); // continues...
  gausswidth += pow(fSigSysSmear * sigmean, 2); gausswidth = sqrt(gausswidth);
  Double_t interpmin = TMath::Max(0.0, gaussmean - 4.0 * gausswidth);
  Double_t interpmax = gaussmean + 4.0 * gausswidth;
  Double_t interpstep = (interpmax - interpmin)/((Double_t)interpn);
//  cout << interpmin << " " << interpmax << " " << interpstep << "\n";
  Double_t gausscorr = 0.5 * (1 + TMath::Erf(gaussmean/sqrt(2)/gausswidth));
  Double_t probtot = 0;
  for (int i_interp = 0; i_interp != interpn; i_interp++) {
    Double_t interpmean = interpmin + (i_interp + 0.5) * interpstep;
    probtot += TMath::Poisson(i, interpmean) //...
      * TMath::Gaus(interpmean, gaussmean, gausswidth, kTRUE) //...
      * interpstep / gausscorr;
  } // end interpolation over gaussian for loop
  return probtot;
}  // End TJFC::PoissonProbSmear

//********* POISSON PROB ************************
Double_t TJFC::PoissonProb(Int_t i, Double_t sigmean, Double_t bkgmean) {
  if (fSeparateSmear == true) return PoissonProbSepSmear(i, sigmean, bkgmean);
  else if (fSmear == true) return PoissonProbSmear(i, sigmean, bkgmean);
  return TMath::Poisson(i, sigmean + bkgmean);
}

//*********  FIND N IN BAND GIVEN MU,B  *********
void TJFC::CalculateBandNoUncertainty(Double_t mu, Double_t b) {
  // FOR TESTING BY JOSH 06.03.10
  if (mu < 0.0) cout << "WTF MU= " << mu << "\n";
  // This finds the "horizontal lines" for confidence belt construction
  Double_t mub; // mu_best
  Double_t ifill;
  std::map<Int_t, Double_t> ratio, prob;
  f_mu = mu; f_b = b;
//  CalcTopBin(); // sets highest bin we will check given mu, b
  CalcTopBinGenerous(); // sets highest bin we will check given mu, b
  Int_t minbin, maxbin;
  ifill = 0;
  ratio.clear(); prob.clear(); // initialize
  for (Int_t i = 0; i != fTopBin+1; ++i) { // calculate probs
    prob[i] = PoissonProb(i, mu, b);
    if (fDoubleSmear) // then
      ratio[i] = PoissonProb(i, mu, b)/PoissonProb(i, MB(mu,b,i), b);
    else // then
      ratio[i] = TMath::Poisson(i, mu+b)/TMath::Poisson(i, MB(mu,b,i)+b);
  }
  // Perform the sorting!!
  if (fVerbose) std::cout << "*** Performing sorting!! ***\n";
  ifill = 0; minbin = ratio.size(); maxbin = 0;
  for (Int_t i = 0; i != ratio.size(); ++i) { // sorting loop
    Int_t tmpmaxpos = -1; // bin of highest ratio
    Double_t tmpmaxval = 0; // value of highest ratio
    for (Int_t j = 0; j != (int)(ratio.size()-1); ++j) {
      if (ratio[j] > tmpmaxval) {
        tmpmaxval = ratio[j];
        tmpmaxpos = j;
      } // end if checking the bin
    } // end for looping over bins
    if (tmpmaxpos == -1) {
      std::cout << "FLAGRANT SYSTEM ERROR NO SYS!\n";
      std::cout << "ifill = " << ifill << "\n";
      std::cout << "mu b i " << mu << " " << b << " " << i << "\n";
      std::cout << "ratio.size fTopBin " << ratio.size() << " " << fTopBin << "\n";
      break;
    }
    ifill += prob[tmpmaxpos]; // increment confidence interval total
    if (minbin > tmpmaxpos) minbin = tmpmaxpos;
    if (maxbin < tmpmaxpos) maxbin = tmpmaxpos;
    // added by jba 070210, let's make sure CalcTopBin is good enough
    if (tmpmaxpos >= fTopBin) { // This really should never happen
      cout << "This is bad!  fTopBin = " << fTopBin << //...
        " and we just selected bin " << tmpmaxpos << "\n";
    } // end the jba 070210 testing section
    ratio[tmpmaxpos] = 0; // so we don't pick it again
    if (ifill > fCL) break; // we have reached the CL!
  } // end sorting loop
  fBandUp = maxbin;
  fBandDn = minbin;
}

//**********  CALCULATE LOWER LIMIT  ************
Double_t TJFC::CalculateLowerLimit(Double_t n, Double_t b) {
  Double_t mutest;    // our mu being tested to find the best value
  Double_t mulast;    // previous value of mutest
  Double_t mustep;    // initial mu step value
  bool initialrun;    // We have yet to cross a boundry, do not shrink mustep
  Int_t lastdir;       // Was our last band too high, low, or non-existent
  Double_t stepratio; // 0.5 for binary search, step ratio for each interation


  //*********** INITIALIZATION *******************************
  // We will use a binary search for this.
  // Choose starting guess for mu
  // The starting guess will also be our starting increment

  mutest = (n-b)/2.0;
  mulast = -999.;
  mustep = mutest;
  stepratio = 0.50;


  if (mutest <= 0) {   // Expect zero to be in our interval. Treat differently.
    mutest = 0.0;
    mustep = 1.0;
  }  // end if for mu ~ 0
  initialrun = true;
  if (fVerbose) cout << "initial mutest = " << mutest << "\n";
  lastdir = 0;  // We have yet to reverse binary search direction
  
  //*************** MAIN LOOP *******************************
  for (Int_t i = 0; i != fMaxIterations; ++i) {
    CalculateBand(mutest, b);

    // Too high, go lower //////////////////
    if (n <= fBandUp) {  // n is in the band
      if (mutest == 0) { // because mu>0... after corrections...
        if (fVerbose) std::cout << "Found it!\n";
        fLowerLimit = 0;
        return fLowerLimit;
      }
      if (!initialrun) mustep = mustep*stepratio; // decrease step
      if (lastdir * -1 == -1) initialrun = false; // a change of direction!
      lastdir = -1;
    }

    // Too low, go higher //////////////////
    if (n > fBandUp) { // n is above the band
      if (!initialrun) mustep = mustep*stepratio; // decrease step
      if (lastdir * 1 == -1) initialrun = false; // a change of direction!
      lastdir = 1;
    }

    // Change things ///////////////////////
    mulast = mutest;
    mutest = mutest + mustep * lastdir;
  
    // Check for stopping point ////////////
    if (sqrt((mutest-mulast)*(mutest-mulast)) <= fPrecision) {
      if (fVerbose) std::cout << "Found it!\n";
      fLowerLimit = mutest;
      return fLowerLimit;
    } // end stopping point check

  }
  std::cout << "Failed to converge\n";
  return -1;
  
}

//**********  CALCULATE UPPER LIMIT  ************
Double_t TJFC::CalculateUpperLimit(Double_t n, Double_t b) {
  Double_t mutest;    // our mu being tested to find the best value
  Double_t mulast;    // previous value of mutest
  Double_t mustep;    // initial mu step value
  bool initialrun;    // We have yet to cross a boundry, do not shrink mustep
  Int_t lastdir;       // Was our last band too high, low, or non-existent
  Double_t stepratio; // 0.5 for binary search, step ratio for each interation

  //*********** INITIALIZATION *******************************
  // We will use a binary search for this.
  // Choose starting guess for mu
  // The starting guess will also be our starting increment
  mutest = n*2.0;
  mulast = -999.;
  mustep = mutest/2.0;
  stepratio = 0.50;

  if (mutest <= 0) {   // We can't allow a zero upper limit...
  mutest = 3.0;
  mustep = 1.0;
  }  // end if for mu ~ 0
  initialrun = true;
  if (fVerbose) cout << "initial mutest = " << mutest << "\n";
  lastdir = 0;  // We are on our first step
  
  //*************** MAIN LOOP *******************************
  for (Int_t i = 0; i != fMaxIterations; ++i) {
    CalculateBand(mutest, b);
    cout << "CalcingBand: mutest, b, fBandDn: " << mutest << "," << b << "," << fBandDn << "\n";

    // Too high, go lower //////////////////
    if (n < fBandDn) {  // n is below the band
      if (!initialrun) mustep = mustep*stepratio; // decrease step
      if (lastdir * -1 == -1) initialrun = false; // a change of direction!
      lastdir = -1;
    }

    // Too low, go higher //////////////////
    if (n >= fBandDn) { // n is in the band
      if (!initialrun) mustep = mustep*stepratio; // decrease step
      if (lastdir * 1 == -1) initialrun = false; // a change of direction!
      lastdir = 1;
    }

    // Change things ///////////////////////
    mulast = mutest;
    mutest = mutest + mustep * lastdir;
  
    // Check for stopping point ////////////
    if (fabs(mutest-mulast) <= fPrecision) {
      Double_t checktest = CheckUpperLimit(mutest, b, n);
      if (checktest < 0.) { // doesn't look like a ripple
        if (fVerbose) std::cout << "Found it!\n";
        fUpperLimit = mutest;
        return fUpperLimit;
      }
      else { // it's a trap!  Initial UL was just a fluctuation
        mutest = checktest;   // Check higher value
        mustep = mustep * 8;  // Back to larger step, begin searching anew
      }
    } // end stopping point check

  }
  std::cout << "Failed to converge\n";
  return -1;
}

//**********  CHECK UPPER LIMIT  ****************
Double_t TJFC::CheckUpperLimit(Double_t mu_i, Double_t b, Double_t n) {
  bool ripple = false;  // innocent until proven guilty
  Double_t newUL = mu_i;
  // Let the step we take to check be the larger of:
  // -  the precision times 2
  // -  0.005
  Double_t stepval = 0.005;
  if (n <= b) { // otherwise it's a waste of time!
    if (fPrecision * 2 > 0.005) stepval = fPrecision*2;
    for (Int_t step = 1; step != fMaxCheckIter; ++step) {
      Double_t stepBandDn = //...
        CalculateBandDn(mu_i + step * stepval, b);
      if (n >= stepBandDn) { // we were just caught on a ripple!
        ripple = true; newUL = mu_i + step * stepval;
      }
    }
  }
  if (ripple) {
    if (fVerbose) std::cout << "rippled!! old mu was : " << mu_i << "\n";
    return newUL; }
  else return -1;
}
