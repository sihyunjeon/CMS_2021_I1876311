// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {
  /// @brief Normalised cross sections for 5.02 TeV for WW diboson production
  class CMS_2021_I1876311 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1876311);


    /// Book histograms and initialise projections before the run
    void init() {

      // final state of all stable particles
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState leptons(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      leptons.acceptTauDecays(false);
      DressedLeptons dressedleptons(photons, leptons, 0.1, Cuts::open(), true); 
      declare(dressedleptons, "DressedLeptons");

      // book histograms
      book(_h["511"], 5, 1, 1); // WW xsection vrs sqrt(s)
      book(_h["517"], 5, 1, 7); // WZ xsection vrs sqrt(s)
      book(_h["5113"], 5, 1, 13); // ZZ xsection vrs sqrt(s)

    }


    /// @brief Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& dressedleptons = apply<DressedLeptons>(event, "DressedLeptons").particlesByPt();
      if (dressedleptons.size() < 2 || dressedleptons.size() > 4) vetoEvent;

      // loop over two leptons
      for (const DressedLepton& l1 : dressedleptons){
          for (const DressedLepton& l2 : dressedleptons){
              // only pick up same flavor lepton pairs when they are not identical
              if (!isSame(l1, l2) && (l1.pid() + l2.pid() == 0)){
                  // if mll < 4GeV pair exists in event, escape
                  if ( (l1.momentum() + l2.momentum()).mass() < 4.0 ){
                      vetoEvent;
                  }
              }
          }
      }

      // WW
      if (dressedleptons.size() == 2){
	  if (dressedleptons[0].charge() + dressedleptons[1].charge() == 0){
              _h["511"]->fill(5.020);
	  }
      }

      // WZ
      if (dressedleptons.size() == 3){
          bool onshellZ = false;
	  if (abs(dressedleptons[0].pid() + dressedleptons[1].pid() + dressedleptons[2].pid()) == 11 || abs(dressedleptons[0].pid() + dressedleptons[1].pid() + dressedleptons[2].pid()) == 13){
              // test all 01, 02, 12 pairs
	      // should not be mutually exclusive (if, if, if), NOT using if/else if/ else
              if (dressedleptons[0].pid() + dressedleptons[1].pid() == 0){
                  double mll = (dressedleptons[0].momentum() + dressedleptons[1].momentum()).mass();
                  if (60 < mll && mll < 120) onshellZ = true;
              }
              if (dressedleptons[0].pid() + dressedleptons[2].pid() == 0){
                  double mll = (dressedleptons[0].momentum() + dressedleptons[2].momentum()).mass();
                  if (60 < mll && mll < 120) onshellZ = true;
              }
              if (dressedleptons[1].pid() + dressedleptons[2].pid() == 0){
                  double mll = (dressedleptons[1].momentum() + dressedleptons[2].momentum()).mass();
                  if (60 < mll && mll < 120) onshellZ = true;
              }
	  }
          if (onshellZ){
              _h["517"]->fill(5.020);
          }
      }

      // ZZ
      if (dressedleptons.size() == 4){
          bool onshellZ = false;
          if (dressedleptons[0].pid() + dressedleptons[1].pid() + dressedleptons[2].pid() + dressedleptons[3].pid() == 0){
              // 01 and 23 combination
              if (dressedleptons[0].pid() + dressedleptons[1].pid() == 0){
                  double mll1 = (dressedleptons[0].momentum() + dressedleptons[1].momentum()).mass();
                  double mll2 = (dressedleptons[2].momentum() + dressedleptons[3].momentum()).mass();
                  if (60 < mll1 && mll1 < 120){
                      if (60 < mll2 && mll2 < 120){
                          onshellZ = true;
                      }
                  }
              }
              // 02 and 13 combination
              if (dressedleptons[0].pid() + dressedleptons[2].pid() == 0){
                  double mll1 = (dressedleptons[0].momentum() + dressedleptons[2].momentum()).mass();
                  double mll2 = (dressedleptons[1].momentum() + dressedleptons[3].momentum()).mass();
                  if (60 < mll1 && mll1 < 120){
                      if (60 < mll2 && mll2 < 120){
                          onshellZ = true;
                      }
                  }
              }
              // 03 and 12 combination
              if (dressedleptons[0].pid() + dressedleptons[3].pid() == 0){
                  double mll1 = (dressedleptons[0].momentum() + dressedleptons[3].momentum()).mass();
                  double mll2 = (dressedleptons[1].momentum() + dressedleptons[2].momentum()).mass();
                  if (60 < mll1 && mll1 < 120){
                      if (60 < mll2 && mll2 < 120){
                          onshellZ = true;
                      }
                  }
              }

          }
          if (onshellZ){
              _h["5113"]->fill(5.020);
          }
      }
    }


    /// @brief Normalise histograms after the run
    void finalize() {
      double norm = (sumOfWeights() != 0) ? crossSection()/picobarn/sumOfWeights() : 1.0;

      scale(_h["511"], norm);///BrachingRatioWW);
      scale(_h["517"], norm);///BrachingRatioWZ);
      scale(_h["5113"], norm);///BrachingRatioZZ);
    }

  protected:
    double BrachingRatioWW;
    double BrachingRatioWZ;
    double BrachingRatioZZ;

  private:
    // Declaration of histograms
    map<string,Histo1DPtr> _h; 
    
  };

  RIVET_DECLARE_PLUGIN(CMS_2021_I1876311);
}

