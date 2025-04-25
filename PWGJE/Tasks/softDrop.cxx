//----------------------------------------------------------------------
/// \file example_softdrop_o2.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::SoftDrop class is used within the O2Physics framework.
//----------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "Framework/runDataProcessing.h" // O2Physics header
#include "Framework/AnalysisTask.h"      // O2Physics header
#include "Framework/AnalysisDataModel.h" // O2Physics header
#include "Framework/HistogramRegistry.h"  // O2Physics header


using namespace o2;
using namespace o2::framework;
using namespace o2::constants::math;
using namespace fastjet; 
using namespace std;

// Forward declaration to make things clearer
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
struct SoftDropTask {
  // Define the input and output data
  HistogramRegistry registry{"registry", {}};

  // Define the parameters for SoftDrop
  Configurable<double> z_cut{"z_cut", 0.10, "Symmetry cut parameter for SoftDrop"};
  Configurable<double> beta{"beta", 2.0, "Beta parameter for SoftDrop"};
  Configurable<double> R{"R", 1.0, "Jet radius parameter"};
  Configurable<double> ptmin{"ptmin", 20.0, "Minimum jet pT"};

  void init(InitContext&)
  {
    // Initialize the SoftDrop groomer
    sd = new contrib::SoftDrop(beta, z_cut);
    cout << "SoftDrop groomer is: " << sd->description() << endl;
    // Initialize histograms
    registry.add("NEvents", "Number of Events", HistType::kTH1F, {{1, 0, 1}});
    registry.add("jet_pt", "Original Jet pT;pT (GeV/c);Counts", HistType::kTH1F, {{100, 0, 200}});
    registry.add("sd_jet_pt", "SoftDropped Jet pT;pT (GeV/c);Counts", HistType::kTH1F, {{100, 0, 200}});
    registry.add("jet_mass", "Original Jet Mass;Mass (GeV/c^2);Counts", HistType::kTH1F, {{100, 0, 200}});
    registry.add("sd_jet_mass", "SoftDropped Jet Mass;Mass (GeV/c^2);Counts", HistType::kTH1F, {{100, 0, 200}});
    registry.add("DeltaR", "Delta R between subjets;Delta R;Counts", HistType::kTH1F, {{100, 0, 1}});
    registry.add("hSymmetry", "Symmetry measure (z);z;Counts", HistType::kTH1F, {{100, 0, 1}});
    registry.add("hMu", "Mass drop (mu);mu;Counts", HistType::kTH1F, {{100, 0, 1}});
    registry.add("hJetPtDeltaR", "Jet pT vs rg ;pT (GeV/c);rg", HistType::kTH2F, {{100, 0, 200}, {100, 0, 1}});
    registry.add("rg", "rg", HistType::kTH1F, {{100, 0, 0.5}});
  }

  void processData(aod::Tracks const& tracks)
  {
    // Increment the event counter
    registry.fill(HIST("NEvents"), 1);

    // Convert O2 tracks to FastJet PseudoJets
    vector<PseudoJet> event;
    for (auto& track : tracks) {
      PseudoJet particle(track.px(), track.py(), track.pz(), track.energy(o2::constants::physics::MassPionCharged));
      event.push_back(particle);
    }
    LOG(info) << "Read an event with " << event.size() << " particles";

    // Define the jet definition and cluster the event
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(event, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));

    // Apply SoftDrop to each jet
    for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
      PseudoJet sd_jet = (*sd)(jets[ijet]);
      
      // Convert PseudoJet to string for logging
      std::stringstream ss_jet;
      ss_jet << jets[ijet];
      LOG(info) << "Original jet: " << ss_jet.str();

      std::stringstream ss_sd_jet;
      ss_sd_jet << sd_jet;
      LOG(info) << "SoftDropped jet: " << ss_sd_jet.str();

      assert(sd_jet != 0); // because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet

      // Fill histograms
      double DeltaR = sd_jet.structure_of<contrib::SoftDrop>().delta_R();
      double rg = DeltaR / R; // Normalisation 
      registry.fill(HIST("jet_pt"), jets[ijet].pt());
      registry.fill(HIST("jet_mass"), jets[ijet].m());
      registry.fill(HIST("sd_jet_pt"), sd_jet.pt());
      registry.fill(HIST("sd_jet_mass"), sd_jet.m());
      registry.fill(HIST("DeltaR"), DeltaR);
      registry.fill(HIST("hSymmetry"), sd_jet.structure_of<contrib::SoftDrop>().symmetry());
      registry.fill(HIST("hMu"), sd_jet.structure_of<contrib::SoftDrop>().mu());
      registry.fill(HIST("rg"), rg);  
      registry.fill(HIST("hJetPtRg"), jets[ijet].pt(), rg);

    }
  }
PROCESS_SWITCH(SoftDropTask, processData, "process task for data particles", true);



void processMC(soa::Join<aod::McParticles, aod::McCollisions> const& mcParticles)

  {
    // Increment the event counter
    registry.fill(HIST("NEvents"), 1);
    
    // Convert O2 MC particles to FastJet PseudoJets
    vector<PseudoJet> event;
    for (auto& particle : mcParticles) {
      PseudoJet pj(particle.px(), particle.py(), particle.pz(), particle.e());
      event.push_back(pj);
    }
    LOG(info) << "Read an event with " << event.size() << " particles";

    // Define the jet definition and cluster the event
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(event, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));

    // Apply SoftDrop to each jet
    for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
      PseudoJet sd_jet = (*sd)(jets[ijet]);
      // Convert PseudoJet to string for logging
      std::stringstream ss_jet;
      ss_jet << jets[ijet];
      LOG(info) << "Original jet: " << ss_jet.str();

      std::stringstream ss_sd_jet;
      ss_sd_jet << sd_jet;
      LOG(info) << "SoftDropped jet: " << ss_sd_jet.str();

      assert(sd_jet != 0); // because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet

      // Fill histograms
      double DeltaR = sd_jet.structure_of<contrib::SoftDrop>().delta_R();
      double rg = DeltaR / R; // Normalisation 
      registry.fill(HIST("jet_pt"), jets[ijet].pt());
      registry.fill(HIST("jet_mass"), jets[ijet].m());
      registry.fill(HIST("sd_jet_pt"), sd_jet.pt());
      registry.fill(HIST("sd_jet_mass"), sd_jet.m());
      registry.fill(HIST("DeltaR"), DeltaR);
      registry.fill(HIST("hSymmetry"), sd_jet.structure_of<contrib::SoftDrop>().symmetry());
      registry.fill(HIST("hMu"), sd_jet.structure_of<contrib::SoftDrop>().mu());
      registry.fill(HIST("rg"), rg);  
      registry.fill(HIST("hJetPtRg"), jets[ijet].pt(), rg);
    }
  }

PROCESS_SWITCH(SoftDropTask, processMC, "process task for MC particles", false);

contrib::SoftDrop* sd;

};


//----------------------------------------------------------------------
/// Overloaded jet info output
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt()
         << " m = " << jet.m()
         << " y = " << jet.rap()
         << " phi = " << jet.phi();
  })
  return ostr;
}

//----------------------------------------------------------------------
// Define the O2Physics workflow
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SoftDropTask>(cfgc)
  };
}
