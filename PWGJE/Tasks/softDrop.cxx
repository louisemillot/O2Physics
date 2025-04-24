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

using namespace std;
using namespace fastjet;
using namespace o2::framework;

// Forward declaration to make things clearer
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
struct SoftDropTask {
  // Define the input and output data
  Inputs inputs{{"inputTracks", "TRACKS", "INPUT TRACKS"}};

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
  }

  void process(aod::Tracks const& tracks)
  {
    // Convert O2 tracks to FastJet PseudoJets
    vector<PseudoJet> event;
    for (auto& track : tracks) {
      PseudoJet particle(track.px(), track.py(), track.pz(), track.e());
      event.push_back(particle);
    }
    cout << "# read an event with " << event.size() << " particles" << endl;

    // Define the jet definition and cluster the event
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(event, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));

    // Apply SoftDrop to each jet
    for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
      PseudoJet sd_jet = (*sd)(jets[ijet]);
      cout << endl;
      cout << "original    jet: " << jets[ijet] << endl;
      cout << "SoftDropped jet: " << sd_jet << endl;

      assert(sd_jet != 0); // because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet

      cout << "  delta_R between subjets: " << sd_jet.structure_of<contrib::SoftDrop>().delta_R() << endl;
      cout << "  symmetry measure(z):     " << sd_jet.structure_of<contrib::SoftDrop>().symmetry() << endl;
      cout << "  mass drop(mu):           " << sd_jet.structure_of<contrib::SoftDrop>().mu() << endl;
    }
  }

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
  }
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
