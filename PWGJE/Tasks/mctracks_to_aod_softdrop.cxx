// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include <iomanip>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::math;
using namespace fastjet; 
using namespace std;

// Fonctions utilitaires pour l'affichage
void print_prongs_with_clustering_info(const PseudoJet &jet, const string &pprefix);
void print_raw_prongs(const PseudoJet &jet);
ostream & operator<<(ostream &, const PseudoJet &);

struct MyCustomTask {

  HistogramRegistry registry{"registry", {}};

  void init(o2::framework::InitContext& /*ic*/)
  {
    registry.add<TH1>("NEvents", "NEvents", HistType::kTH1F, {{1, 0, 1}}, false);
    registry.add<TH1>("pt", "pt", HistType::kTH1F, {{100, 0, 100}}, false);
    registry.add<TH1>("eta", "eta", HistType::kTH1F, {{100, -1, 1}}, false);
    registry.add<TH1>("phi", "phi", HistType::kTH1F, {{100, 0, 2 * PI}}, false);
    registry.add<TH1>("jet_pt", "jet pt", HistType::kTH1F, {{100, 0, 100}}, false);
    registry.add<TH1>("jet_mass", "jet mass", HistType::kTH1F, {{100, 0, 20}}, false);
    registry.add<TH1>("jet_eta", "jet eta", HistType::kTH1F, {{100, -1, 1}}, false); // peut etre pas besoin 
    registry.add<TH1>("jet_phi", "jet phi", HistType::kTH1F, {{100, 0, 2 * PI}}, false);// peut etre pas besoin 
    registry.add<TH1>("rsd_jet_pt", "RSD jet pt", HistType::kTH1F, {{100, 0, 100}}, false);
    registry.add<TH1>("rsd_jet_mass", "RSD jet mass", HistType::kTH1F, {{100, 0, 20}}, false);
    registry.add<TH1>("zg", "zg", HistType::kTH1F, {{100, 0, 0.5}}, false);
    registry.add<TH1>("thetag", "thetag", HistType::kTH1F, {{100, 0, 0.5}}, false);
    registry.add<TH1>("kT", "kT", HistType::kTH1F, {{100, 0, 10}}, false);
    registry.add<TH2>("h2_lnkt_vs_lnthetag", "ln(kT) vs ln(1/#theta_{g}); ln(1/#theta_{g}); ln(kT)", HistType::kTH2F, {{100, -0.5, 7}, {100, -5, 7}},false);
  }

  void process(aod::McCollision const&, aod::McParticles const& mcParticles)
  {
    registry.fill(HIST("NEvents"), 0.5);
    
    // Convert MC particles to FastJet format
    vector<PseudoJet> particles;
    for (auto& mcparticle : mcParticles) {
      
      registry.fill(HIST("pt"), mcparticle.pt());
      registry.fill(HIST("eta"), mcparticle.eta());
      registry.fill(HIST("phi"), mcparticle.phi());
      
      // Create PseudoJet from MC particle
      particles.emplace_back(mcparticle.px(), mcparticle.py(), mcparticle.pz(), mcparticle.e());//rempli le pseudojet avec px py pz e
      particles.back().set_user_index(mcparticle.globalIndex()); //Attache l'index original de la particule MC au PseudoJet
    }

    if (particles.size() == 0) return; //Vérifie si le vecteur particles est vide après la conversion des mcParticles en PseudoJet

    // Jet definition
    double R = 1.0;
    double ptmin = 100.0;
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));

    // Recursive Soft Drop parameters
    double z_cut = 0.4;
    double beta = 0.2;
    int n = 4; // infinite recursion
    
    contrib::RecursiveSoftDrop rsd(beta, z_cut, n, R);
    rsd.set_verbose_structure(true);
    rsd.set_dynamical_R0();

    // Print infos in Recursive SD
    LOGF(info, "RecursiveSoftDrop groomer is: %s", rsd.description().c_str());
    LOGF(info, "Processing event with %d particles", particles.size());

    for (const auto& jet : jets) {
      registry.fill(HIST("jet_pt"), jet.pt());
      registry.fill(HIST("jet_mass"), jet.m());
      registry.fill(HIST("jet_eta"), jet.eta());
      registry.fill(HIST("jet_phi"), jet.phi_std());

      // Apply Recursive Soft Drop
      PseudoJet rsd_jet = rsd(jet);
      
      if (rsd_jet == 0) continue; // skip if grooming failed
      
      registry.fill(HIST("rsd_jet_pt"), rsd_jet.pt());
      registry.fill(HIST("rsd_jet_mass"), rsd_jet.m());

      // Print of jet information 
      std::stringstream ss_jet;
      ss_jet << jet;
      LOGF(info, "\nOriginal jet: %s", ss_jet.str().c_str());
      std::stringstream ss_rsd_jet;
      ss_rsd_jet << rsd_jet;
      LOGF(info, "RecursiveSoftDropped jet: %s", ss_rsd_jet.str().c_str());
      
      // Print of prongs structure 
      LOGF(info, "\nProngs with clustering information\n----------------------------------");
      print_prongs_with_clustering_info(rsd_jet, " ");
      
      LOGF(info, "\nProngs without clustering information\n-------------------------------------");
      print_raw_prongs(rsd_jet);

      // Get grooming information
      const auto& rsd_struct = rsd_jet.structure_of<contrib::RecursiveSoftDrop>();
      vector<pair<double, double>> ztg = rsd_struct.sorted_zg_and_thetag(); //ztg is a vector with pair (zg, θg)

      LOGF(info, "Groomed prongs information:");
      LOGF(info, "index            zg        thetag         pT         kT");
      
      for (unsigned int i=0; i<ztg.size(); ++i) {
        double kT = ztg[i].first * ztg[i].second * rsd_jet.pt();
        LOGF(info, "%5d%14.4f%14.4f%14.4f%14.4f", 
          i+1, ztg[i].first, ztg[i].second, rsd_jet.pt(), kT);
        registry.fill(HIST("zg"), ztg[i].first); //first mean zg
        registry.fill(HIST("thetag"), ztg[i].second); //second mean thetag
        registry.fill(HIST("kT"), kT);

        if (ztg[i].second > 0 && kT > 0) { // Vérification to avoid ln(0) or ln(negatif)
          double ln_kt = TMath::Log(kT);
          double ln_inv_thetag = TMath::Log(1./ztg[i].second);
          registry.fill(HIST("h2_lnkt_vs_lnthetag"), ln_inv_thetag, ln_kt);
        }

      }
    }
  }
};

// print the prongs inside the jet, showing the clustering info
void print_prongs_with_clustering_info(const PseudoJet &jet, const string &prefix){
  stringstream output;
  if (prefix.size() == 1){
    output << " " << setw(14) << " "
           << setw(8) << "branch" << setw(14) << "branch"
           << setw(10) << "N_groomed"
           << setw(11) << "max loc"
           << setw(22) << "substructure" << endl;
    output << " " << setw(14) << " "
           << setw(8) << "pt" << setw(14) << "mass"
           << setw(5) << "loc"
           << setw(5) << "tot"
           << setw(11) << "zdrop"
           << setw(11) << "zg"
           << setw(11) << "thetag"<< endl;
  }
  const contrib::RecursiveSoftDrop::StructureType &structure = jet.structure_of<contrib::RecursiveSoftDrop>();
  double dR = structure.delta_R();
  output << " " << left << setw(14) << (prefix.substr(0, prefix.size()-1)+"+--> ") << right
         << setw(8) << jet.pt() << setw(14) << jet.m()
         << setw(5) << structure.dropped_count(false)
         << setw(5) << structure.dropped_count()
         << setw(11) << structure.max_dropped_symmetry(false);
  
  if (structure.has_substructure()){
    output << setw(11) << structure.symmetry()
           << setw(11) << structure.delta_R();
  }
  output << endl;
  
  if (dR>=0){
    vector<PseudoJet> pieces = jet.pieces();
    assert(pieces.size()==2);
    print_prongs_with_clustering_info(pieces[0], prefix+" |");
    print_prongs_with_clustering_info(pieces[1], prefix+"  ");
  }
  LOGF(info, "%s", output.str().c_str());
}

//----------------------------------------------------------------------
// print all the prongs inside the jet (no clustering info)
void print_raw_prongs(const PseudoJet &jet){
  stringstream output;
  output << "(Raw) list of prongs:" << endl;
  if (!jet.has_structure_of<contrib::RecursiveSoftDrop>()){
    output << "  None (bad structure)" << endl;
    LOGF(info, "%s", output.str().c_str());
    return;
  }
  
  output << setw(5) << " " << setw(11) << "pt" << setw(14) << "mass" << endl;

  vector<PseudoJet> prongs = contrib::recursive_soft_drop_prongs(jet);
  for (unsigned int iprong=0; iprong<prongs.size(); ++iprong){
    const PseudoJet & prong = prongs[iprong];
    const contrib::RecursiveSoftDrop::StructureType &structure = prong.structure_of<contrib::RecursiveSoftDrop>();
    output << setw(5) << iprong << setw(11) << prong.pt() << setw(14) << prong.m() << endl;
  
    assert(!structure.has_substructure());
  }
  output << endl;
  LOGF(info, "%s", output.str().c_str());
}

//----------------------------------------------------------------------
/// overloaded jet info output
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

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<MyCustomTask>(cfgc)};
}
