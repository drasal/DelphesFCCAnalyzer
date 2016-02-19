//////////////////////////////////////////////////////////////////////////////////////
//
// Analyze simulation done by FCCSW: Pythia + Delphes
//
//  Z. Drasal, CERN
//
//  February 15th 2016
//
//////////////////////////////////////////////////////////////////////////////////////
#include "Global_constants.h"

// System
#include <cmath>
#include <string>
#include <vector>

// BOOST
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

// PODIO
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

// FCCEDM
#include "datamodel/MCParticle.h"
#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertex.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/Particle.h"
#include "datamodel/ParticleCollection.h"
#include "datamodel/GenJet.h"
#include "datamodel/GenJetCollection.h"
#include "datamodel/Tag.h"
#include "datamodel/TagCollection.h"
#include "datamodel/MET.h"
#include "datamodel/METCollection.h"
#include "datamodel/ParticleMCParticleAssociation.h"
#include "datamodel/GenJetParticleAssociation.h"
#include "datamodel/GenJetTagAssociation.h"

// ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"

int main(int argc, char* argv[]) {

  std::string usage("Usage: ");
  usage += argv[0];
  usage += " [options]";

  std::string inFileName  = "";
  std::string outFileName = "";

  unsigned nEvents     = 0;

  //
  // Program options
  boost::program_options::options_description help("Program options");
  help.add_options()
    ("help,h"        , "Display help")
    ("input-file,if" , boost::program_options::value<std::string>(), "Specify name of input FCC-EDM ROOT file")
    ("output-file,of", boost::program_options::value<std::string>(), "Specify name of output ROOT file, where results of analysis will be saved)")
    ("n-events,n"    , boost::program_options::value<int>()        , "Specify number of processed events, n<=number of events in the input file (optional)")
    ;

  // Read user input
  boost::program_options::variables_map varMap;
  try {

    // Parse user defined options
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(help).run(), varMap);
    boost::program_options::notify(varMap);

    // Check
    if      (!varMap.count("input-file") && !varMap.count("help")) throw boost::program_options::error("Forgot to define FCC-EDM ROOT input file???");
    else {

      // Write help
      if (varMap.count("help")) {
        std::cout << usage << std::endl << help << std::endl;
        return 0;
      }

      // Open FCC input ROOT file
      inFileName = varMap["input-file"].as<std::string>();

      std::vector<std::string> infos;
      boost::algorithm::split(infos,  inFileName, boost::algorithm::is_any_of("."));

      // Check that ROOT file given
      if (infos[infos.size()-1]=="root" ) {

        std::cout << "> Will process input-file: " << inFileName << std::endl;
      }
      else {

        std::cerr << "Provide ROOT file on input!!!" << std::endl << std::endl;
        std::cout << usage                           << std::endl << help << std::endl;
        return EXIT_FAILURE;
      }

      // Number of events
      if (varMap.count("n-events") && varMap["n-events"].as<int>()>0) nEvents = varMap["n-events"].as<int>();

      // ROOT output file
      if (varMap.count("output-file")) {

        outFileName = varMap["output-file"].as<std::string>();
        std::cout << "> Output will be written to: " << outFileName << std::endl;
      }
      else {
        std::cout << "> No output file provided!!! " << std::endl;
      }
    }
  } catch(boost::program_options::error& e) {

    // Display error type
    std::cerr << "\nERROR: " << e.what() << std::endl << std::endl;
    std::cout << usage                   << std::endl << help << std::endl;
    return EXIT_FAILURE;
  }

  // Analyze FCC
  auto reader = podio::ROOTReader();
  auto store  = podio::EventStore();

  // Set reader
  try {
    reader.openFile(inFileName);
  }
  catch(std::runtime_error& e) {

    std::cerr << "\nERROR: " << e.what() << std::endl << std::endl;
    std::cout << "Couldn't open the input file: " << inFileName << std::endl;
    return EXIT_FAILURE;
  }
  store.setReader(&reader);

  // Set n events
  if (reader.getEntries()<=nEvents) nEvents = reader.getEntries();

  // Open output file
  TFile* outFile = nullptr;
  if (outFileName!="") outFile = new TFile(outFileName.c_str(), "RECREATE");

  // Process events
  for (unsigned iEvent=0; iEvent<nEvents; ++iEvent) {

    if(iEvent%100==0 && iEvent!=0) std::cout << "Reading event: "<< iEvent << std::endl;

    // Read MC particles
    const fcc::MCParticleCollection* colGenParticles(nullptr);
    const fcc::GenVertexCollection*  colGenVertices(nullptr);
    const fcc::ParticleCollection*   colRecMuons(nullptr);
    const fcc::ParticleCollection*   colRecElectrons(nullptr);
    const fcc::ParticleCollection*   colRecCharged(nullptr);
    const fcc::ParticleCollection*   colRecNeutral(nullptr);
    const fcc::ParticleCollection*   colRecPhotons(nullptr);
    const fcc::GenJetCollection*     colRecJets(nullptr);
    const fcc::TagCollection*        colRecBTags(nullptr);
    const fcc::TagCollection*        colRecTauTags(nullptr);
    const fcc::METCollection*        colRecMETs(nullptr);

    const fcc::ParticleMCParticleAssociationCollection* acolRecMuonsToMC(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecElectronsToMC(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecChargedToMC(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecNeutralToMC(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecPhotonsToMC(nullptr);
    const fcc::GenJetParticleAssociationCollection*     acolRecJetsToMC(nullptr);
    const fcc::GenJetTagAssociationCollection*          acolRecJetsToBTags(nullptr);
    const fcc::GenJetTagAssociationCollection*          acolRecJetsToTauTags(nullptr);

    bool colGenParticlesOK = store.get("genParticles" , colGenParticles);
    bool colGenVerticesOK  = store.get("genVertices"  , colGenVertices);
    bool colRecMuonsOK     = store.get("recMuons"     , colRecMuons);
    bool colRecElectronsOK = store.get("recElectrons" , colRecElectrons);
    bool colRecChargedOK   = store.get("recCharged"   , colRecCharged);
    bool colRecNeutralOK   = store.get("recNeutral"   , colRecNeutral);
    bool colRecPhotonsOK   = store.get("recPhotons"   , colRecPhotons);
    bool colRecJetsOK      = store.get("recJets"      , colRecJets);
    bool colRecBTagsOK     = store.get("recBTags"     , colRecBTags);
    bool colRecTauTagsOK   = store.get("recTauTags"   , colRecTauTags);
    bool colRecMETsOK      = store.get("recMETs"      , colRecMETs);

//    bool acolRecMuonsToMCOK     = store.get("recMuonsToMC"     , acolRecMuonsToMC);
//    bool acolRecElectronsToMCOK = store.get("recElectronsToMC" , acolRecElectronsToMC);
//    bool acolRecChargedToMCOK   = store.get("recChargedToMC"   , acolRecElectronsToMC);
//    bool acolRecNeutralToMCOK   = store.get("recNeutralToMC"   , acolRecElectronsToMC);
//    bool acolRecPhotonsToMCOK   = store.get("recPhotonsToMC"   , acolRecElectronsToMC);
//    bool acolRecJetsToMCOK      = store.get("recJetsToMC"      , acolRecElectronsToMC);
//    bool acolRecJetsToBTagsOK   = store.get("recJetsToBTags"   , acolRecElectronsToMC);
//    bool acolRecJetsToTauTagsOK = store.get("recJetsToTauTags" , acolRecElectronsToMC);

    if (colGenParticlesOK &&
        colGenVerticesOK &&
        colRecMuonsOK &&
        colRecElectronsOK &&
        colRecChargedOK &&
        //colRecNeutralOK &&
        colRecPhotonsOK &&
        colRecJetsOK &&
        colRecBTagsOK &&
        colRecTauTagsOK &&
        colRecMETsOK //&&
//        acolRecMuonsToMCOK &&
//        acolRecElectronsToMCOK &&
//        acolRecChargedToMCOK &&
//        acolRecNeutralToMCOK &&
//        acolRecPhotonsToMCOK &&
//        acolRecJetsToMCOK &&
//        acolRecJetsToBTagsOK &&
//        acolRecJetsToTauTagsOK
        ) {

      std::cout << "Event: "             << iEvent << std::endl;
      std::cout << " -> #GenParticles: " << colGenParticles->size()<< std::endl;
      std::cout << " -> #GenVertices:  " << colGenVertices->size() << std::endl;
      std::cout << " -> #RecMuons:     " << colRecMuons->size()    << std::endl;
      std::cout << " -> #RecElectrons: " << colRecElectrons->size()<< std::endl;
      std::cout << " -> #RecCharged:   " << colRecCharged->size()  << std::endl;
      //std::cout << " -> #RecNeutral:   " << colRecNeutral->size()  << std::endl;
      std::cout << " -> #RecJets:      " << colRecJets->size()     << std::endl;
      std::cout << " -> #RecBTags:     " << colRecBTags->size()    << std::endl;
      std::cout << " -> #RecTuaTags:   " << colRecTauTags->size()  << std::endl;
      std::cout << " -> #RecMETs:      " << colRecMETs->size()     << std::endl;

      //std::cout << " -> #RelMuons:     " << acolRecMuonsToMC->size()     << std::endl;

      std::cout << std::endl;
      std::cout << "GenParticles: " << std::endl;
      int idPart = 0;
      for (auto& iPart=colGenParticles->begin(); iPart!=colGenParticles->end(); ++iPart) {

        idPart++;
        double partE = sqrt(iPart->Core().P4.Px  *iPart->Core().P4.Px +
                            iPart->Core().P4.Py  *iPart->Core().P4.Py +
                            iPart->Core().P4.Pz  *iPart->Core().P4.Pz +
                            iPart->Core().P4.Mass*iPart->Core().P4.Mass);

        std::cout << " MCParticle: "
                  << " Id: "       << std::setw(3)  << idPart
                  << " Pdg: "      << std::setw(5)  << iPart->Core().Type
                  << " Stat: "     << std::setw(2)  << iPart->Core().Status
                  << " Bits: "     << std::setw(2)  << iPart->Core().Bits
                  << std::scientific
                  << " Px: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Px
                  << " Py: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Py
                  << " Pz: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Pz
                  << " E: "        << std::setprecision(2) << std::setw(9) << partE
                  << " M: "        << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Mass;
        std::cout << std::endl;
      }
      std::cout << std::endl;

    }
    else {

      std::cout << "Some collection not available" << std::endl;
    }

    // Prepare for the next event
    store.clear();
    reader.endOfEvent();
  }

  // Close output file
  if (outFile) {

    outFile->Close();
    delete outFile;
    outFile=nullptr;
  }

  // OK
  return EXIT_SUCCESS;
}
