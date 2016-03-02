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
#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/ParticleCollection.h"
#include "datamodel/GenJetCollection.h"
#include "datamodel/JetCollection.h"
#include "datamodel/IntTagCollection.h"
#include "datamodel/TagCollection.h"
#include "datamodel/METCollection.h"
#include "datamodel/ParticleMCParticleAssociationCollection.h"
#include "datamodel/ParticleTagAssociationCollection.h"
#include "datamodel/GenJetParticleAssociationCollection.h"
#include "datamodel/GenJetIntTagAssociationCollection.h"
#include "datamodel/JetParticleAssociationCollection.h"
#include "datamodel/JetIntTagAssociationCollection.h"
#include "datamodel/JetTagAssociationCollection.h"

// ROOT
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"

/** @mainpage DelphesFCCAnalyzer
 *
 *  Short program to be used to read in the FCC-EDM based root file and process the collections
 *  saved in the file. The program uses standard boost input parameters to control input to the
 *  program. Use -h option to get the help.
 *
 *  @author: Z. Drasal (CERN)
 *
 */
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

    // Read FCC-EDM containers
    const fcc::MCParticleCollection* colGenParticles(nullptr);
    const fcc::GenVertexCollection*  colGenVertices(nullptr);
    const fcc::GenJetCollection*     colGenJets(nullptr);
    const fcc::IntTagCollection*     colGenJetsFlavor(nullptr);
    const fcc::ParticleCollection*   colRecMuons(nullptr);
    const fcc::TagCollection*        colRecITagMuons(nullptr);
    const fcc::ParticleCollection*   colRecElectrons(nullptr);
    const fcc::TagCollection*        colRecITagElectrons(nullptr);
    const fcc::ParticleCollection*   colRecCharged(nullptr);
    const fcc::ParticleCollection*   colRecNeutral(nullptr);
    const fcc::ParticleCollection*   colRecPhotons(nullptr);
    const fcc::TagCollection*        colRecITagPhotons(nullptr);
    const fcc::JetCollection*        colRecJets(nullptr);
    const fcc::ParticleCollection*   colRecJetParts(nullptr);
    const fcc::IntTagCollection*     colRecJetsFlavor(nullptr);
    const fcc::TagCollection*        colRecBTags(nullptr);
    const fcc::TagCollection*        colRecTauTags(nullptr);
    const fcc::METCollection*        colRecMETs(nullptr);

    const fcc::GenJetParticleAssociationCollection*     acolGenJetsToMC(nullptr);
    const fcc::GenJetIntTagAssociationCollection*       acolGenJetsToFlavor(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecMuonsToMC(nullptr);
    const fcc::ParticleTagAssociationCollection*        acolRecMuonsToITags(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecElectronsToMC(nullptr);
    const fcc::ParticleTagAssociationCollection*        acolRecElectronsToITags(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecChargedToMC(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecNeutralToMC(nullptr);
    const fcc::ParticleMCParticleAssociationCollection* acolRecPhotonsToMC(nullptr);
    const fcc::ParticleTagAssociationCollection*        acolRecPhotonsToITags(nullptr);
    const fcc::JetParticleAssociationCollection*        acolRecJetsToParts(nullptr);
    const fcc::JetIntTagAssociationCollection*          acolRecJetsToFlavor(nullptr);
    const fcc::JetTagAssociationCollection*             acolRecJetsToBTags(nullptr);
    const fcc::JetTagAssociationCollection*             acolRecJetsToTauTags(nullptr);

    bool colGenParticlesOK     = store.get("genParticles" , colGenParticles);
    bool colGenVerticesOK      = store.get("genVertices"  , colGenVertices);
    bool colGenJetsOK          = store.get("genJets"      , colGenJets);
    bool colGenJetsFlavorOK    = store.get("genJetsFlavor", colGenJetsFlavor);
    bool colRecMuonsOK         = store.get("muons"        , colRecMuons);
    bool colRecITagMuonsOK     = store.get("muonITags"    , colRecITagMuons);
    bool colRecElectronsOK     = store.get("electrons"    , colRecElectrons);
    bool colRecITagElectronsOK = store.get("electronITags", colRecITagElectrons);
    bool colRecChargedOK       = store.get("charged"      , colRecCharged);
    bool colRecNeutralOK       = store.get("neutral"      , colRecNeutral);
    bool colRecPhotonsOK       = store.get("photons"      , colRecPhotons);
    bool colRecITagPhotonsOK   = store.get("photonITags"  , colRecITagPhotons);
    bool colRecJetsOK          = store.get("jets"         , colRecJets);
    bool colRecJetPartsOK      = store.get("jetParts"     , colRecJetParts);
    bool colRecJetsFlavorOK    = store.get("jetsFlavor"   , colRecJetsFlavor);
    bool colRecBTagsOK         = store.get("bTags"        , colRecBTags);
    bool colRecTauTagsOK       = store.get("tauTags"      , colRecTauTags);
    bool colRecMETsOK          = store.get("met"          , colRecMETs);

    bool acolGenJetsToMCOK         = store.get("genJetsToMC"     , acolGenJetsToMC);
    bool acolGenJetsToFlavorOK     = store.get("genJetsToFlavor" , acolGenJetsToFlavor);
    bool acolRecMuonsToMCOK        = store.get("muonsToMC"       , acolRecMuonsToMC);
    bool acolRecMuonsToITagsOK     = store.get("muonsToITags"    , acolRecMuonsToITags);
    bool acolRecElectronsToMCOK    = store.get("electronsToMC"   , acolRecElectronsToMC);
    bool acolRecElectronsToITagsOK = store.get("electronsToITags", acolRecElectronsToITags);
    bool acolRecChargedToMCOK      = store.get("chargedToMC"     , acolRecChargedToMC);
    bool acolRecNeutralToMCOK      = store.get("neutralToMC"     , acolRecNeutralToMC);
    bool acolRecPhotonsToMCOK      = store.get("photonsToMC"     , acolRecPhotonsToMC);
    bool acolRecPhotonsToITagsOK   = store.get("photonsToITags"  , acolRecPhotonsToITags);
    bool acolRecJetsToPartsOK      = store.get("jetsToParts"     , acolRecJetsToParts);
    bool acolRecJetsToFlavorOK     = store.get("jetsToFlavor"    , acolRecJetsToFlavor);
    bool acolRecJetsToBTagsOK      = store.get("jetsToBTags"     , acolRecJetsToBTags);
    bool acolRecJetsToTauTagsOK    = store.get("jetsToTauTags"   , acolRecJetsToTauTags);

    if (colGenParticlesOK &&
        colGenVerticesOK  &&
        colGenJetsOK      && acolGenJetsToMCOK     && colGenJetsFlavorOK     && acolGenJetsToFlavorOK     &&
        colRecMuonsOK     && colRecITagMuonsOK     && acolRecMuonsToMCOK     && acolRecMuonsToITagsOK     &&
        colRecElectronsOK && colRecITagElectronsOK && acolRecElectronsToMCOK && acolRecElectronsToITagsOK &&
        colRecChargedOK   && acolRecChargedToMCOK  &&
        colRecNeutralOK   && acolRecNeutralToMCOK  &&
        colRecPhotonsOK   && colRecITagPhotonsOK   && acolRecPhotonsToMCOK   && acolRecPhotonsToITagsOK   &&
        colRecJetsOK      && colRecJetPartsOK      && colRecJetsFlavorOK     && colRecBTagsOK             && colRecTauTagsOK           &&
                             acolRecJetsToPartsOK  && acolRecJetsToFlavorOK  && acolRecJetsToBTagsOK      && acolRecJetsToTauTagsOK    &&
        colRecMETsOK ) {

      std::cout << "Event: "                 << iEvent << std::endl;
      std::cout << " Collections: "          << std::endl;
      std::cout << " -> #GenParticles:     " << colGenParticles->size()    << std::endl;
      std::cout << " -> #GenVertices:      " << colGenVertices->size()     << std::endl;
      std::cout << " -> #GenJets:          " << colGenJets->size()         << std::endl;
      std::cout << " -> #GenJetsFlavor:    " << colGenJetsFlavor->size()   << std::endl;
      std::cout << " -> #RecMuons:         " << colRecMuons->size()        << std::endl;
      std::cout << " -> #RecITagMuons:     " << colRecITagMuons->size()    << std::endl;
      std::cout << " -> #RecElectrons:     " << colRecElectrons->size()    << std::endl;
      std::cout << " -> #RecITagElectrons: " << colRecITagElectrons->size()<< std::endl;
      std::cout << " -> #RecCharged:       " << colRecCharged->size()      << std::endl;
      std::cout << " -> #RecNeutral:       " << colRecNeutral->size()      << std::endl;
      std::cout << " -> #RecPhotons:       " << colRecPhotons->size()      << std::endl;
      std::cout << " -> #RecITagPhotons:   " << colRecITagPhotons->size()  << std::endl;
      std::cout << " -> #RecJets:          " << colRecJets->size()         << std::endl;
      std::cout << " -> #RecJetParts:      " << colRecJetParts->size()     << std::endl;
      std::cout << " -> #RecJetsFlavor:    " << colRecJetsFlavor->size()   << std::endl;
      std::cout << " -> #RecBTags:         " << colRecBTags->size()        << std::endl;
      std::cout << " -> #RecTuaTags:       " << colRecTauTags->size()      << std::endl;
      std::cout << " -> #RecMETs:          " << colRecMETs->size()         << std::endl;

      std::cout << " Relations: "               << std::endl;
      std::cout << " -> #RelGenJetsToMC:      " << acolGenJetsToMC->size()         << std::endl;
      std::cout << " -> #RelGenJetsToFlavor:  " << acolGenJetsToFlavor->size()     << std::endl;
      std::cout << " -> #RelMuonsToMC:        " << acolRecMuonsToMC->size()        << std::endl;
      std::cout << " -> #RelMuonsToITags:     " << acolRecMuonsToITags->size()     << std::endl;
      std::cout << " -> #RelElectronsToMC:    " << acolRecElectronsToMC->size()    << std::endl;
      std::cout << " -> #RelElectronsToITags: " << acolRecElectronsToITags->size() << std::endl;
      std::cout << " -> #RelChargedToMC:      " << acolRecChargedToMC->size()      << std::endl;
      std::cout << " -> #RelNeutronsToMC:     " << acolRecNeutralToMC->size()      << std::endl;
      std::cout << " -> #RelPhotonsToMC:      " << acolRecPhotonsToMC->size()      << std::endl;
      std::cout << " -> #RelPhotonsToITags:   " << acolRecPhotonsToITags->size()   << std::endl;
      std::cout << " -> #RelJetsToParts:      " << acolRecJetsToParts->size()      << std::endl;
      std::cout << " -> #RelJetsToFlavor:     " << acolRecJetsToFlavor->size()     << std::endl;
      std::cout << " -> #RelJetsToBTags:      " << acolRecJetsToBTags->size()      << std::endl;
      std::cout << " -> #RelJetsToTauTags:    " << acolRecJetsToTauTags->size()    << std::endl;

      std::cout << std::endl;
//      std::cout << "GenParticles: " << std::endl;
//      int idPart = 0;
//      for (auto& iPart=colGenParticles->begin(); iPart!=colGenParticles->end(); ++iPart) {
//
//        idPart++;
//        double partE = sqrt(iPart->Core().P4.Px  *iPart->Core().P4.Px +
//                            iPart->Core().P4.Py  *iPart->Core().P4.Py +
//                            iPart->Core().P4.Pz  *iPart->Core().P4.Pz +
//                            iPart->Core().P4.Mass*iPart->Core().P4.Mass);
//
//        std::cout << " MCParticle: "
//                  << " Id: "       << std::setw(3)  << idPart
//                  << " Pdg: "      << std::setw(5)  << iPart->Core().Type
//                  << " Stat: "     << std::setw(2)  << iPart->Core().Status
//                  << " Bits: "     << std::setw(2)  << iPart->Core().Bits
//                  << std::scientific
//                  << " Px: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Px
//                  << " Py: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Py
//                  << " Pz: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Pz
//                  << " E: "        << std::setprecision(2) << std::setw(9) << partE
//                  << " M: "        << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Mass;
//        if (iPart->StartVertex().isAvailable()) {
//          std::cout << " Vx: "       << std::setprecision(2) << std::setw(9) << iPart->StartVertex().Position().X
//                    << " Vy: "       << std::setprecision(2) << std::setw(9) << iPart->StartVertex().Position().Y
//                    << " Vz: "       << std::setprecision(2) << std::setw(9) << iPart->StartVertex().Position().Z
//                    << " T: "        << std::setprecision(2) << std::setw(9) << iPart->StartVertex().Ctau();
//        }
//        if (iPart->EndVertex().isAvailable()) {
//          std::cout << " Vx: "       << std::setprecision(2) << std::setw(9) << iPart->EndVertex().Position().X
//                    << " Vy: "       << std::setprecision(2) << std::setw(9) << iPart->EndVertex().Position().Y
//                    << " Vz: "       << std::setprecision(2) << std::setw(9) << iPart->EndVertex().Position().Z
//                    << " T: "        << std::setprecision(2) << std::setw(9) << iPart->EndVertex().Ctau();
//        } // GenParticles
//        std::cout << std::endl;
//      }
//      std::cout << std::endl;
//
//      std::cout << "GenJets: " << std::endl;
//      idPart = 0;
//      for (auto& iPart=colGenJets->begin(); iPart!=colGenJets->end(); ++iPart) {
//
//        idPart++;
//        double recE = sqrt(iPart->Core().P4.Px  *iPart->Core().P4.Px +
//                           iPart->Core().P4.Py  *iPart->Core().P4.Py +
//                           iPart->Core().P4.Pz  *iPart->Core().P4.Pz +
//                           iPart->Core().P4.Mass*iPart->Core().P4.Mass);
//
//        std::cout << " GenJet: "
//                  << " Id: "       << std::setw(3)  << idPart
//                  << " Bits: "     << std::setw(2)  << iPart->Core().Bits
//                  << std::scientific
//                  << " Px: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Px
//                  << " Py: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Py
//                  << " Pz: "       << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Pz
//                  << " E: "        << std::setprecision(2) << std::setw(9) << recE
//                  << " M: "        << std::setprecision(2) << std::setw(9) << iPart->Core().P4.Mass;
//        std::cout << std::endl;
//
//        // Relation jet to MC
//        double totSimE= 0;
//        int    refID  = 0;
//        for (auto& iRel=acolGenJetsToMC->begin(); iRel!=acolGenJetsToMC->end(); ++iRel) {
//
//          if (iRel->Jet().getObjectID().index==iPart->getObjectID().index) {
//
//            refID++;
//            double simE   = sqrt(iRel->Particle().Core().P4.Px  *iRel->Particle().Core().P4.Px +
//                                 iRel->Particle().Core().P4.Py  *iRel->Particle().Core().P4.Py +
//                                 iRel->Particle().Core().P4.Pz  *iRel->Particle().Core().P4.Pz +
//                                 iRel->Particle().Core().P4.Mass*iRel->Particle().Core().P4.Mass);
//            totSimE += simE;
//            std::cout << "  RefId: " << std::setw(3)            << refID
//                      << " Rel E: "  << std::setprecision(2)
//                                     << std::scientific
//                                     << std::setw(9) << simE    << " "
//                                     << std::setw(9) << totSimE << " <-> "
//                                     << std::setw(9) << recE
//                                     << std::fixed
//                                     << std::endl;
//          }
//        } // Relation jets to MC
//
//        // Relation jet to BTag
//        refID  = 0;
//        std::cout << std::endl;
//        for (auto& iRel=acolRecJetsToBTags->begin(); iRel!=acolRecJetsToBTags->end(); ++iRel) {
//
//          if (iRel->Jet().getObjectID().index==iPart->getObjectID().index) {
//
//            refID++;
//            std::cout << "  RefId: " << std::setw(3)         << refID
//                      << " BTag: "   << std::setprecision(2) << iRel->Tag().Value()
//                                     << std::endl;
//          }
//        } // Relation jets to BTag
//
//        // Relation jet to TauTag
//        refID  = 0;
//        std::cout << std::endl;
//        for (auto& iRel=acolRecJetsToTauTags->begin(); iRel!=acolRecJetsToTauTags->end(); ++iRel) {
//
//          if (iRel->Jet().getObjectID().index==iPart->getObjectID().index) {
//
//            refID++;
//            std::cout << "  RefId: " << std::setw(3)         << refID
//                      << " TauTag: " << std::setprecision(2) << iRel->Tag().Value()
//                                     << std::endl;
//          }
//        } // Relation jets to TauTag
//        std::cout << std::endl;
//
//      } // RecJets
    }
    else {

      if (!colGenParticlesOK)         std::cout << "Missing genParticles collection"                 << std::endl;
      if (!colGenVerticesOK)          std::cout << "Missing genVertices collection"                  << std::endl;
      if (!colGenJetsOK)              std::cout << "Missing genJets collection"                      << std::endl;
      if (!colGenJetsFlavorOK)        std::cout << "Missing genJetsFlavor collection"                << std::endl;
      if (!acolGenJetsToMCOK)         std::cout << "Missing genJetsToMC association collection"      << std::endl;
      if (!acolGenJetsToFlavorOK)     std::cout << "Missing genJetsToFlavor association collection"  << std::endl;
      if (!colRecMuonsOK)             std::cout << "Missing muons collection"                        << std::endl;
      if (!colRecITagMuonsOK)         std::cout << "Missing muonITags collection"                    << std::endl;
      if (!acolRecMuonsToMCOK)        std::cout << "Missing muonsToMC association collection"        << std::endl;
      if (!acolRecMuonsToITagsOK)     std::cout << "Missing muonsToITags association collection"     << std::endl;
      if (!colRecElectronsOK)         std::cout << "Missing electrons collection"                    << std::endl;
      if (!colRecITagElectronsOK)     std::cout << "Missing electronITags collection"                << std::endl;
      if (!acolRecElectronsToMCOK)    std::cout << "Missing electronsToMC association collection"    << std::endl;
      if (!acolRecElectronsToITagsOK) std::cout << "Missing electronsToITags association collection" << std::endl;
      if (!colRecChargedOK)           std::cout << "Missing charged collection"                      << std::endl;
      if (!acolRecChargedToMCOK)      std::cout << "Missing chargedToMC association collection"      << std::endl;
      if (!colRecNeutralOK)           std::cout << "Missing neutral collection"                      << std::endl;
      if (!acolRecNeutralToMCOK)      std::cout << "Missing neutralToMC association collection"      << std::endl;
      if (!colRecPhotonsOK)           std::cout << "Missing photons collection"                      << std::endl;
      if (!colRecITagPhotonsOK)       std::cout << "Missing photonITags collection"                  << std::endl;
      if (!acolRecPhotonsToMCOK)      std::cout << "Missing photonsToMC association collection"      << std::endl;
      if (!acolRecPhotonsToITagsOK)   std::cout << "Missing photonsToITags association collection"   << std::endl;
      if (!colRecJetsOK)              std::cout << "Missing jets collection"                         << std::endl;
      if (!colRecJetPartsOK)          std::cout << "Missing jetParts collection"                     << std::endl;
      if (!colRecJetsFlavorOK)        std::cout << "Missing jetsFlavor collection"                   << std::endl;
      if (!colRecBTagsOK)             std::cout << "Missing bTags collection"                        << std::endl;
      if (!colRecTauTagsOK)           std::cout << "Missing tauTags collection"                      << std::endl;
      if (!acolRecJetsToPartsOK)      std::cout << "Missing jetsToParts association collection"      << std::endl;
      if (!acolRecJetsToFlavorOK)     std::cout << "Missing jetsToFlavor association collection"     << std::endl;
      if (!acolRecJetsToBTagsOK)      std::cout << "Missing jetsToBTags association collection"      << std::endl;
      if (!acolRecJetsToTauTagsOK)    std::cout << "Missing jetsToTauTags association collection"    << std::endl;
      if (!colRecMETsOK)              std::cout << "Missing met collection"                          << std::endl;
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
