// Tau analysis module by J. Hewes <jhewes15@fnal.gov>
// Modified by C. Sarasty <csarasty@fnal.gov>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Services/Optional/TFileService.h"
//LArSorft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
//ROOT includes 
#include "TH1.h"
#include "TTree.h"
///c++ includes
#include <vector>
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"



namespace dune {

  class TrueEnergy : public art::EDAnalyzer {

  public:

    explicit TrueEnergy(fhicl::ParameterSet const& pset);
    //virtual ~TrueEnergy();

    void analyze(const art::Event& evt);

  private:

    //TH1* fTrueEnergy;
    TTree* fTree;

    double fTrueE, fTvisible_E;
    double fTnu_directionX, fTnu_directionY, fTnu_directionZ;
    double fTvisible_directionX, fTvisible_directionY, fTvisible_directionZ;
    double fCalo_Energy;
  }; // class dune::TrueEnergy
} // namespace dune

dune::TrueEnergy::TrueEnergy(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset)
{
  // Get a handle to the TFile service
  art::ServiceHandle<art::TFileService> tfs;
  // Create tree 
   fTree = tfs->make<TTree>("TauTree", "TauTree");
   fTree->Branch("TrueE", &fTrueE, "TrueE/D");
   fTree->Branch("Tvisible_E",&fTvisible_E,"Tvisible_E/D");
   fTree->Branch("Tnu_directionX", &fTnu_directionX,"Tnu_directionX/D");
   fTree->Branch("Tnu_directionY", &fTnu_directionY,"Tnu_directionY/D");
   fTree->Branch("Tnu_directionZ", &fTnu_directionZ,"Tnu_directionZ/D");
   fTree->Branch("Tvisible_directionX",&fTvisible_directionX,"Tvisible_directionX/D");
   fTree->Branch("Tvisible_directionY",&fTvisible_directionY,"Tvisible_directionY/D");
   fTree->Branch("Tvisible_directionZ",&fTvisible_directionZ,"Tvisible_directionZ/D");
   fTree->Branch("CaloE", &fCalo_Energy, "CaloE/D");
} // dune::TrueEnergy::TrueEnergy()

void dune::TrueEnergy::analyze(const art:: Event& evt)
{
  //double EvtNum = evt.id().event();
  
  //if(EvtNum == 65864){
  
  // First, get an art handle to MC truth information
  art::Handle<std::vector<simb::MCTruth>> mct;
  std::vector<art::Ptr<simb::MCTruth>> truthlist;
  if (evt.getByLabel("generator", mct))
    art::fill_ptr_vector(truthlist, mct);
  simb::MCTruth truth(*(truthlist[0]));
 // std::cerr<<"[Carlos mct ]" << truthlist.size() << std::endl;
  // Get  an art handle to  MCShowers information  
  art::Handle<std::vector<sim::MCShower>> mcs;
  if (!evt.getByLabel("mcreco", mcs))
      throw std::runtime_error("Could not get MCShower information!");
  std::vector<sim::MCShower> const& mcshower(*mcs);
  //std::cerr<<"[Carlos showers ]" << mcshower.size() << std::endl;
 // Get an art handle to  MC track information 
 art::Handle<std::vector<sim::MCTrack>> mctr;
  if (!evt.getByLabel("mcreco", mctr))
    throw std::runtime_error("Could not get MCShower information!");
  std::vector<sim::MCTrack> const& mctrack(*mctr);
  //std::cerr<<"[Carlos tracks ]" << mctrack.size() << std::endl;
  //Calorimetric Energy
  art::Handle<dune::EnergyRecoOutput> handle_energy_reco_calo;
  evt.getByLabel("energynutau", handle_energy_reco_calo);
  fCalo_Energy = handle_energy_reco_calo->fNuLorentzVector.E();
 // std::cerr<< "[Carlos ECalo ] " << fCalo_Energy << std::endl;
  
  double fEnergy_Shower, fEnergy_track;
  TVector3 visible_direction(0, 0, 0);
  TVector3 visible_directionUnit(0, 0, 0);
//loop over MCShowers
  for(auto &k : mcshower){
    fEnergy_Shower += k.Start().E()-k.End().E();
    //MC shower momentum
    TVector3 start_momentum = k.Start().Momentum().Vect();
    TVector3 end_momentum = k.End().Momentum().Vect();
    visible_direction += start_momentum - end_momentum;
  }
//std::cerr<<"[Carlos Energy shower]" << fEnergy_Shower<<std::endl;

//loop over MCTracks 
  for(auto &k : mctrack){
    // Energy from tracks 
    fEnergy_track += k.Start().E()-k.End().E();
    // Track momentum 
    TVector3 start_momentum = k.Start().Momentum().Vect();
    TVector3 end_momentum = k.End().Momentum().Vect();
    visible_direction += start_momentum - end_momentum;
  }
//std::cerr<<"[Carlos Energy track ]" << fEnergy_track<<std::endl;
//visible true energy 
fTvisible_E = fEnergy_track + fEnergy_Shower;
//std::cerr<<"[Carlos Visible Energy ]" << fTvisible_E << std::endl;
//visible direction 
//std::cerr<<"[Carlos antes ] " << visible_direction.X() << std::endl;
visible_directionUnit = visible_direction.Unit();
//std::cerr<<"[Carlos Unitario ]" << visible_directionUnit.X() << std::endl;

fTvisible_directionX = visible_directionUnit.X();
fTvisible_directionY = visible_directionUnit.Y();
fTvisible_directionZ = visible_directionUnit.Z();

  simb::MCNeutrino nu = truthlist[0]->GetNeutrino();
  fTrueE = nu.Nu().E();
  //direction of the incoming neutrino 
  TVector3 nu_direction(truth.GetNeutrino().Nu().Momentum().Vect().Unit());
  fTnu_directionX = nu_direction.X();
  fTnu_directionY = nu_direction.Y();
  fTnu_directionZ = nu_direction.Z();
 // std::cerr<<"[Carlos neutrino direction x, y, z ] " << nu_direction.X()<< " "
  //<<nu_direction.X()<< " " << nu_direction.Z() << std::endl;  

  //Store info in tree
  fTree->Fill();
 // }  // if EvtNum
 
} // dune::TrueEnergy::analyze()

namespace dune {

  DEFINE_ART_MODULE(TrueEnergy)

}
