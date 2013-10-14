// -*- C++ -*-
//
// Package:    OccupancyAnalyzer
// Class:      OccupancyAnalyzer
// 
/**\class OccupancyAnalyzer OccupancyAnalyzer.cc Occupancy/OccupancyAnalyzer/src/OccupancyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yusuf Erdogan,32 4-B20,+41227676487
//         Created:  Fri Sep 20 09:53:14 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/PhysObjectMatcher.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "TMath.h"
#include "TH2F.h"

//
// class declaration
//

class OccupancyAnalyzer : public edm::EDAnalyzer {
   public:
      explicit OccupancyAnalyzer(const edm::ParameterSet&);
      ~OccupancyAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      
      std::map<std::string, TH1*> histos_;
	  edm::Service<TFileService> theFileService;
      
};

//
// constants, enums and typedefs
//

Double_t const  kPI        = TMath::Pi();
Double_t const  kTWOPI     = 2.*kPI;

typedef math::XYZTLorentzVector LorentzVector;
typedef edm::Ref<reco::GenJetCollection> ObjectRef;

//
// static data member definitions
//


//
// constructors and destructor
//
OccupancyAnalyzer::OccupancyAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   	histos_["gen_status"] = theFileService->make<TH1F> ("gen_status", "gen_status", 10, -0.5, 9.5);
	histos_["gen_pdgid"] = theFileService->make<TH1F> ("gen_pdgid", "gen_pdgid", 12000, -6000.5, 5999.5);
	histos_["n_muons"] = theFileService->make<TH1F> ("n_muons", "n_muons", 20, -0.5, 19.5);
	histos_["phi_eta_muons"] = theFileService->make<TH2F> ("phi_eta_muons", "phi_eta_muons", 100, -kPI, kPI, 100, -1.3, 1.3);
	histos_["phi_eta_parts"] = theFileService->make<TH2F> ("phi_eta_parts", "phi_eta_parts", 100, -kPI, kPI, 100, -1.32, 1.32);
	
	histos_["const_size"] = theFileService->make<TH1F> ("const_size", "const_size", 10, -0.5, 9.5);

	histos_["n_ak5GenJets"] = theFileService->make<TH1F> ("n_ak5GenJets", "n_ak5GenJets", 20, -0.5, 19.5);
	
	histos_["e_diff"] = theFileService->make<TH1F> ("e_diff", "e_diff", 400, -2., 2.);
	histos_["e_part"] = theFileService->make<TH1F> ("e_part", "e_part", 200, 0., 2.);
	histos_["e_part_below5"] = theFileService->make<TH1F> ("e_part_below10", "e_part_below10", 200, 0., 2.);
	histos_["e_part_5to10"] = theFileService->make<TH1F> ("e_part_10to20", "e_part_10to20", 200, 0., 2.);
	histos_["e_part_10to20"] = theFileService->make<TH1F> ("e_part_10to20", "e_part_10to20", 200, 0., 2.);
	histos_["e_part_20to30"] = theFileService->make<TH1F> ("e_part_20to30", "e_part_20to30", 200, 0., 2.);
	histos_["e_part_30to40"] = theFileService->make<TH1F> ("e_part_30to40", "e_part_30to40", 200, 0., 2.);
	histos_["e_part_above40"] = theFileService->make<TH1F> ("e_part_above40", "e_part_above40", 200, 0., 2.);
	
	histos_["e_punchthr"] = theFileService->make<TH1F> ("e_punchthr", "e_punchthr", 1000, 0., 1000.);
	
	histos_["e_reco"] = theFileService->make<TH1F> ("e_reco", "e_reco", 100, 0., 100.);
	histos_["e_gen"] = theFileService->make<TH1F> ("e_gen", "e_gen", 100, 0., 100.);
	
	histos_["dR"] = theFileService->make<TH1F> ("dR", "dR", 100, 0., 0.1);
	histos_["mu_hit_id"] = theFileService->make<TH1F> ("mu_hit_id", "mu_hit_id", 100, -0.5, 99.5);
}


OccupancyAnalyzer::~OccupancyAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
OccupancyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
	
	Handle<vector<reco::GenParticle>> genparts;
	iEvent.getByLabel("genParticles", genparts);
	
	int n_mu = 0;
	int n_part = 0;
	
	for(vector<reco::GenParticle>::const_iterator gen_i = genparts->begin();  gen_i != genparts->end(); gen_i++){
		histos_["gen_status"]->Fill((*gen_i).status());
		if((*gen_i).status() != 1) continue;
		
		histos_["gen_pdgid"]->Fill((*gen_i).pdgId());
		n_part++;
		const LorentzVector part_mom = (*gen_i).p4();
		histos_["phi_eta_parts"]->Fill(part_mom.phi(),part_mom.eta());
		if(abs((*gen_i).pdgId()) == 13){
			n_mu++;
			if(part_mom.e() >= 6.){
				histos_["phi_eta_muons"]->Fill(part_mom.phi(),part_mom.eta());
			}
		}
	}
	
	Handle<vector<reco::GenJet>> ak5GenJets;
	iEvent.getByLabel("ak5GenJets", ak5GenJets);
	
	Handle<View<reco::PFJet> > ak5PFJets;
	iEvent.getByLabel("ak5PFJetsJEC", ak5PFJets);
	
	Handle<Association<reco::GenJetCollection>> genJetMatch;
	iEvent.getByLabel("JetGenJetMatch", genJetMatch);

	for (unsigned int idx = 0; idx < ak5PFJets->size(); idx++) {

		edm::RefToBase<reco::PFJet> jetRef = ak5PFJets->refAt(idx);
		reco::PFJet ajet(*(jetRef.get()));
		
		reco::GenJetRef genjetref = (*genJetMatch)[jetRef];
		
    	if (genjetref.isNonnull() && genjetref.isAvailable()) {
    		reco::GenJet genjet = *(genjetref.get());
    		
    		//double gentheta = -99.;
			//double genangle_rad = -99.;
			//gentheta = 2*atan(exp(-1*genjet.eta()));
			//genangle_rad = 1.5708 - gentheta;
			
			//double recotheta = -99.;
			//double recoangle_rad = -99.;
			//recotheta = 2*atan(exp(-1*ajet.eta()));
			//recoangle_rad = 1.5708 - recotheta;
			
			histos_["e_reco"]->Fill( ajet.energy() );
			histos_["e_gen"]->Fill( genjet.energy() );
			
			double deta = ajet.eta()-genjet.eta();
			double dphi = ajet.phi()-genjet.phi();
			while (dphi >= kPI) dphi -= kTWOPI;
   			while (dphi < -kPI) dphi += kTWOPI;
   			
			double dR = sqrt( deta*deta+dphi*dphi );
    		
    		histos_["dR"]->Fill(dR);
    		
    		histos_["e_diff"]->Fill( (ajet.energy() - genjet.energy()) / ajet.energy() );
    		histos_["e_part"]->Fill( ajet.energy() / genjet.energy() );
    		
    		if(ajet.energy() <= 5.) histos_["e_part_below5"]->Fill( ajet.energy() / genjet.energy() );
    		if(ajet.energy() > 5. && ajet.energy() <= 10.) histos_["e_part_5to10"]->Fill( ajet.energy() / genjet.energy() );
			if(ajet.energy() > 10. && ajet.energy() <= 20.) histos_["e_part_10to20"]->Fill( ajet.energy() / genjet.energy() );
			if(ajet.energy() > 20. && ajet.energy() <= 30.) histos_["e_part_20to30"]->Fill( ajet.energy() / genjet.energy() );
			if(ajet.energy() > 30. && ajet.energy() <= 40.) histos_["e_part_30to40"]->Fill( ajet.energy() / genjet.energy() );
			if(ajet.energy() > 40.) histos_["e_part_above40"]->Fill( ajet.energy() / genjet.energy() );
			
			if( (ajet.energy() / genjet.energy()) > 1.){
				histos_["e_punchthr"]->Fill(ajet.energy());	
			}
		}
	}
	
	histos_["n_ak5GenJets"]->Fill(ak5GenJets->size());
	histos_["n_muons"]->Fill(n_mu);

	Handle<vector<PSimHit>> muonSimHits;
	iEvent.getByLabel("MuonSimHits","MuonDTHits", muonSimHits);
	
	for(vector<PSimHit>::const_iterator simhit_i = muonSimHits->begin();  simhit_i != muonSimHits->end(); simhit_i++){
		histos_["mu_hit_id"]->Fill(abs(simhit_i->particleType()));
	}


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
OccupancyAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
OccupancyAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
OccupancyAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
OccupancyAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
OccupancyAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
OccupancyAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
OccupancyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OccupancyAnalyzer);
