#include <iomanip>

// -- artg4tk includes
#include "MyG4HadronPhysicsQGSP_BERT_ArHP.hh"
#include "MyG4NeutronPHPBuilder.hh"

//#include "G4HadronPhysicsQGSP_BERT_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4NeutronBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
#include "G4BertiniNeutronBuilder.hh"
//#include "G4NeutronPHPBuilder.hh"

#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

#include "G4HadronicParameters.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//

namespace Blip
{
  G4_DECLARE_PHYSCONSTR_FACTORY(MyG4HadronPhysicsQGSP_BERT_ArHP);

  MyG4HadronPhysicsQGSP_BERT_ArHP::MyG4HadronPhysicsQGSP_BERT_ArHP(G4int)
      :  MyG4HadronPhysicsQGSP_BERT_ArHP("hInelastic MyQGSP_BERT_ArHP")
  {}

  MyG4HadronPhysicsQGSP_BERT_ArHP::MyG4HadronPhysicsQGSP_BERT_ArHP(const G4String& name, G4bool /*quasiElastic */ )
      :  G4HadronPhysicsQGSP_BERT(name)
  {
      minBERT_neutron = 19.9*MeV;
  }

  void MyG4HadronPhysicsQGSP_BERT_ArHP::Neutron()
  {
    auto neu = new G4NeutronBuilder( true ); // Fission on
    AddBuilder(neu);
    auto qgs = new G4QGSPNeutronBuilder(QuasiElasticQGS);
    AddBuilder(qgs);
    qgs->SetMinEnergy(minQGSP_neutron);
    neu->RegisterMe(qgs);
    auto ftf = new G4FTFPNeutronBuilder(QuasiElasticFTF);
    AddBuilder(ftf);
    ftf->SetMinEnergy(minFTFP_neutron);
    ftf->SetMaxEnergy(maxFTFP_neutron);
    neu->RegisterMe(ftf);
    auto bert = new G4BertiniNeutronBuilder;
    AddBuilder(bert);
    bert->SetMinEnergy(minBERT_neutron);
    bert->SetMaxEnergy(maxBERT_neutron);
    neu->RegisterMe(bert);
    auto hp = new MyG4NeutronPHPBuilder;
    AddBuilder(hp);
    neu->RegisterMe(hp);
    neu->Build();
  }

  void MyG4HadronPhysicsQGSP_BERT_ArHP::ExtraConfiguration()
  {
    // --- Neutrons ---
    const G4ParticleDefinition* neutron = G4Neutron::Neutron();
    G4HadronicProcess* capture = G4PhysListUtil::FindCaptureProcess(neutron);
    if (capture) {
      G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture();
      theNeutronRadCapture->SetMinEnergy( minBERT_neutron );
      capture->RegisterMe( theNeutronRadCapture );
    }
    G4HadronicProcess* fission = G4PhysListUtil::FindFissionProcess(neutron);
    if (fission) {
      G4LFission* theNeutronLEPFission = new G4LFission();
      theNeutronLEPFission->SetMinEnergy( minBERT_neutron );
      theNeutronLEPFission->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
      fission->RegisterMe( theNeutronLEPFission );
    }
  }
}
