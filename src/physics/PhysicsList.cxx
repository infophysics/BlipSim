/**
 * @file PhysicsList.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-27
 */
#include "PhysicsList.hh"

namespace Blip
{
    PhysicsList::PhysicsList()
    : G4VModularPhysicsList()
    {   
        SetVerboseLevel(1);

        new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
        new G4UnitDefinition( "mm2/g",  "mm2/g", "Surface/Mass", mm2/g);
        new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);
        
        // // Standard EM Physics
        // RegisterPhysics(new G4EmStandardPhysics());

        // // Synchroton Radiation & GN Physics
        // RegisterPhysics(new G4EmExtraPhysics());

        // // Optical Physics
        // RegisterPhysics(new G4OpticalPhysics());

        // // Decay Physics
        // RegisterPhysics(new G4DecayPhysics());
        // RegisterPhysics(new G4RadioactiveDecayPhysics());

        // // Stopping Physics
        // RegisterPhysics(new G4StoppingPhysics());

        // // Ion Physics
        // RegisterPhysics(new G4IonPhysics());

        // // Hadron Elastic scattering
        // RegisterPhysics(new G4HadronElasticPhysicsHP());

        // // Hadron Physics
        // RegisterPhysics(new G4HadronPhysicsQGSP_BERT_HP());

        // // Step limiter
        // RegisterPhysics(new G4StepLimiterPhysics());

        // Neutron Physics
        RegisterPhysics(new NeutronHPPhysics("neutronHP"));
    }

    PhysicsList::~PhysicsList()
    {
    }

    void PhysicsList::ConstructParticle()
    {
        G4BosonConstructor  pBosonConstructor;
        pBosonConstructor.ConstructParticle();

        G4LeptonConstructor pLeptonConstructor;
        pLeptonConstructor.ConstructParticle();

        G4MesonConstructor pMesonConstructor;
        pMesonConstructor.ConstructParticle();

        G4BaryonConstructor pBaryonConstructor;
        pBaryonConstructor.ConstructParticle();

        G4IonConstructor pIonConstructor;
        pIonConstructor.ConstructParticle();

        G4ShortLivedConstructor pShortLivedConstructor;
        pShortLivedConstructor.ConstructParticle();  
    }

    void PhysicsList::SetCuts()
    {
        SetCutValue(0 * mm, "proton");
        // SetCutValue(1 * um, "gamma");
        // SetCutValue(1 * um, "e-");
        // SetCutValue(1 * um, "e+");
        SetCutValue(0 * mm, "neutron");
    }

    void PhysicsList::PrintPhysicsLists()
    {
        // G4cout << "Enabled Physics Lists:" << G4endl;
        // for(size_t ii = 0; ii < mPhysicsLists.size(); ii++)
        // {
        //     G4cout << "\t[" << ii << "]: " << mPhysicsLists[ii] << G4endl;
        // }
    }
}