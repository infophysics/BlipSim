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
        RegisterPhysics(new NeutronHPPhysics("neutronHP"));
    }

    PhysicsList::PhysicsList(YAML::Node config)
    : mConfig(config)
    {
        SetVerboseLevel(1);
        new G4UnitDefinition( "millielectronVolt", "meV", "Energy", 1.e-3*eV);   
        new G4UnitDefinition( "mm2/g",  "mm2/g", "Surface/Mass", mm2/g);
        new G4UnitDefinition( "um2/mg", "um2/mg","Surface/Mass", um*um/mg);
        if(mConfig["physics"]["register_em_standard"]) {
            if(mConfig["physics"]["register_em_standard"].as<G4bool>()) {
                RegisterPhysics(new G4EmStandardPhysics());
            }
        }
        if(mConfig["physics"]["register_em_extra"]) {
            if(mConfig["physics"]["register_em_extra"].as<G4bool>()) {
                RegisterPhysics(new G4EmExtraPhysics());
            }
        }
        if(mConfig["physics"]["register_optical"]) {
            if(mConfig["physics"]["register_optical"].as<G4bool>()) {
                RegisterPhysics(new G4OpticalPhysics());
            }
        }
        if(mConfig["physics"]["register_decay"]) {
            if(mConfig["physics"]["register_decay"].as<G4bool>()) {
                RegisterPhysics(new G4DecayPhysics());
            }
        }
        if(mConfig["physics"]["register_radioactive_decay"]) {
            if(mConfig["physics"]["register_radioactive_decay"].as<G4bool>()) {
                RegisterPhysics(new G4RadioactiveDecayPhysics());
            }
        }
        if(mConfig["physics"]["register_stopping"]) {
            if(mConfig["physics"]["register_stopping"].as<G4bool>()) {
                RegisterPhysics(new G4StoppingPhysics());
            }
        }
        if(mConfig["physics"]["register_ion"]) {
            if(mConfig["physics"]["register_ion"].as<G4bool>()) {
                RegisterPhysics(new G4IonPhysics());
            }
        }
        if(mConfig["physics"]["register_hadron_elastic_hp"]) {
            if(mConfig["physics"]["register_hadron_elastic_hp"].as<G4bool>()) {
                RegisterPhysics(new G4HadronElasticPhysicsHP());
            }
        }
        if(mConfig["physics"]["register_hadron_qgsp_bert_hp"]) {
            if(mConfig["physics"]["register_hadron_qgsp_bert_hp"].as<G4bool>()) {
                RegisterPhysics(new G4HadronPhysicsQGSP_BERT_HP());
            }
        }
        // if(mConfig["physics"]["register_hadron_inclxx"]) {
        //     if(mConfig["physics"]["register_hadron_inclxx"].as<G4bool>()) {
        //         RegisterPhysics(new G4HadronPhysicsINCLXX());
        //     }
        // }
        if(mConfig["physics"]["register_step_limiter"]) {
            if(mConfig["physics"]["register_step_limiter"].as<G4bool>()) {
                RegisterPhysics(new G4StepLimiterPhysics());
            }
        }
        if(mConfig["physics"]["register_fast_simulation"]) {
            if(mConfig["physics"]["register_fast_simulation"].as<G4bool>()) {
                RegisterPhysics(new G4FastSimulationPhysics("fastSimPhys"));
            }
        }
        // if(mConfig["physics"]["register_lar_scintillation"]) {
        //     if(mConfig["physics"]["register_lar_scintillation"].as<G4bool>()) {
        //         G4OpticalPhoton::Definition();
        //         ThermalElectron::Definition();
        //         ReplacePhysics(new LArScintillationPhysics());
        //     }
        // }  
        if(mConfig["physics"]["register_neutron_hp"]) {
            if(mConfig["physics"]["register_neutron_hp"].as<G4bool>()) {  
                RegisterPhysics(new NeutronHPPhysics("neutronHP"));
            }
        }
        if(mConfig["physics"]["neutron"]["skip_missing_isotopes"]) {                                                                                
            G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes(
                mConfig["physics"]["neutron"]["skip_missing_isotopes"].as<G4bool>()
            );
        }
        if(mConfig["physics"]["neutron"]["do_not_adjust_final_state"]) {                                                                                  
            G4ParticleHPManager::GetInstance()->SetDoNotAdjustFinalState(
                mConfig["physics"]["neutron"]["do_not_adjust_final_state"].as<G4bool>()
            );
        }
        if(mConfig["physics"]["neutron"]["use_only_photo_evaporation"]) { 
            G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation(
                mConfig["physics"]["neutron"]["use_only_photo_evaporation"].as<G4bool>()
            );
        }
        if(mConfig["physics"]["neutron"]["neglect_doppler"]) { 
            G4ParticleHPManager::GetInstance()->SetNeglectDoppler(
                mConfig["physics"]["neutron"]["neglect_doppler"].as<G4bool>()
            );
        }
        if(mConfig["physics"]["neutron"]["produce_fission_fragments"]) { 
            G4ParticleHPManager::GetInstance()->SetProduceFissionFragments(
                mConfig["physics"]["neutron"]["produce_fission_fragments"].as<G4bool>()
            );
        }
        if(mConfig["physics"]["neutron"]["use_wendt_fission_model"]) { 
            G4ParticleHPManager::GetInstance()->SetUseWendtFissionModel(
                mConfig["physics"]["neutron"]["use_wendt_fission_model"].as<G4bool>()
            );
        }
        if(mConfig["physics"]["neutron"]["use_nresp71_model"]) { 
            G4ParticleHPManager::GetInstance()->SetUseNRESP71Model(
                mConfig["physics"]["neutron"]["use_nresp71_model"].as<G4bool>()
            );
        }
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
        // SetCutValue(0 * mm, "proton");
        // SetCutValue(1 * um, "gamma");
        // SetCutValue(1 * um, "e-");
        // SetCutValue(1 * um, "e+");
        SetCutValue(0 * mm, "neutron");
    }
}