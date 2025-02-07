/**
 * @file NeutronHPPhysics.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-27
 */
#pragma once
#include <iostream>
#include <fstream>
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticProcess.hh"

#ifdef BLIP_GEANT_10
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#else
#include "G4NeutronFissionProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#endif

#include "G4HadronInelasticProcess.hh"

#include "G4ParticleHPThermalScatteringData.hh"
#include "G4ParticleHPThermalScattering.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPFission.hh"

#include "G4NeutronBuilder.hh"
// #include "G4LEPNeutronBuilder.hh"
#include "G4QGSPNeutronBuilder.hh"
// #include "G4BertiniNeutronBuilder.hh"
#include "G4INCLXXNeutronBuilder.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "Argon.hh"

namespace Blip
{
    class NeutronHPPhysics : public G4VPhysicsConstructor
    {
    public:
        NeutronHPPhysics(const G4String& name="neutron");
        virtual ~NeutronHPPhysics();

        virtual void ConstructParticle();
        virtual void ConstructProcess();
        
        void SetThermalPhysics(G4bool flag) {mThermal = flag;};  
        
    private:
        G4bool mThermal = {true};
        // void CreateModels();
        // G4NeutronBuilder * theNeutrons;
        // // G4LEPNeutronBuilder * theLEPNeutron;
        // G4QGSPNeutronBuilder * theQGSPNeutron;
        // // G4BertiniNeutronBuilder * theBertiniNeutron;
        // G4INCLXXNeutronBuilder * theInclAblaNeutron;
        
    };
}