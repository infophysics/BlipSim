/**
 * @file BlipPrimaryGeneratorAction.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12-13
 */
#include "BlipPrimaryGeneratorAction.hh"

namespace Blip
{
    BlipPrimaryGeneratorAction::BlipPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
    {
        mParticleGun = new G4ParticleGun();
        G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        mParticleGun->SetNumberOfParticles(1);
        mParticleGun->SetParticleEnergy(57*keV);
        mParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
        mParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
        mParticleGun->SetParticleDefinition(particle);
    }

    BlipPrimaryGeneratorAction::~BlipPrimaryGeneratorAction()
    {
    }

    BlipPrimaryGeneratorAction::BlipPrimaryGeneratorAction(YAML::Node config)
    : G4VUserPrimaryGeneratorAction()
    , mConfig(config)
    {
        mTRandom3 = new TRandom3();
        mParticleGun = new G4ParticleGun();
        G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        mParticleGun->SetNumberOfParticles(1);
        mParticleGun->SetParticleEnergy(57*keV);
        mParticleGun->SetParticlePosition(G4ThreeVector(0.,0., mTZeroLocation));
        mParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
        mParticleGun->SetParticleDefinition(particle);
    }

    void BlipPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
    {
        mParticleGun->SetNumberOfParticles(1);
        G4double BeamEnergy = SampleBeamEnergy();
        G4double deltaTOF = SampleTOF(BeamEnergy);
        G4ThreeVector StartPosition = SampleBeamProfile(mTZeroLocation);
        G4ThreeVector StartMomentum = SampleBeamMomentum(StartPosition);
        mParticleGun->SetParticleTime(deltaTOF);
        mParticleGun->SetParticlePosition(StartPosition);
        mParticleGun->SetParticleMomentumDirection(StartMomentum);
        mParticleGun->SetParticleEnergy(BeamEnergy);
        mParticleGun->GeneratePrimaryVertex(event);
    }
}