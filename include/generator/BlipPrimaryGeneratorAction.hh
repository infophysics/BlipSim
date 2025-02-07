/**
 * @file BlipPrimaryGeneratorAction.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-27
 */
#pragma once
#include <memory>
#include <random>
#include <array>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "EventManager.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"

#include "Core.hh"
#include "EventManager.hh"

#include "yaml-cpp/yaml.h"

namespace Blip
{
    class BlipPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
    {   
    public:
        BlipPrimaryGeneratorAction();
        ~BlipPrimaryGeneratorAction();

        virtual void GeneratePrimaries(G4Event* event);

        BlipPrimaryGeneratorAction(YAML::Node config);
        YAML::Node Config() const { return mConfig; }

    private:
        G4ParticleGun* mParticleGun;

        G4String mParticleName;
        G4double mParticleMomentum;
        G4double mParticleEnergy;
        G4ThreeVector mParticlePosition;
        G4ThreeVector mParticleMomentumDirection;

        TRandom3* mTRandom3 = {0};

        YAML::Node mConfig;
    };
}