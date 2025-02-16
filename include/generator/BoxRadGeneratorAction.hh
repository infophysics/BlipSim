#pragma once
#include <memory>
#include <random>
#include <array>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "EventManager.hh"
#include "G4ThreeVector.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"

#include "SpectrumReader.hh" // -- Reads radiological energy spectrum

#include "Core.hh"
#include "EventManager.hh"

#include "yaml-cpp/yaml.h"

namespace Blip
{
    class BoxRadGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
    public:
        BoxRadGeneratorAction();
        BoxRadGeneratorAction(G4LogicalVolume* volume, G4String radName="Argon_39");
        BoxRadGeneratorAction(YAML::Node config); // -- constructor that takes the radName and the logical volume from the config
        
        ~BoxRadGeneratorAction();

        virtual void GeneratePrimaries(G4Event* event) override;

        YAML::Node Config() const {return mConfig; }

        G4ThreeVector GetRandomPositionInVolume(G4LogicalVolume* volume);

        private:
            G4ParticleGun* mParticleGun;
            SpectrumReader* mSpectrumReader; // -- Reader for radiological
            G4LogicalVolume* mVolume; // -- Logical volume for the rads
            G4String mRadName; // -- name of radiological (e.g. Ar39, Ar42, etc)
            G4String mSpectrumPath; // -- path of radioIsotope file
            G4double mFixedEnergy; // -- fixed energy for now

            G4String mParticleName;
            G4double mParticleMomentum;
            G4double mParticleEnergy;
            G4ThreeVector mParticlePosition;
            G4ThreeVector mParticleMomentumDirection;

            TRandom3* mTRandom3 = {0};

            YAML::Node mConfig;

            int mNBetas;
    };
}