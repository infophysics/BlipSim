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
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "Randomize.hh"
#include "G4EventManager.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleGun.hh"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

#include "TLorentzVector.h"

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
    struct particle {
        int pdgID;
        int genID;
        G4ThreeVector position;
        TLorentzVector momentum4Vector;
    };

    class BoxRadGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
    public:
        BoxRadGeneratorAction();
        BoxRadGeneratorAction(G4LogicalVolume* volume, G4String radName="Argon_39");
        BoxRadGeneratorAction(YAML::Node config); // -- constructor that takes the radName and the logical volume from the config
        
        ~BoxRadGeneratorAction();

        virtual void GeneratePrimaries(G4Event* event) override;
        void AddDecayParticle(G4Event* event, const int pdg, const G4ThreeVector& pos, const double time, const double energy, const G4ThreeVector& momentum);

        YAML::Node Config() const {return mConfig; }

        G4ThreeVector GetRandomPositionInVolume(G4LogicalVolume* volume);
        TLorentzVector GetDirection(double& p, double& m);

        private:
            G4ParticleGun* mParticleGun;
            SpectrumReader* mSpectrumReader;    // -- Reader for radiological
            G4LogicalVolume* mLogicalVolume;    // -- Logical volume for the rads
            G4Box* mBox;                        // -- box of argon
            G4double mVolume;                   // -- volume of argon box in units of mm^3
            std::string mRadName;               // -- name of radiological (e.g. Argon_39, Argon_42, etc)
            G4double mRateInBq;                 // -- Bequerel
            std::string mSpectrumPath;          // -- path of radioIsotope file
            G4double mFixedEnergy;              // -- fixed energy for now
            G4double mStartT;                   // -- start time in ns
            G4double mEndT;                     // -- end time in ns

            G4String mParticleName;
            G4double mParticleMomentum;
            G4double mParticleEnergy;
            G4ThreeVector mParticlePosition;
            G4ThreeVector mParticleMomentumDirection;

            // -- Random Engine
            CLHEP::HepRandomEngine* mEngine = {0};

            //TRandom3* mTRandom3 = {0};
            TRandom3* mTRandom3;
            //RandPoissonQ* mPoissonQ;

            YAML::Node mConfig;

            long mNBetas;
    };
}