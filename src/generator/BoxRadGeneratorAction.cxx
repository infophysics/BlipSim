/**
 * @file BoxRadGeneratorAction.cxx
 * @author David Rivera [rivera@lanl.gov]
 * @brief Radiological generator action implementation
 * @version 0.1
 * @date 2025-02-12
 */

 #include "BoxRadGeneratorAction.hh"
 #include "BlipDetectorConstruction.hh"

 #include "G4RunManager.hh"
 #include "G4VUserDetectorConstruction.hh"
 #include "G4LogicalVolume.hh"
 #include "G4PhysicalConstants.hh" 
 #include "G4SystemOfUnits.hh"
 #include "G4Electron.hh"
 #include "G4Box.hh"
 #include "G4ThreeVector.hh"
 #include "G4PrimaryVertex.hh"

 #include <fstream>
 #include <sys/stat.h>
 #include <string>


 namespace Blip
 {
    // -- DUMMY function for now
    BoxRadGeneratorAction::BoxRadGeneratorAction()
    {
        const G4VUserDetectorConstruction* detectorConstruction = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
        const BlipDetectorConstruction* LArBox = dynamic_cast<const BlipDetectorConstruction*>(detectorConstruction);
        mLogicalVolume = LArBox->GetArgonLogicalVolume(); // -- need to get this from the detector construction
        mRadName = "Argon_39";
        mTRandom3 = new TRandom3();
        mParticleGun = new G4ParticleGun(1); // -- single rads at a time... I guess
        mParticleGun->SetParticleDefinition(G4Electron::Electron()); // -- the beta
        mParticleGun->SetParticleEnergy(500 * keV); // -- default energy
        mParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }    

    BoxRadGeneratorAction::BoxRadGeneratorAction(G4LogicalVolume* volume, G4String radName)
    : mLogicalVolume(volume), mRadName(radName), mSpectrumReader(new SpectrumReader(radName, "./rads/","Argon_39.root"))
    {
        const G4VUserDetectorConstruction* detectorConstruction = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
        const BlipDetectorConstruction* LArBox = dynamic_cast<const BlipDetectorConstruction*>(detectorConstruction);
        mLogicalVolume = LArBox->GetArgonLogicalVolume(); // -- need to get this from the detector construction
        mTRandom3 = new TRandom3();
        mParticleGun = new G4ParticleGun(1); // -- single rads at a time... I guess
        mParticleGun->SetParticleDefinition(G4Electron::Electron()); // -- the beta
        mParticleGun->SetParticleEnergy(500 * keV); // -- default energy
        mParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }

    BoxRadGeneratorAction::~BoxRadGeneratorAction()
    {
        delete mParticleGun;
        delete mSpectrumReader;
    }

    // -- DUMMY function for now
    BoxRadGeneratorAction::BoxRadGeneratorAction(YAML::Node config)
    : mConfig(config)
    {

        // -- read in and check the config. TODO: add more thorough checks..
        if (mConfig["generator"]["radiological"])           { mRadName = mConfig["generator"]["radiological"].as<std::string>() ; }
        if (mConfig["generator"]["rateInBqPerCC"])          { rateInBq = mConfig["generator"]["rateInBqPerCC"].as<G4double>(); } // -- decays per sec per cm^3
        if (mConfig["generator"]["nDecays"])                { mNBetas = mConfig["generator"]["nDecays"].as<int>() ; }
        if (mConfig["generator"]["spectrumPath"])           { mSpectrumPath = mConfig["generator"]["spectrumPath"].as<std::string>() ; }
        if (mConfig["generator"]["energy"])                 { mFixedEnergy = mConfig["generator"]["energy"].as<G4double>() * keV; }

        G4bool fileExists = false;
        std::string fullPath = mSpectrumPath + mRadName + ".root";
        std::ifstream f(fullPath.c_str());
        fileExists = ( f.good() ? true : false );

        if (fileExists)
        {
            std::cout << "Using file: " << fullPath << " to sample rads" << std::endl;
            mSpectrumReader = new SpectrumReader(mRadName, mSpectrumPath, mRadName+".root");
            //<--auto const & mEnergies = mSpectrumRead->GetRandomEnergy(mNBetas);
        } else {
            std::cout << "Cannot find file : " << fullPath << ". Will generate fixed energy betas" << std::endl;
            mRadName = "Argon_39";
            mFixedEnergy = 500 * keV;
        }


        const G4VUserDetectorConstruction* detectorConstruction = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
        const BlipDetectorConstruction* LArBox = dynamic_cast<const BlipDetectorConstruction*>(detectorConstruction);
        mLogicalVolume = LArBox->GetArgonLogicalVolume(); // -- need to get this from the detector construction
        
        // -- Get the box corresponding to the logical volume
        mBox = dynamic_cast<G4Box*>(mLogicalVolume->GetSolid());
        if (mBox) { mVolume = mBox->GetCubicVolume(); } // -- volume in mm^3
        std::cout << "Argon box volume = " << mVolume << std::endl;

        mTRandom3 = new TRandom3();

        // -- Calculate the number of decays to generate
        G4double mVolume_cm3 = mVolume / 1000.; // going from mm^3 to cm^3
        G4double deltaT_s = 3.E-3; // -- 15 ms in seconds to match the units of Bq
        //double meanNumberOfDecays = (rateInBq * mVolume_cm3 * deltaT_s)/(1.E6);
        double meanNumberOfDecays = (rateInBq * mVolume_cm3 * deltaT_s);

        // -- sample from a Poisson distribution
        mNBetas = mTRandom3->Poisson(meanNumberOfDecays);
        std::cout << "Generating " << mNBetas << " beta decays" << std::endl;

        // -- test
        //mNBetas = 100;

        // -- sample the decays
        //bool test = mSpectrumReader->SampleTheRadiological(mNBetas);
        bool test = mSpectrumReader->GenNDecays(mRadName, mNBetas);
        if (test) { std::cout << "Successfully sampled the spectrum " << mNBetas << " times!\n"; }

        //mParticleGun = new G4ParticleGun(1); // -- Multiple decays
        //mParticleGun->SetParticleDefinition(G4Electron::Electron()); // -- the beta
        //mParticleGun->SetParticleEnergy(mFixedEnergy); // -- default energy
        //mParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }    

    void BoxRadGeneratorAction::AddDecayParticle(G4Event* event, const int pdg, const G4ThreeVector& pos, const double time, const double energy, const G4ThreeVector& mom)
    {
        auto vertex = new G4PrimaryVertex{pos, time};
        vertex->SetPrimary(new G4PrimaryParticle{pdg, mom.x(), mom.y(), mom.z(), energy});
        event->AddPrimaryVertex(vertex);
    }

    void BoxRadGeneratorAction::GeneratePrimaries(G4Event* event)
    {
        // loop over the beta energies
        std::vector<double> theKEs = mSpectrumReader->getKineticEnergies();

        size_t NBetas = theKEs.size();

        // -- get the mass
        G4double mass = G4Electron::Electron()->GetPDGMass();
        std::cout << "Electron mass = " << mass << std::endl;
        std::cout << "Producing " << NBetas << " decays" << std::endl;

        for (size_t i = 0; i<NBetas; i++)
        {
            // -- sample the energy from the ROOT spectrum
            //<--G4double energy = mSpectrumReader->GetRandomEnergy(); // -- TODO: get the spectrum files and fix this
            G4double T = theKEs.at(i);
            G4double energy = T + mass;
            G4double p = energy * energy - mass * mass;
            if (p>=0.) { p=TMath::Sqrt(p); }
            else { p=0.; }

            TLorentzVector p4 = GetDirection(mass, p);
            G4ThreeVector pVec{p4.Px(), p4.Py(), p4.Pz()};

            // -- sample random position from within the logical volume
            G4ThreeVector position = GetRandomPositionInVolume(mLogicalVolume);

            std::cout << "Generating beta with decay energy of : " << T << ", mass of " << mass << ", and total energy of: " << energy << std::endl;
            AddDecayParticle(event, 11, position, 0.0, energy, pVec);
        }
    }    

    // -- Box sampling
    G4ThreeVector BoxRadGeneratorAction::GetRandomPositionInVolume(G4LogicalVolume* volume)
    {
        G4Box* box = dynamic_cast<G4Box*>(volume->GetSolid());
        double box_XhL=0.0, box_YhL=0.0, box_ZhL=0.0;

        G4ThreeVector pos(0.,0.,0.);
        if (box) {
            box_XhL = box->GetXHalfLength();
            box_YhL = box->GetYHalfLength();
            box_ZhL = box->GetZHalfLength();
            //std::cout << "Box half lengths [mm] (x,y,z) : " 
            //                    << box_XhL << ", " 
            //                    << box_YhL << ", " 
            //                    << box_ZhL << std::endl;

            
            pos.setX(mTRandom3->Uniform(-1.0*box_XhL, box_XhL));
            pos.setY(mTRandom3->Uniform(-1.0*box_YhL, box_YhL));
            pos.setZ(mTRandom3->Uniform(-1.0*box_ZhL, box_ZhL));
        }
        return pos;
    }
    
    TLorentzVector BoxRadGeneratorAction::GetDirection(double& p, double& m)
    {
        // -- make it isotropic
        double cosTheta = (2. * mTRandom3->Uniform(0., 1.) - 1.);
        if (cosTheta < -1.0) cosTheta = -1.;
        if (cosTheta > 1.0) cosTheta = 1.;
        double const sinTheta = sqrt(1. - cosTheta * cosTheta);
        double const phi = 2. * TMath::Pi() * mTRandom3->Uniform(0., 1.);
        
        TLorentzVector theDirection(p*sinTheta*std::cos(phi), p*sinTheta*std::sin(phi), p*cosTheta, std::sqrt(p*p + m*m));

        return theDirection;
    }
 }