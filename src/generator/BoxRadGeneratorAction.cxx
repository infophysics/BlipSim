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
        mVolume = LArBox->GetArgonLogicalVolume(); // -- need to get this from the detector construction
        mRadName = "Argon_39";
        mTRandom3 = new TRandom3();
        mParticleGun = new G4ParticleGun(1); // -- single rads at a time... I guess
        mParticleGun->SetParticleDefinition(G4Electron::Electron()); // -- the beta
        mParticleGun->SetParticleEnergy(500 * keV); // -- default energy
        mParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }    

    BoxRadGeneratorAction::BoxRadGeneratorAction(G4LogicalVolume* volume, G4String radName)
    : mVolume(volume), mRadName(radName), mSpectrumReader(new SpectrumReader("rads/Argon_39.root"))
    {
        const G4VUserDetectorConstruction* detectorConstruction = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
        const BlipDetectorConstruction* LArBox = dynamic_cast<const BlipDetectorConstruction*>(detectorConstruction);
        mVolume = LArBox->GetArgonLogicalVolume(); // -- need to get this from the detector construction
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
        if (mConfig["generator"]["radiological"])           { mRadName = mConfig["generator"]["number"].as<G4String>() ; }
        if (mConfig["generator"]["nDecays"])                { mNBetas = mConfig["generator"]["number"].as<int>() ; }
        if (mConfig["generator"]["spectrumPath"])           { mSpectrumPath = mConfig["generator"]["filePath"].as<G4String>() ; }
        if (mConfig["generator"]["energy"])                 { mFixedEnergy = mConfig["generator"]["energy"].as<G4double>() ; }

        G4bool fileExists = false;
        std::string fullPath = mSpectrumPath + mRadName + ".root";
        ifstream f(fullPath.c_str());
        fileExists = ( f.good() ? true : false );

        if (fileExisist)
        {
            std::cout << "Using file: " << fullPath << " to sample rads" << std::endl;
            mSpectrumRead = new SpectrumReader(fullPath);
            auto const & mEnergies = mSpectrumRead->GetRandomEnergy(mNBetas);
        } else {
            std::cout << "Cannot find file : " << fullPath << ". Will generate fixed energy betas" << std::endl;
            mRadName = "Argon_39";
            mFixedEnergy = 500 * keV;
        }


        const G4VUserDetectorConstruction* detectorConstruction = G4RunManager::GetRunManager()->GetUserDetectorConstruction();
        const BlipDetectorConstruction* LArBox = dynamic_cast<const BlipDetectorConstruction*>(detectorConstruction);
        mVolume = LArBox->GetArgonLogicalVolume(); // -- need to get this from the detector construction

        mTRandom3 = new TRandom3();

        mParticleGun = new G4ParticleGun(mNBetas); // -- Multiple decays
        mParticleGun->SetParticleDefinition(G4Electron::Electron()); // -- the beta
        mParticleGun->SetParticleEnergy(mFixedEnergy); // -- default energy
        mParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }    

    void BoxRadGeneratorAction::GeneratePrimaries(G4Event* event)
    {
        // -- sample the energy from the ROOT spectrum
        //<--G4double energy = mSpectrumReader->GetRandomEnergy(); // -- TODO: get the spectrum files and fix this
        G4double energy = 500 * keV; 
        mParticleGun->SetParticleEnergy(energy);

        // -- sample random position from within the logical volume
        G4ThreeVector position = GetRandomPositionInVolume(mVolume);
        std::cout << "Volume name for radiological: " << mVolume->GetName() << std::endl;
        mParticleGun->SetParticlePosition(position);

        // -- Generate the radiological decay products
        mParticleGun->GeneratePrimaryVertex(event);

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
 }