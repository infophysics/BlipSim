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
 #include "G4Neutron.hh"
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
    : mConfig(config),
    mRadName("Argon_39"),
    mRateInBq(0.00141),
    mFixedEnergy(500*keV),
    mStartT(0.),
    mEndT(3.E6)
    {
        CLHEP::RandFlat mFlat(*mEngine);
        CLHEP::RandPoisson mPoisson(*mEngine);
        // -- read in and check the config. TODO: add more thorough checks..
        if (mConfig["generator"]["radiological"])           { mRadName = mConfig["generator"]["radiological"].as<std::string>() ; }
        if (mConfig["generator"]["rateInBqPerCC"])          { mRateInBq = mConfig["generator"]["rateInBqPerCC"].as<G4double>(); } // -- decays per sec per cm^3
        if (mConfig["generator"]["spectrumPath"])           { mSpectrumPath = mConfig["generator"]["spectrumPath"].as<std::string>() ; }
        if (mConfig["generator"]["energy"])                 { mFixedEnergy = mConfig["generator"]["energy"].as<G4double>() * keV; }
        if (mConfig["generator"]["start_time_ns"])          { mStartT = mConfig["generator"]["start_time_ns"].as<G4double>() * ns; }
        if (mConfig["generator"]["end_time_ns"])            { mEndT = mConfig["generator"]["end_time_ns"].as<G4double>() * ns; }

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

        mTRandom3 = new TRandom3(1094);

        // -- Calculate the number of decays to generate
        G4double mVolume_cm3 = mVolume / 1000.; // going from mm^3 to cm^3
        //G4double deltaT_s = 3.E-3; // -- 3 ms in seconds to match the units of Bq
        G4double deltaT_s = (mEndT - mStartT)/(1.E9); // -- seconds to match units of Bq
        //double meanNumberOfDecays = (mRateInBq * mVolume_cm3 * deltaT_s)/(1.E6);
        double meanNumberOfDecays = (mRateInBq * mVolume_cm3 * deltaT_s);

        bool overrideNDecays = false;
        if (mConfig["generator"]["overrideNDecays"]) { overrideNDecays = mConfig["generator"]["overrideNDecays"].as<bool>(); }

        // -- some fixed number of decays to make it simpler
        if (overrideNDecays) {
            mNBetas = mConfig["generator"]["nDecays"].as<int>() ; 
        } else {
            // -- sample from a Poisson distribution
            //mPoissonQ = new RandPoissonQ();
            //mNBetas = mPoissonQ->Poisson(meanNumberOfDecays);
            //mNBetas = mPoisson.fire(meanNumberOfDecays);
            mNBetas = mTRandom3->Poisson(meanNumberOfDecays);
        }
        std::cout << "Generating " << mNBetas << " beta decays" << ", (mean number should be " << meanNumberOfDecays << ")" << std::endl;

        // -- sample the decays
        //bool test = mSpectrumReader->SampleTheRadiological(mNBetas);
        bool test = mSpectrumReader->GenNDecays(mRadName, mNBetas);
        if (test) { std::cout << "Successfully sampled the spectrum " << mNBetas << " times!\n"; }
    }    

    void BoxRadGeneratorAction::AddDecayParticle(G4Event* event, const int pdg, const G4ThreeVector& pos, const double time, const double energy, const G4ThreeVector& mom)
    {
        auto vertex = new G4PrimaryVertex{pos, time};
        vertex->SetPrimary(new G4PrimaryParticle{pdg, mom.x(), mom.y(), mom.z(), energy});
        event->AddPrimaryVertex(vertex);
    }

    void BoxRadGeneratorAction::GeneratePrimaries(G4Event* event)
    {
        // -- Add a thermal neutron in here for now. The neutron will always be the first particle
        G4double thermalEnergy = (1/40.)/(1.0E6); // -- 1/40th of an eV in MeV
        G4ThreeVector neutronPos = GetRandomPositionInVolume(mLogicalVolume);
        G4double neutronMass = G4Neutron::Neutron()->GetPDGMass();
        G4double neutronEnergy = thermalEnergy + neutronMass; // -- in MeV
        G4double neutronP = neutronEnergy * neutronEnergy - neutronMass * neutronMass;
        if (neutronP>=0.) { neutronP=TMath::Sqrt(neutronP); }
        else { neutronP=0.; }

        TLorentzVector neutronP4 = GetDirection(neutronMass, neutronP);
        G4ThreeVector neutronPVec{neutronP4.Px(), neutronP4.Py(), neutronP4.Pz()};

        // -- It's the same procedure for adding a neutrons as the decay particles
        AddDecayParticle(event, 2112, neutronPos, 0.0, neutronEnergy, neutronPVec);

        // loop over the beta energies
        std::vector<double> theKEs = mSpectrumReader->getKineticEnergies();

        size_t NBetas = theKEs.size();

        // -- get the mass
        G4double mass = G4Electron::Electron()->GetPDGMass();
        G4cout << "Electron mass = " << mass << "\n";
        G4cout << "Producing " << NBetas << " decays" << "\n";

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

            // -- sample random time between the start and end time of the simulation
            G4double t_i = mTRandom3->Uniform(mStartT, mEndT);

            //<--std::cout << "Generating beta with decay energy of : " << T << ", mass of " << mass << ", and total energy of: " << energy << std::endl;
            AddDecayParticle(event, 11, position, t_i, energy, pVec);
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