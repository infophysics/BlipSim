/**
 * @file    main.cc
 * @brief   The main function for the neutron simulation
 * @ingroup main
 * @author  Nicholas Carrara [nmcarrara@ucdavis.edu],
**/
#include <iostream>
#include <vector>
#include <string>

// Geant4 includes 
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4PhysListFactory.hh"
// Neutron HP
#include "G4ParticleHPManager.hh"
#include "G4Types.hh"
#include "G4GDMLParser.hh"

// #include "QGSP_INCLXX_HP.hh"
// #include "G4HadronPhysicsINCLXX.hh"

#include "Analysis.hh"
#include "BlipActionInitialization.hh"
#include "BlipDetectorConstruction.hh"
#include "Core.hh"
#include "EventManager.hh"
#include "PhysicsList.hh"

int main(int argc, char** argv)
{
    // load in the config file
    G4String ConfigFile = "";
    if(argc == 2) {
        ConfigFile = argv[1];
    }
    else if(argc == 1) {
        ConfigFile = "blip_sim_v1.yaml";
    }
    else 
    {
        std::cout << "ERROR! This program takes only one";
        std::cout << " command line argument for now.";
        std::cout << "  Must pass in a config file!" << std::endl;
        exit(0);
    }
    auto Config = Blip::ConfigParser(ConfigFile);
    Blip::EventManager::GetEventManager()->SetConfig(Config.GetConfig());

    // choose the Random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    
    // set the seeds
    if(Config.GetConfig()["manager"]["seed"])
    {
        G4Random::setTheSeed(Config.GetConfig()["manager"]["seed"].as<G4int>());
    }

    // create a user session
    G4UIExecutive* UIExecutive = 0;

#ifdef G4MULTITHREADED
    G4MTRunManager* RunManager = new G4MTRunManager();
    RunManager->SetNumberOfThreads(
        Blip::EventManager::GetEventManager()->NumberOfThreads()
    );
#else
    G4RunManager* RunManager = new G4RunManager;
#endif
    
    // apply the detector, physics list and initialization
    RunManager->SetUserInitialization(
        new Blip::BlipDetectorConstruction(Config.GetConfig())
    );
    RunManager->SetUserInitialization(
        new Blip::PhysicsList()
    );
    RunManager->SetUserInitialization(
        new Blip::BlipActionInitialization(Config.GetConfig())
    );

    if(Config.GetConfig()["manager"]["mode"].as<std::string>() == "interactive")
    {
        UIExecutive = new G4UIExecutive(argc, argv);
        auto VisManager = std::make_shared<G4VisExecutive>();
        VisManager->Initialize();
        RunManager->Initialize();
        G4UImanager* UIManager = G4UImanager::GetUIpointer();
        if(Config.GetConfig()["manager"]["macro"]) {
            UIManager->ApplyCommand("/control/execute " + 
                Config.GetConfig()["manager"]["macro"].as<std::string>()
            );
        }
        else {
            UIManager->ApplyCommand("/control/execute macros/vis.mac");
        }
        UIManager->ApplyCommand("/run/verbose 1");
        UIManager->ApplyCommand("/event/verbose 0");
        UIExecutive->SessionStart();
    }
    else if(Config.GetConfig()["manager"]["mode"].as<std::string>() == "macro")
    {
        UIExecutive = new G4UIExecutive(argc, argv);
        RunManager->Initialize();
        G4UImanager* UIManager = G4UImanager::GetUIpointer();
        if(Config.GetConfig()["manager"]["macro"]) {
            UIManager->ApplyCommand("/control/execute " + 
                Config.GetConfig()["manager"]["macro"].as<std::string>()
            );
        }
        UIManager->ApplyCommand("/run/verbose 1");
        UIManager->ApplyCommand("/event/verbose 0");
        UIExecutive->SessionStart();
    }
    else if(Config.GetConfig()["manager"]["mode"].as<std::string>() == "terminal")
    {
        auto VisManager = std::make_shared<G4VisExecutive>();
        VisManager->Initialize();
        RunManager->Initialize();
        G4UImanager* UIManager = G4UImanager::GetUIpointer();
#ifdef ARTIE_USE_UI_TCSH
        auto UITerminal = std::make_shared<G4UIterminal>(new G4UItcsh);
#else   
        auto UITerminal = std::make_shared<G4UIterminal>();
#endif
        UITerminal->SessionStart();
    }
    else
    {
        RunManager->Initialize();
        G4UImanager* UIManager = G4UImanager::GetUIpointer();
        UIManager->ApplyCommand("/run/verbose 1");
        UIManager->ApplyCommand("/event/verbose 0");
        G4String command = "/run/beamOn ";
        for(G4int run = 0; run < Config.GetConfig()["manager"]["number_of_runs"].as<G4int>(); run++)
        {
            UIManager->ApplyCommand(
                command + Config.GetConfig()["manager"]["number_of_events"].as<std::string>()
            );
        }
    }
}