/**
 * @file Core.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-09-12
 */
#include "Core.hh"

namespace Blip
{

    ArrayManager::ArrayManager()
    {
    }

    ArrayManager* ArrayManager::sInstance = nullptr;

    ArrayManager* ArrayManager::GetInstance()
    {
        if(sInstance == nullptr) {
            sInstance = new ArrayManager();
        }
        return sInstance;
    }

    G4String GetVolumeMaterialName(const G4Step* step)
    {
        G4String name = "none";
        const G4StepPoint *preStepPoint = step->GetPreStepPoint();
        const G4StepPoint *postStepPoint = step->GetPostStepPoint();
        
        G4LogicalVolume* preVolume = 0;
        G4VPhysicalVolume* prePhysical = 0;
        G4Material* preMaterial = 0;
        const G4TouchableHandle& preTouchable = preStepPoint->GetTouchableHandle();
        if(preTouchable) {
            prePhysical = preTouchable->GetVolume();
        }
        if(prePhysical) {
            preVolume = prePhysical->GetLogicalVolume();
        }
        if(preVolume) { 
            preMaterial = preVolume->GetMaterial();
        }
        if(preMaterial) {
            return preMaterial->GetName();
        }

        G4LogicalVolume* postVolume = 0;
        G4VPhysicalVolume* postPhysical = 0;
        G4Material* postMaterial = 0;
        const G4TouchableHandle& postTouchable = postStepPoint->GetTouchableHandle();
        if(postTouchable) {
            postPhysical = postTouchable->GetVolume();
        }
        if(postPhysical) {
            postVolume = postPhysical->GetLogicalVolume();
        }
        if(postVolume) { 
            postMaterial = postVolume->GetMaterial();
        }
        if(postMaterial) {
            return postMaterial->GetName();
        }
        return name;
    }

    G4String GetVolumeName(const G4Step* step)
    {
        G4String name = "none";
        const G4StepPoint *preStepPoint = step->GetPreStepPoint();
        const G4StepPoint *postStepPoint = step->GetPostStepPoint();
        
        G4LogicalVolume* preVolume = 0;
        G4VPhysicalVolume* prePhysical = 0;
        const G4TouchableHandle& preTouchable = preStepPoint->GetTouchableHandle();
        if(preTouchable) {
            prePhysical = preTouchable->GetVolume();
        }
        if(prePhysical) {
            preVolume = prePhysical->GetLogicalVolume();
        }
        if(preVolume) { 
            return preVolume->GetName();
        }

        G4LogicalVolume* postVolume = 0;
        G4VPhysicalVolume* postPhysical = 0;
        const G4TouchableHandle& postTouchable = postStepPoint->GetTouchableHandle();
        if(postTouchable) {
            postPhysical = postTouchable->GetVolume();
        }
        if(postPhysical) {
            postVolume = postPhysical->GetLogicalVolume();
        }
        if(postVolume) { 
            return postVolume->GetName();
        }

        return name;
    }

    G4String GetPostProcessName(const G4Step* step)
    {
        G4String name = "none";
        const G4StepPoint *postStepPoint = step->GetPostStepPoint();
        const G4VProcess* process = postStepPoint->GetProcessDefinedStep();
        if(process) {
            return process->GetProcessName();
        }
        return name;
    }

#ifdef BLIP_YAML
    ConfigParser::ConfigParser(G4String configurationFile)
    : mConfigurationFile(configurationFile)
    {
        mConfig = YAML::LoadFile(mConfigurationFile);
        if(!mConfig["manager"]) 
        {
            std::cout << "CONFIG ERROR! the section ";
            std::cout << "'manager' does not exist!" << std::endl;
            exit(0);
        }
        if(!mConfig["generator"]) 
        {
            std::cout << "CONFIG ERROR! the section ";
            std::cout << "'generator' does not exist!" << std::endl;
            exit(0);
        }
        if(!mConfig["manager"]["number_of_threads"]) 
        {
            std::cout << "CONFIG ERROR! the parameter ";
            std::cout << "'manager:number_of_threads' does not exist!" << std::endl;
            exit(0);
        }
        if(!mConfig["manager"]["number_of_runs"]) 
        {
            std::cout << "CONFIG ERROR! the parameter ";
            std::cout << "'manager:number_of_runs' does not exist!" << std::endl;
            exit(0);
        }
        if(!mConfig["manager"]["number_of_events"]) 
        {
            std::cout << "CONFIG ERROR! the parameter ";
            std::cout << "'manager:number_of_events' does not exist!" << std::endl;
            exit(0);
        }
        if(!mConfig["manager"]["mode"]) 
        {
            std::cout << "CONFIG ERROR! the parameter ";
            std::cout << "'manager:mode' does not exist!" << std::endl;
            exit(0);
        }
        if(!mConfig["manager"]["output_filename"]) 
        {
            std::cout << "CONFIG ERROR! the parameter ";
            std::cout << "'manager:output_filename' does not exist!" << std::endl;
            exit(0);
        }
        if(!mConfig["manager"]["save_particle_maps"]) 
        {
            std::cout << "CONFIG ERROR! the parameter ";
            std::cout << "'manager:save_particle_maps' does not exist!" << std::endl;
            exit(0);
        }
        if(!mConfig["manager"]["save_primary_info"]) 
        {
            std::cout << "CONFIG ERROR! the parameter ";
            std::cout << "'manager:save_primary_info' does not exist!" << std::endl;
            exit(0);
        }
    }

    ConfigParser::~ConfigParser()
    {
    }
#endif 
}