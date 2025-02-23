/**
 * @file BlipActionInitialization.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12-13
 */
#pragma once
#include "G4VUserActionInitialization.hh"

#include "BlipPrimaryGeneratorAction.hh"
#include "BoxRadGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "Core.hh"

#ifdef BLIP_YAML
#include "yaml-cpp/yaml.h"
#endif

namespace Blip
{
    class BlipActionInitialization : public G4VUserActionInitialization
    {
    public:
        BlipActionInitialization();
        ~BlipActionInitialization();

#ifdef BLIP_YAML
        BlipActionInitialization(YAML::Node config);
#endif

        virtual void Build() const;
        virtual void BuildForMaster() const;

    private:
#ifdef BLIP_YAML
        YAML::Node mConfig;
#endif
    };
}