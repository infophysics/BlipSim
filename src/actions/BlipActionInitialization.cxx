/**
 * @file BlipActionInitialization.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-27
 */
#include "BlipActionInitialization.hh"

namespace Blip
{
    BlipActionInitialization::BlipActionInitialization()
    {
    }

    BlipActionInitialization::~BlipActionInitialization()
    {
    }

#ifdef BLIP_YAML
    BlipActionInitialization::BlipActionInitialization(YAML::Node config)
    : mConfig(config)
    {
    }
#endif

    void BlipActionInitialization::Build() const
    {
#ifdef BLIP_YAML
        SetUserAction(new BlipPrimaryGeneratorAction(mConfig));
#else
        SetUserAction(new BlipPrimaryGeneratorAction());
#endif 
        SetUserAction(new RunAction());
        SetUserAction(new EventAction());
        SetUserAction(new SteppingAction());
        SetUserAction(new TrackingAction());
    }

    void BlipActionInitialization::BuildForMaster() const
    {
        SetUserAction(new RunAction());
    }
}