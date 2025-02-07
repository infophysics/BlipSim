/**
 * @file SteppingAction.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12-13
 */
#include "SteppingAction.hh"

namespace Blip
{
    SteppingAction::SteppingAction()
    : G4UserSteppingAction()
    {
    }

    SteppingAction::~SteppingAction()
    {}

    void SteppingAction::UserSteppingAction(const G4Step *step)
    {   
        if(step->GetTrack()->GetGlobalTime() >= EventManager::GetEventManager()->EventMaxTime()) {
            step->GetTrack()->SetTrackStatus(fStopAndKill);
        }
        auto Manager = EventManager::GetEventManager();
        Manager->AddParticleInfoFromStep(step);
        Manager->AddEnergyDepositInfoFromStep(step);
    }
}