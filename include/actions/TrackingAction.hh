/**
 * @file TrackingAction.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12-13
 */
#pragma once
#include "G4UserTrackingAction.hh"
#include "G4Run.hh"
#include "G4GenericMessenger.hh"
#ifdef BLIP_GEANT_10
#include "g4root.hh"
#else
#include "G4AnalysisManager.hh"
#endif

#include "EventManager.hh"

namespace Blip
{
    class TrackingAction : public G4UserTrackingAction
    {
    public:
        TrackingAction();
        ~TrackingAction();
        
        virtual void PreUserTrackingAction(const G4Track*);
        virtual void PostUserTrackingAction(const G4Track*);
        
    private:

    };
}