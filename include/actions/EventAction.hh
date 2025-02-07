/**
 * @file EventAction.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12-13
 */
#pragma once
#include <memory>
#include <map>
#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4Event.hh"
#ifdef BLIP_GEANT_10
#include "g4root.hh"
#else
#include "G4AnalysisManager.hh"
#endif
#include "G4RunManager.hh"

#include "RunAction.hh"
#include "EventManager.hh"

namespace Blip
{
    class EventAction : public G4UserEventAction
    {
    public:
        EventAction();
        ~EventAction();
        
        virtual void BeginOfEventAction(const G4Event*);
        virtual void EndOfEventAction(const G4Event*);

    private:

    };
}