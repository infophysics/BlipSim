/**
 * @file SteppingAction.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12-13
 */
#pragma once
#include <memory>
#include "G4UserSteppingAction.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

#include "EventManager.hh"

namespace Blip
{
    class SteppingAction : public G4UserSteppingAction
    {
    public:
        SteppingAction();
        ~SteppingAction();
        
        virtual void UserSteppingAction(const G4Step*);
        
    private:

    };
}