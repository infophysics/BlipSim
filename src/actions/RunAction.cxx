/**
 * @file RunAction.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-12-13
 */
#include "RunAction.hh"

namespace Blip
{
    RunAction::RunAction()
    : G4UserRunAction()
    {
        EventManager::GetEventManager()->CreateTuples();
    }

    RunAction::~RunAction()
    {}

    void RunAction::BeginOfRunAction(const G4Run* run)
    {
        auto Manager = EventManager::GetEventManager();
        Manager->ResetProfiling();
        Manager->EvaluateRunBegin();
        Manager->OpenOutputFile(run->GetRunID());
    }

    void RunAction::EndOfRunAction(const G4Run* run)
    {
        auto Manager = EventManager::GetEventManager();
        Manager->EvaluateRunEnd();
        Manager->CloseOutputFile(run->GetRunID());
    }
}