/**
 * @file Analysis.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-27
 */
#pragma once
#include "EventManager.hh"

namespace Blip
{
    struct AnalysisTuple
    {
        G4double x_pos = 0;
        G4double y_pos = 0;
        G4double z_pos = 0;
        G4double n_energy = 0;
        G4int N_x = 0;
        G4int N_y = 0;
        G4int N_z = 0;
        G4double n_tof = 0;
        G4int n_track_id = 0;
        G4double x_tracking = 0;
        G4double y_tracking = 0;
        G4double z_tracking = 0;
    };

    void AnalysisFunctionRunBegin();
    void AnalysisFunctionRunEnd();
    void AnalysisFunctionEventBegin();
    void AnalysisFunctionEventEnd();
}