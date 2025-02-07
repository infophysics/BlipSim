/**
 * @file EnergyDeposit.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.0
 * @date 2022-10-13
 */
#pragma once

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

namespace Blip
{
    struct EnergyDeposit
    {
        G4int track_id = {0};
        G4String name = {"none"};
        G4String volume = {"none"};
        G4int pdg = {0};
        G4double local_time = {0};
        G4double global_time = {0};
        G4ThreeVector pre_step_position = {0,0,0};
        G4ThreeVector post_step_position = {0,0,0};
        G4double energy = {0};

        EnergyDeposit(
            G4int _track_id, G4String _volume, G4String _name,
            G4int _pdg, G4double _local_time, G4double _global_time,
            G4ThreeVector _pre_step_position, G4ThreeVector _post_step_position,
            G4double _energy
        )
        : track_id(_track_id), volume(_volume), name(_name)
        , pdg(_pdg), local_time(_local_time), global_time(_global_time)
        , pre_step_position(_pre_step_position), post_step_position(_post_step_position)
        , energy(_energy)
        {
        }
    };
}