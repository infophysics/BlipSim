/**
 * @file Profile.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-27
 */
#pragma once

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

namespace Blip
{
    struct ProfileEventData
    {
        G4String name = {"default"};
        G4String volume = {"volume"};
        G4int exit = {0};
        G4double z = {0};
        std::vector<G4double> x = {};
        std::vector<G4double> y = {};
    };
}