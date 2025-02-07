/**
 * @file Material.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-30
 */
#pragma once
#include <memory>

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"

#include "Argon.hh"
#include "EventManager.hh"

namespace Blip
{
    G4Material* CreateMaterial(G4String material, G4String volume);
}