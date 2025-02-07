/**
 * @file BlipDetectorConstruction.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.0
 * @details 
 *  Change log:
 *      2022/09/20 - Initial creation of the file.
 * @date 2022-09-23
 */
#include "BlipDetectorConstruction.hh"

namespace Blip
{
    BlipDetectorConstruction::BlipDetectorConstruction()
    : G4VUserDetectorConstruction()
    {
        DefineMaterials();
    }

#ifdef BLIP_YAML
    BlipDetectorConstruction::BlipDetectorConstruction(YAML::Node config)
    : G4VUserDetectorConstruction()
    , mConfig(config)
    {

        if(mConfig["hall"]["world_material"])   { mWorldMaterialName = mConfig["hall"]["world_material"].as<std::string>(); }
        if(mConfig["hall"]["world_x"])          { mExperimentalHallX = mConfig["hall"]["world_x"].as<G4double>() * m; }
        if(mConfig["hall"]["world_y"])          { mExperimentalHallY = mConfig["hall"]["world_y"].as<G4double>() * m; }
        if(mConfig["hall"]["world_z"])          { mExperimentalHallZ = mConfig["hall"]["world_z"].as<G4double>() * m; }

        DefineMaterials();
    }
#endif

    void BlipDetectorConstruction::DefineMaterials()
    {
        mWorldMaterial = CreateMaterial(mWorldMaterialName, "World");
    }

    BlipDetectorConstruction::~BlipDetectorConstruction()
    {
    }

    G4VPhysicalVolume* BlipDetectorConstruction::Construct()
    {
        DefineMaterials();
        // create the world volume
        mSolidExperimentalHall = new G4Box(
            "Solid_BlipExperimentalHall", 
            mExperimentalHallX, 
            mExperimentalHallY, 
            mExperimentalHallZ
        );
        mLogicalExperimentalHall = new G4LogicalVolume(
            mSolidExperimentalHall, 
            mWorldMaterial, 
            "Logical_BlipExperimentalHall"
        );
        mPhysicalExperimentalHall = new G4PVPlacement(
            0, 
            G4ThreeVector(0., 0., 0.),
            mLogicalExperimentalHall,
            "Physical_BlipExperimentalHall",
            0, 
            false,
            0
        );

        return mPhysicalExperimentalHall;
    }

    void BlipDetectorConstruction::ConstructSDandField()
    {
    }
}