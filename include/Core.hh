/**
 * @file Core.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-27
 */
#pragma once
#include <string>
#include <map>
#include <chrono>
#include <iostream>
#include <sys/time.h>
#include <ctime>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Step.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#ifdef BLIP_YAML
#include "yaml-cpp/yaml.h"
#endif

namespace Blip
{

    inline G4double NeutronMassSI()    { return 1.674927498e-27 * kg; }
    inline G4double NeutronMassMeV()   { return 939.56542052 * MeV; }
    inline G4double keVToJoules(G4double energy) { return energy * 1.60218e-16; }
    inline G4double SpeedOfLight()  { return 299792458 * m / s; }

    inline G4int GetTimeInMilliseconds() { 
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
        ).count(); 
    }

    class ArrayManager
    {
    private:
        static ArrayManager* sInstance;

    protected:
        ArrayManager();

    public:
        ArrayManager(ArrayManager& other) = delete;
        void operator=(const ArrayManager&) = delete;

        static ArrayManager* GetInstance();
    };


    struct PrimaryGeneration
    {
        G4ParticleDefinition* definition = {nullptr};
        G4String name = "geantino";
        G4int Z = 0;
        G4int A = 0;
        G4bool ion = false;
        G4double charge = 0.0*eplus;
        G4double time = 0.;
        G4ThreeVector position = {0.,0.,0.};
        G4double energy = 0.;
        G4ThreeVector momentum_direction = {0.,0.,0.};

        PrimaryGeneration(
            G4String particle,
            G4double t,
            G4ThreeVector pos,
            G4double e, 
            G4ThreeVector mom_dir
        )
        : name(particle), time(t), position(pos)
        , energy(e), momentum_direction(mom_dir)
        {
            auto ParticleTable = G4ParticleTable::GetParticleTable();
            definition = ParticleTable->FindParticle(particle);
            charge = definition->GetPDGCharge();
        }
        PrimaryGeneration(
            G4String particle,
            G4int z,
            G4int a,
            G4double q,
            G4double t,
            G4ThreeVector pos,
            G4double e,
            G4ThreeVector mom_dir
        )
        : name(particle), Z(z), A(a), charge(q), time(t)
        , position(pos), energy(e), momentum_direction(mom_dir)
        {
            ion = true;
        }
    };

    struct PrimaryData
    {
        G4String name = "geantino";
        G4int track_id = -1;
        G4int pdg = 0;
        G4double init_t = -1;
        G4double init_x = 0;
        G4double init_y = 0;
        G4double init_z = 0;
        G4double init_energy = -1;
        G4double final_t = -1;
        G4double final_x = 0;
        G4double final_y = 0;
        G4double final_z = 0;
        G4double final_energy = -1;
        G4int number_of_daughters = 0;
        G4double total_daughter_init_energy = 0;
        G4double total_daughter_final_energy = 0;
        G4double total_daughter_edep = 0;
        G4int number_of_edeps = 0;
        G4double total_edep = 0;
        G4int number_of_hits = 0;

        PrimaryData(
            G4String _name, G4int _track_id, G4int _pdg,
            G4double _init_t, G4double _init_x, G4double _init_y,
            G4double _init_z, G4double _init_energy
        )
        : name(_name), track_id(_track_id), pdg(_pdg)
        , init_t(_init_t), init_x(_init_x), init_y(_init_y)
        , init_z(_init_z), init_energy(_init_energy)
        {
        }

        PrimaryData(){}
        PrimaryData(G4int _track_id)
        : track_id(_track_id)
        {}

        void AddFinalTrackData(
            G4double _final_t, G4double _final_x, G4double _final_y,
            G4double _final_z, G4double _final_energy
        )
        {
            final_t = _final_t;
            final_x = _final_x;
            final_y = _final_y;
            final_z = _final_z;
            final_energy = _final_energy;
        }
    };

    G4String GetVolumeName(const G4Step* step);
    G4String GetPostProcessName(const G4Step* step);
    G4String GetVolumeMaterialName(const G4Step* step);

#ifdef BLIP_YAML
    class ConfigParser
    {
    public:
        ConfigParser(G4String configurationFile);
        ~ConfigParser();

        YAML::Node GetConfig() const { return mConfig; }

    private:
        G4String mConfigurationFile = {""};
        YAML::Node mConfig;
    };
#endif 

    // function for computing the nominal tof
    // from energy and length
    G4double GetNominalTOF(G4double energy, G4double length);

}