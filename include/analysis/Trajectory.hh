/**
 * @file Trajectory.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.0
 * @date 2022-10-13
 */
#pragma once
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

namespace Blip
{
    struct Trajectory
    {
        std::vector<G4double> t = {};
        std::vector<G4double> x = {};
        std::vector<G4double> y = {};
        std::vector<G4double> z = {};
        std::vector<G4double> E = {};
        std::vector<G4double> px = {};
        std::vector<G4double> py = {};
        std::vector<G4double> pz = {};

        G4int size() const { 
            return t.size();
        }

        void AddTrajectoryPoint(
            G4double _t, G4ThreeVector _position,
            G4double _E, G4ThreeVector _momentum
        )
        {
            t.emplace_back(_t);
            x.emplace_back(_position[0]);
            y.emplace_back(_position[1]);
            z.emplace_back(_position[2]);
            E.emplace_back(_E);
            px.emplace_back(_momentum[0]);
            py.emplace_back(_momentum[1]);
            pz.emplace_back(_momentum[2]);
        }

        Trajectory()
        {
        }

    };
}