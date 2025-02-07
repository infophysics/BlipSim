/**
 * @file Particle.hh
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.0
 * @date 2022-10-13
 */
#pragma once
#include <string>
#include <memory>
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "Trajectory.hh"

namespace Blip
{
    class Particle
    {
    public:
        Particle();
        ~Particle();

        Particle(
            G4int track_id, 
            G4int pdg,
            G4int parent_track_id,
            G4String creation_process,
            G4double global_time
        );

        Particle(
            G4int track_id, 
            G4int pdg,
            G4int parent_track_id,
            G4String creation_process,
            G4double global_time,
            G4double t,
            G4ThreeVector position,
            G4double E,
            G4ThreeVector momentum
        );

        G4bool CheckTrajectoryIndex(const G4int ii) const;
        void AddTrajectoryPoint(
            G4double t, G4ThreeVector position,
            G4double E, G4ThreeVector momentum
        )
        { return mTrajectory.AddTrajectoryPoint(t, position, E, momentum); }

        void SetTrackID(G4int track_id)     { mTrackID = track_id; }
        void SetPDG(G4int pdg)              { mPDG = pdg; }
        void SetParentTrackID(G4int parent_track_id)    { mParentTrackID = parent_track_id; }
        void SetCreationProcess(G4String creation_process)  { mCreationProcess = creation_process; }
        void SetEndProcess(G4String end_process)            { mEndProcess = end_process; }

        G4int TrackID() const   { return mTrackID; }
        G4int PDG() const       { return mPDG; }
        G4int ParentTrackID() const         { return mParentTrackID; }
        G4String CreationProcess() const    { return mCreationProcess; }
        G4String EndProcess() const         { return mEndProcess; }

        G4int DaughterTrackID(const G4int ii=0) const;

        G4double GetGlobalCreationTime() const { return mGlobalCreationTime; }

        Trajectory GetTrajectory() const {return mTrajectory;}

        G4double GetT(const G4int ii=0) const;
        G4double GetX(const G4int ii=0) const;
        G4double GetY(const G4int ii=0) const;
        G4double GetZ(const G4int ii=0) const;
        G4double GetE(const G4int ii=0) const;
        G4double GetPx(const G4int ii=0) const;
        G4double GetPy(const G4int ii=0) const;
        G4double GetPz(const G4int ii=0) const;

    private:
        G4int mTrackID = {0};
        G4int mPDG = {0};
        G4int mParentTrackID = {0};
        G4String mCreationProcess = {""};
        G4String mEndProcess = {""};
        G4double mGlobalCreationTime = {0};
        std::vector<G4int> mDaughterTrackIDs = {};
        Trajectory mTrajectory;

    };
}