/**
 * @file Particle.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.0
 * @date 2022-10-13
 */
#include "Particle.hh"

namespace Blip
{
    Particle::Particle()
    {
    }

    Particle::~Particle()
    {
    }

    Particle::Particle(
        G4int track_id, 
        G4int pdg,
        G4int parent_track_id,
        G4String creation_process,
        G4double global_time
    )
    : mTrackID(track_id)
    , mPDG(pdg)
    , mParentTrackID(parent_track_id)
    , mCreationProcess(creation_process)
    , mGlobalCreationTime(global_time)
    {
    }

    Particle::Particle(
        G4int track_id, 
        G4int pdg,
        G4int parent_track_id,
        G4String creation_process,
        G4double global_time,
        G4double t,
        G4ThreeVector position,
        G4double E,
        G4ThreeVector momentum
    )
    : mTrackID(track_id)
    , mPDG(pdg)
    , mParentTrackID(parent_track_id)
    , mCreationProcess(creation_process)
    , mGlobalCreationTime(global_time)
    {
        mTrajectory.AddTrajectoryPoint(t, position, E, momentum);
    }

    G4int Particle::DaughterTrackID(const G4int ii) const
    {
        if(ii >= mDaughterTrackIDs.size()) {
            throw std::invalid_argument(
                "Index " + std::to_string(ii) + " out of bounds for "
                + " daughter track id list with " + std::to_string(mDaughterTrackIDs.size())
                + " entries."
            );
        }
        return mDaughterTrackIDs[ii];
    }

    G4bool Particle::CheckTrajectoryIndex(const G4int ii) const
    {
        if(ii >= mTrajectory.size()) {
            throw std::invalid_argument(
                "Index " + std::to_string(ii) + " out of bounds for "
                + " trajectory with " + std::to_string(mTrajectory.size())
                + " entries."
            );
            return false;
        }
        return true;
    }

    G4double Particle::GetT(const G4int ii) const
    {
        if(CheckTrajectoryIndex(ii)) {
            return mTrajectory.t[ii];
        }
    }
    G4double Particle::GetX(const G4int ii) const
    {
        if(CheckTrajectoryIndex(ii)) {
            return mTrajectory.x[ii];
        }
    }
    G4double Particle::GetY(const G4int ii) const
    {
        if(CheckTrajectoryIndex(ii)) {
            return mTrajectory.y[ii];
        }
    }
    G4double Particle::GetZ(const G4int ii) const
    {
        if(CheckTrajectoryIndex(ii)) {
            return mTrajectory.z[ii];
        }
    }
    G4double Particle::GetE(const G4int ii) const
    {
        if(CheckTrajectoryIndex(ii)) {
            return mTrajectory.E[ii];
        }
    }
    G4double Particle::GetPx(const G4int ii) const
    {
        if(CheckTrajectoryIndex(ii)) {
            return mTrajectory.px[ii];
        }
    }
    G4double Particle::GetPy(const G4int ii) const
    {
        if(CheckTrajectoryIndex(ii)) {
            return mTrajectory.py[ii];
        }
    }
    G4double Particle::GetPz(const G4int ii) const
    {
        if(CheckTrajectoryIndex(ii)) {
            return mTrajectory.pz[ii];
        }
    }
}