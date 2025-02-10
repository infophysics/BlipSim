#pragma once

#include "ArCaptureGammas.hh"

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ReactionProductVector.hh"
#include "G4ParticleHPNames.hh"
#include "G4ParticleHPPhotonDist.hh"
#include "G4ParticleHPEnAngCorrelation.hh"


namespace Blip
{
	class ArParticleHPCaptureFS : public G4ParticleHPFinalState
	{
	public:
		bool useArCapGamma = true;

		ArParticleHPCaptureFS()
		{
			hasXsec = false;
			hasExactMF6 = false;
			targetMass = 0;
		}

		~ArParticleHPCaptureFS()
		{
		}

		void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType, G4ParticleDefinition* );
		G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack);
		G4ParticleHPFinalState * New()
		{
			ArParticleHPCaptureFS * theNew = new ArParticleHPCaptureFS;
			return theNew;
		}

	private:

		G4double targetMass;

		G4ParticleHPPhotonDist theFinalStatePhotons;
		ArCaptureGammas       theFinalgammas;

		G4ParticleHPEnAngCorrelation theMF6FinalState;
		G4bool hasExactMF6;

		G4ParticleHPNames theNames;

	//  G4double theCurrentA;
	//  G4double theCurrentZ;
	};
}