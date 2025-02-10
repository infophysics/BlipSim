#pragma once

#include "globals.hh"
#include "G4ParticleHPChannel.hh"
#include "G4HadronicInteraction.hh"
#include "G4NucleiProperties.hh"

namespace Blip
{
	class ArParticleHPCapture : public G4HadronicInteraction
	{
	public:

		ArParticleHPCapture();

		~ArParticleHPCapture();

		G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);

		virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

	public:
		G4int GetVerboseLevel() const;
		void SetVerboseLevel( G4int );
		void BuildPhysicsTable(const G4ParticleDefinition&);
		virtual void ModelDescription(std::ostream& outFile) const;

	private:

		std::vector<G4ParticleHPChannel*>* theCapture;
		G4String dirName;
		G4int numEle;

		G4HadFinalState theResult;
	};
}