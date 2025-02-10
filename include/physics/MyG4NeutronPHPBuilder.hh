#pragma once

// -- artg4tk includes
#include "ArParticleHPCapture.hh"

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4NeutronFissionProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPFission.hh"
#include "G4ParticleHPFissionData.hh"
//#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4NucleiProperties.hh"

namespace Blip
{

	class MyG4NeutronPHPBuilder : public G4VNeutronBuilder
	{
	public:
		MyG4NeutronPHPBuilder();
		virtual ~MyG4NeutronPHPBuilder() {}

	public:
		virtual void Build(G4HadronElasticProcess * aP) final override;
		virtual void Build(G4NeutronFissionProcess * aP) final override;
		virtual void Build(G4NeutronCaptureProcess * aP) final override;
		virtual void Build(G4HadronInelasticProcess * aP) final override;

		virtual void SetMinEnergy(G4double aM) final override
		{
		theMin=aM;
		theIMin = theMin;
		}
		void SetMinInelasticEnergy(G4double aM)
		{
		theIMin=aM;
		}
		virtual void SetMaxEnergy(G4double aM) final override
		{
		theIMax = aM;
		theMax=aM;
		}
		void SetMaxInelasticEnergy(G4double aM)
		{
		theIMax = aM;
		}

		using G4VNeutronBuilder::Build; //Prevent compiler warning

	private:

		G4double theMin;
		G4double theIMin;
		G4double theMax;
		G4double theIMax;

		G4ParticleHPElastic * theHPElastic;
		G4ParticleHPElasticData * theHPElasticData;
		G4ParticleHPInelastic * theHPInelastic;
		G4ParticleHPInelasticData * theHPInelasticData;
		G4ParticleHPFission * theHPFission;
		G4ParticleHPFissionData * theHPFissionData;
		ArParticleHPCapture * theHPCapture;
		G4ParticleHPCaptureData * theHPCaptureData;

	};
}
