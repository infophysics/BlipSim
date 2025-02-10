#pragma once

#include "G4HadronPhysicsQGSP_BERT.hh"


namespace Blip
{
	class MyG4HadronPhysicsQGSP_BERT_ArHP : public G4HadronPhysicsQGSP_BERT
	{
	public:
		MyG4HadronPhysicsQGSP_BERT_ArHP(G4int verbose =1);
		MyG4HadronPhysicsQGSP_BERT_ArHP(const G4String& name, G4bool quasiElastic=true);
		virtual ~MyG4HadronPhysicsQGSP_BERT_ArHP() {}

	protected:
		virtual void Neutron() override;
		virtual void ExtraConfiguration();
	};
}
