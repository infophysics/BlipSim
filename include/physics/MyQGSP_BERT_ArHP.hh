#pragma once

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

namespace Blip
{

	class MyQGSP_BERT_ArHP: public G4VModularPhysicsList
	{
	public:
		MyQGSP_BERT_ArHP(G4int ver=1);
		virtual ~MyQGSP_BERT_ArHP()=default;

		MyQGSP_BERT_ArHP(const MyQGSP_BERT_ArHP &) = delete;
		MyQGSP_BERT_ArHP & operator=(const MyQGSP_BERT_ArHP &)=delete;

	// SetCuts()
	virtual void SetCuts();
	};

}