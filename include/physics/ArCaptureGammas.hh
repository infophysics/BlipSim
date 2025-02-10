#pragma once

#include "G4ReactionProductVector.hh"
#include <vector>


namespace Blip
{
	class ArCaptureGammas
	{
	public:

		ArCaptureGammas();
		~ArCaptureGammas();
		G4ReactionProductVector * GetGammas ();
		std::vector<double>  Initialize ();
		std::vector<double> CapAr40();
		std::vector<double> continuum();


	public:
		double Elevel;
		double  xint[4][750];

	};
}
