#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4UImanager.hh"

#include "Particle.hh"

#include <vector>

namespace Blip {

    class SpectrumReader {
    public:
        SpectrumReader(const std::string& radName, const std::string& rootFilePath, const std::string& fileName);
        virtual ~SpectrumReader();

        bool SampleTheRadiological(const G4double rate); // -- determine the number of decays to be produced
        bool GenNDecays(const std::string& radName, long i); // -- generate the kinematics for the decays
        bool GenNDecays(const std::string& radName, long i, std::vector<double>& outputKEs); // -- generate the kinematics for the decays
        //std::vector<double> getKineticEnergies() { return move(mDecayKEs); }
        std::vector<double> getKineticEnergies() { return mDecayKEs; }

    private:
        //std::vector<Particle*>* GetNDecays(long N); // -- generate the kinematics for the decays
        TLorentzVector GetDirection(double p, double m); // -- get the decay product's direction
        //double GetRandomEnergy(const std::string& radName, int& pdg, double& KE, double& mass, double& momentum); // -- sample from the spectrum randomly
        double GetRandomEnergy(const std::string& radName);

        //TFile* mRootFile;
        //TGraph* mEnergySpectra; // -- TGraph inside file
        TRandom3* mRandom; // -- rng

        //std::vector<Particle*>* primDecaysParticles;

        std::unique_ptr<TH1D> mSamplingHist; // -- Histogram to be sampled from
        //std::unique_ptr<std::vector<double>> mDecayKEs;
        std::vector<double> mDecayKEs;

        double sIntegral; // -- integral for sampling hist
        int mDecays; // -- number of decays to be generated
    };

}