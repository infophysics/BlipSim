#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom3.h"

namespace Blip {

    class SpectrumReader {
    public:
        SpectrumReader(const std::string& rootFilePath);
        virtual ~SpectrumReader();

        double GetRandomEnergy(); // -- sample from the spectrum randomly

    private:
        TFile* mRootFile;
        TH1D* mEnergySpectra; // -- histo inside file
        TRandom3* mRandom; // -- rng
    };

}