#include "TFile.h"
#include "TGraph.h"
#include "TRandom3.h"

#include <vector>

namespace Blip {

    class SpectrumReader {
    public:
        SpectrumReader(const std::string& rootFilePath);
        virtual ~SpectrumReader();

        double GetRandomEnergy(); // -- sample from the spectrum randomly
        std::vector<double> GetRandomEnergy(int nDecays); // -- sample n times

    private:
        TFile* mRootFile;
        TGraph* mEnergySpectra; // -- TGraph inside file
        TRandom3* mRandom; // -- rng
    };

}