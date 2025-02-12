#include "SpectrumReader.hh"

namespace Blip
{
    SpectrumReader::SpectrumReader(const std::string& rootFilePath)
    {
        mRootFile = TFile::Open("");
        if (mRootFile && mRootFile->IsOpen()) {
            mEnergySpectra = (TH1D*)mRootFile->Get("Ar30_TH1D"); // -- FIXME: Need to actually get the right names from the LArSoft files
            mRandom = new TRandom3();
        }
    }

    SpectrumReader::~SpectrumReader()
    {
        delete mRandom;
        if (mRootFile) mRootFile->Close();
    }

    double SpectrumReader::GetRandomEnergy()
    {
        if (!mEnergySpectra) return 0.0;

        // -- sample form the hist
        return mEnergySpectra->GetRandom();
    }
}