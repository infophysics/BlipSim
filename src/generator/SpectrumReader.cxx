#include "SpectrumReader.hh"

namespace Blip
{


    SpectrumReader::SpectrumReader(const std::string& radName, const std::string& rootFilePath, const std::string& fileName)
    {
        std::string theFullFileName = rootFilePath + fileName;
        std::cout << "Reading in spectrum for " << radName.c_str() << " from : " << rootFilePath.c_str() << "\n";

        Bool_t addStatus = TH1::AddDirectoryStatus();
        TH1::AddDirectory(kFALSE);

        TFile mRootFile(theFullFileName.c_str(), "READ");
        TGraph* mEnergySpectra = new TGraph();
        if (mRootFile.IsOpen()) {
            mEnergySpectra = (TGraph*)mRootFile.Get("Betas"); // -- TGraph name within the spectrum root file (for beta decay isotopes)
            std::cout << "Successfully opened the root file" << "\n";
            mRandom = new TRandom3();
        } else {
            G4cerr << "Could not open the file!\n";
            return;
        }

        int numPoints = mEnergySpectra->GetN();
        double* y = mEnergySpectra->GetY();

        // -- The TGraphs from LArSoft are the rates vs. KE and the y-axis has units of: #/(keV * s)
        // -- Note: that G4 uses MeV and ns, so we'll have to convert
        std::string title=radName + " Beta decay spectrum";
        auto samplingHist = std::make_unique<TH1D>(radName.c_str(), title.c_str(), numPoints, 0, numPoints);
        for (int i=0; i<numPoints; i++)
        {
            samplingHist->SetBinContent(i+1, y[i]);
            samplingHist->SetBinError(i+1, 0.);
        }
        sIntegral = samplingHist->Integral();
        mSamplingHist = std::move(samplingHist);
        assert(!samplingHist);

        // -- close file
        delete mEnergySpectra;
        mRootFile.Close();
        
        TH1::AddDirectory(addStatus);
    }

    SpectrumReader::~SpectrumReader()
    {
        delete mRandom;
        //if (mRootFile) mRootFile->Close();
    }


    bool SpectrumReader::SampleTheRadiological(const G4double rate)
    {
        bool success = false;

        // -- calculate the number of decays to be produced
        double meanNumber = rate * 1.0;
        mDecays = mRandom->Poisson(meanNumber);


        return success;
    }

    //double SpectrumReader::GetRandomEnergy(const std::string& radName, int& pdg, double& KE, double& mass, double& momentum)
    double SpectrumReader::GetRandomEnergy(const std::string& radName)
    {
        // -- will eventually add more rads
        if (radName != "Argon_39") return 0.0;

        if (mSamplingHist == nullptr) return 0.0;

        // -- sample form the TGraph
        double t = (*mSamplingHist).GetRandom(); // -- keV
        double KE = t / 1.0E3; // -- convert to MeV, the default unit of G4

        return KE;
    }

    bool SpectrumReader::GenNDecays(const std::string& radName, long N)
    {
        // -- will eventually add more rads
        if (radName != "Argon_39") return false;

        if (mSamplingHist == nullptr) return false;

        double nEntries = (*mSamplingHist).GetEntries();
        std::cout << "TH1D has " << nEntries << std::endl;
        for (int i=0; i<N; i++) {
            // -- sample form the TGraph
            double t = (*mSamplingHist).GetRandom(); // -- keV
            double KE = t / 1.0E3; // -- convert to MeV, the default unit of G4
            mDecayKEs.push_back(KE);
            std::cout << "Generating beta with decay energy of : " << t << " [keV]" << std::endl;
        }
        return true;
    }

    bool SpectrumReader::GenNDecays(const std::string& radName, long N, std::vector<double>& outputKEs)
    {
        // -- will eventually add more rads
        if (radName != "Argon_39") return false;

        if (mSamplingHist == nullptr) return false;

        double nEntries = (*mSamplingHist).GetEntries();
        std::cout << "TH1D has " << nEntries << std::endl;
        //std::vector<double> mDecayKEs(N, 0.0);
        for (int i=0; i<N; i++) {
            // -- sample form the TGraph
            double t = (*mSamplingHist).GetRandom(); // -- keV
            double KE = t / 1.0E3; // -- convert to MeV, the default unit of G4
            outputKEs.push_back(KE);
            std::cout << "Generating beta with decay energy of : " << t << " [keV] " << std::endl;
        }
        return true;
    }

    TLorentzVector SpectrumReader::GetDirection(double p, double m)
    {
        // -- make it isotropic
        double cosTheta = (2. * mRandom->Uniform(0., 1.) - 1.);
        if (cosTheta < -1.0) cosTheta = -1.;
        if (cosTheta > 1.0) cosTheta = 1.;
        double const sinTheta = sqrt(1. - cosTheta * cosTheta);
        double const phi = 2. * TMath::Pi() * mRandom->Uniform(0., 1.);
        
        TLorentzVector theDirection(p*sinTheta*std::cos(phi), p*sinTheta*std::sin(phi), p*cosTheta, std::sqrt(p*p + m*m));

        return theDirection;
    }
}