/**
 * @file NeutronHPPhysics.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-27
 */
#include "NeutronHPPhysics.hh"

namespace Blip
{
    NeutronHPPhysics::NeutronHPPhysics(const G4String& name)
    :  G4VPhysicsConstructor(name)
    , mThermal(true)
    {
        std::filesystem::create_directory(".physics");
    }

    NeutronHPPhysics::~NeutronHPPhysics()
    {
    }

    void NeutronHPPhysics::ConstructParticle()
    {
        G4Neutron::Neutron();
    }

    void NeutronHPPhysics::ConstructProcess()
    {
        G4ParticleDefinition* neutron = G4Neutron::Neutron();
        G4ProcessManager* pManager = neutron->GetProcessManager();
        G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

        // delete all neutron processes if already registered
        G4VProcess* process = 0;
        process = processTable->FindProcess("hadElastic", neutron);      
        if (process) {
            pManager->RemoveProcess(process);
        }
        process = processTable->FindProcess("neutronInelastic", neutron);
        if (process) {
            pManager->RemoveProcess(process);
        }
        process = processTable->FindProcess("nCapture", neutron);      
        if (process) {
            pManager->RemoveProcess(process);
        }
        process = processTable->FindProcess("nFission", neutron);      
        if (process) {
            pManager->RemoveProcess(process);
        }
                        
        // Elastic Processes
        G4HadronElasticProcess* process1 = new G4HadronElasticProcess();
        G4ParticleHPElastic* model1a = new G4ParticleHPElastic();
        G4ParticleHPElasticData* dataSet1a = new G4ParticleHPElasticData();

        // output the elastic cross section data
        dataSet1a->BuildPhysicsTable(*neutron);
        std::size_t numberOfElements = G4Element::GetNumberOfElements();
        static G4ThreadLocal G4ElementTable *theElementTable = 0 ;
        if (!theElementTable) theElementTable= G4Element::GetElementTable();
        auto ElasticCrossSections = G4ParticleHPManager::GetInstance()->GetElasticCrossSections();
        for (std::size_t ii = 0; ii < numberOfElements; ii++)
        {
            std::ofstream ElasticOutputFile;
            ElasticOutputFile.open(".physics/" + (*theElementTable)[ii]->GetName() + "_elastic.csv");
            for (G4int ie = 0 ; ie < 10000 ; ++ie )
            {
                G4double eKinetic = 1.0e-5 * G4Pow::GetInstance()->powA ( 10.0 , ie/1000.0 ) *eV;
                G4bool outOfRange = false;

                if ( eKinetic < 20*MeV )
                {
                    ElasticOutputFile << eKinetic/eV << ",";
                    ElasticOutputFile << (*((*ElasticCrossSections)(ii))).GetValue(eKinetic, outOfRange)/barn << "\n";
                }
            }
            ElasticOutputFile.close();
            std::ofstream ElasticElementFile;
            ElasticElementFile.open(".physics/" + (*theElementTable)[ii]->GetName() + "_properties.csv");
            ElasticElementFile << "Property,Value\n";
            ElasticElementFile << "Symbol," << (*theElementTable)[ii]->GetSymbol() << "\n";
            ElasticElementFile << "Z," << std::to_string((*theElementTable)[ii]->GetZ()) << "\n";
            ElasticElementFile << "N," << std::to_string((*theElementTable)[ii]->GetN()) << "\n";
            ElasticElementFile << "AtomicMassAmu," << std::to_string((*theElementTable)[ii]->GetAtomicMassAmu()) << "\n";
            ElasticElementFile << "A," << std::to_string((*theElementTable)[ii]->GetA()) << "\n";
            ElasticElementFile << "NbOfAtomicShells," << std::to_string((*theElementTable)[ii]->GetNbOfAtomicShells()) << "\n";
            G4IsotopeVector* IsotopeVector = (*theElementTable)[ii]->GetIsotopeVector();
            G4double* RelativeAbundanceVector = (*theElementTable)[ii]->GetRelativeAbundanceVector();
            for(G4int jj = 0; jj < IsotopeVector->size(); jj++)
            {
                ElasticElementFile << "Name," << (*IsotopeVector)[jj]->GetName() << "\n";
                ElasticElementFile << "Z," << std::to_string((*IsotopeVector)[jj]->GetZ()) << "\n";
                ElasticElementFile << "N," << std::to_string((*IsotopeVector)[jj]->GetN()) << "\n";
                ElasticElementFile << "A," << std::to_string((*IsotopeVector)[jj]->GetA()) << "\n";
                ElasticElementFile << "RelativeAbundance," << std::to_string(RelativeAbundanceVector[jj]) << "\n";
            }
        }

        G4ParticleHPThermalScattering* model1b = new G4ParticleHPThermalScattering();
        G4ParticleHPThermalScatteringData* dataSet1b = new G4ParticleHPThermalScatteringData();
        dataSet1b->BuildPhysicsTable(*neutron);

        pManager->AddDiscreteProcess(process1);   
        process1->RegisterMe(model1a);
        process1->AddDataSet(dataSet1a);
        if (mThermal) 
        {
            model1a->SetMinEnergy(4*eV);
            process1->RegisterMe(model1b);
            process1->AddDataSet(dataSet1b);
        }
        

        // Inelastic Processes
#ifdef ARTIE_GEANT_10
        G4HadronInelasticProcess* process2 = new G4HadronInelasticProcess(
            "neutronInelastic", neutron
        );
#else
        G4HadronInelasticProcess* process2 = new G4HadronInelasticProcess(
            "neutronInelastic"
        );
#endif
        G4ParticleHPInelasticData* dataSet2 = new G4ParticleHPInelasticData();
        dataSet2->BuildPhysicsTable(*neutron);

        G4ParticleHPInelastic* model2 = new G4ParticleHPInelastic();
        pManager->AddDiscreteProcess(process2);   
        process2->AddDataSet(dataSet2);  
        process2->RegisterMe(model2);    

        // Neutron Capture
#ifdef ARTIE_GEANT_10   
        G4HadronCaptureProcess* process3 = new G4HadronCaptureProcess();
#else
        G4NeutronCaptureProcess* process3 = new G4NeutronCaptureProcess();
#endif
        G4ParticleHPCaptureData* dataSet3 = new G4ParticleHPCaptureData();
        dataSet3->BuildPhysicsTable(*neutron);
        auto CaptureCrossSections = G4ParticleHPManager::GetInstance()->GetCaptureCrossSections();
        for (std::size_t ii = 0; ii < numberOfElements; ii++)
        {
            std::ofstream CaptureFile;
            CaptureFile.open(".physics/" + (*theElementTable)[ii]->GetName() + "_capture.csv");
            for (G4int ie = 0 ; ie < 10000 ; ++ie )
            {
                G4double eKinetic = 1.0e-5 * G4Pow::GetInstance()->powA ( 10.0 , ie/1000.0 ) *eV;
                G4bool outOfRange = false;

                if ( eKinetic < 20*MeV )
                {
                    CaptureFile << eKinetic/eV << ",";
                    CaptureFile << (*((*CaptureCrossSections)(ii))).GetValue(eKinetic, outOfRange)/barn << "\n";
                }
            }
            CaptureFile.close();
        }

        G4ParticleHPCapture* model3 = new G4ParticleHPCapture();
        pManager->AddDiscreteProcess(process3);  
        process3->AddDataSet(dataSet3);
        process3->RegisterMe(model3);
        
        // Neutron Fission  
#ifdef ARTIE_GEANT_10
        G4HadronFissionProcess* process4 = new G4HadronFissionProcess();
#else
        G4NeutronFissionProcess* process4 = new G4NeutronFissionProcess();
#endif
        G4ParticleHPFissionData* dataSet4 = new G4ParticleHPFissionData();
        G4ParticleHPFission* model4 = new G4ParticleHPFission();
        pManager->AddDiscreteProcess(process4);
        process4->AddDataSet(dataSet4);     
        process4->RegisterMe(model4);
    }
}