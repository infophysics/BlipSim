/**
 * @file EventManager.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 */
#include "EventManager.hh"

namespace Blip
{
    std::shared_ptr<EventManager> EventManager::sInstance;
    std::mutex EventManager::sMutex;

    EventManager::EventManager()
    {
    }

    EventManager::~EventManager()
    {
    }

    void EventManager::SetConfig(YAML::Node config)
    {
        mConfig = config;
        if(mConfig["manager"])
        {
            if(mConfig["manager"]["number_of_threads"]) { 
                sNumberOfThreads = mConfig["manager"]["number_of_threads"].as<G4int>(); 
            }
            if(mConfig["manager"]["number_of_runs"])    { 
                sNumberOfRuns = mConfig["manager"]["number_of_runs"].as<G4int>(); 
            }
            if(mConfig["manager"]["number_of_events"])  { 
                sNumberOfEvents = mConfig["manager"]["number_of_events"].as<G4int>(); 
            }
            if(mConfig["manager"]["output_filename"])   { 
                sOutputFileName = mConfig["manager"]["output_filename"].as<std::string>(); 
            }
            if(mConfig["manager"]["save_particle_maps"]) { 
                sSaveParticleMaps = mConfig["manager"]["save_particle_maps"].as<G4bool>(); 
            }
            if(mConfig["manager"]["save_primary_info"]) { 
                sSavePrimaryInfo = mConfig["manager"]["save_primary_info"].as<G4bool>(); 
            }
            if(mConfig["manager"]["save_particle_info"]) { 
                sSaveParticleInfo = mConfig["manager"]["save_particle_info"].as<G4bool>(); 
            }
            if(mConfig["manager"]["save_energy_deposits"]) { 
                sSaveEnergyDeposits = mConfig["manager"]["save_energy_deposits"].as<G4bool>(); 
            }
            if(mConfig["manager"]["save_profile_data"]) { 
                sSaveProfileData = mConfig["manager"]["save_profile_data"].as<G4bool>(); 
            }
            if(mConfig["argon"]["use_g4_definition"])   { 
                mUseG4Definition = mConfig["argon"]["use_g4_definition"].as<G4bool>(); 
            }
            if(mConfig["argon"]["argon_36_ratio"])      { 
                mArgon36Ratio = mConfig["argon"]["argon_36_ratio"].as<G4double>(); 
            }
            if(mConfig["argon"]["argon_38_ratio"])      { 
                mArgon38Ratio = mConfig["argon"]["argon_38_ratio"].as<G4double>(); 
            }
            if(mConfig["argon"]["argon_40_ratio"])      { 
                mArgon40Ratio = mConfig["argon"]["argon_40_ratio"].as<G4double>(); 
            }
            if(mConfig["argon"]["lar_density"])         { 
                mLArDensity = mConfig["argon"]["lar_density"].as<G4double>() * g/cm3; 
            }
            if(mConfig["argon"]["lar_temperature"])     { 
                mLArTemperature = mConfig["argon"]["lar_temperature"].as<G4double>() * kelvin; 
            }
            if(mConfig["argon"]["lar_pressure"])        { 
                mLArPressure = mConfig["argon"]["lar_pressure"].as<G4double>() * atmosphere; 
            }
        }
        if(mConfig["hall"]["world_x"]) { mHallX = mConfig["hall"]["world_x"].as<G4double>() * m; }
        if(mConfig["hall"]["world_y"]) { mHallY = mConfig["hall"]["world_y"].as<G4double>() * m; }
        if(mConfig["hall"]["world_z"]) { mHallZ = mConfig["hall"]["world_z"].as<G4double>() * m; }

        G4GDMLParser* mGDMLParser;

    }

    void EventManager::SaveGDML()
    {
        G4GDMLParser* mGDMLParser;
        mGDMLParser->SetRegionExport(true);
        mGDMLParser->SetOutputFileOverwrite(true);
        mGDMLParser->Write(
            "gdml/" + OutputFileName() + ".gdml", 
            G4TransportationManager::GetTransportationManager()
            ->GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume()
        );
    }

    G4int EventManager::GetIndex(G4String tuple)
    {
        for(size_t ii = 0; ii < sTuples.size(); ii++)
        {
            if(sTuples[ii].name == tuple) {
                return ii;
            }
        }
        sCurrentTupleIndex += 1;
        sTuples.emplace_back(Tuple(tuple, sCurrentTupleIndex));
        return sCurrentTupleIndex;
    }

    void EventManager::OpenOutputFile(G4int RunID)
    {
        std::filesystem::create_directory("outputs");
        auto AnalysisManager = G4AnalysisManager::Instance();
        G4String OutputFile = "outputs/" + OutputFileName() + "_" + std::to_string(RunID) + ".root";
        G4bool fileopen = AnalysisManager->OpenFile(OutputFile);
        if (!fileopen) {
            G4cout << "File - " + OutputFile 
                << " - not opened!" << G4endl;
        }
        else {
            G4cout << "File - " + OutputFile
                << " - opened successfully." << G4endl;
        }
    }

    void EventManager::CloseOutputFile(G4int RunID)
    {
        auto AnalysisManager = G4AnalysisManager::Instance();
        AnalysisManager->Write();
        AnalysisManager->CloseFile();

#ifdef ARTIE_PROFILING
        std::ofstream ProfilingFile;
        std::filesystem::create_directory(".logs");
        auto ThreadID = G4Threading::G4GetThreadId();
        if(ThreadID < 0) { 
            return;
        }
        ProfilingFile.open(
            ".logs/" + OutputFileName() + "_" + std::to_string(RunID)
            + "_t" + std::to_string(ThreadID) + ".profile"  
        );
        ProfilingFile << "function,number_of_calls,total_time[ms]\n"; 
        
        auto Profiles = GetFunctionProfiles();
        for(auto const& [key, val] : Profiles)
        {
            ProfilingFile << key << ","; 
            ProfilingFile << val.calls << ","; 
            ProfilingFile << val.time << "\n";
        }
        ProfilingFile.close();
#endif
    }

    void EventManager::CreateTuples()
    {
        StartFunctionProfile();
        auto AnalysisManager = G4AnalysisManager::Instance();
#ifdef ARTIE_GEANT_10
#else
        AnalysisManager->SetDefaultFileType("root");
#endif
        AnalysisManager->SetVerboseLevel(0);
        AnalysisManager->SetNtupleMerging(true);

#ifdef ARTIE_YAML
        G4int index = GetIndex("Configuration");
        AnalysisManager->CreateNtuple("Configuration", "Configuration");
        AnalysisManager->CreateNtupleSColumn("category");
        AnalysisManager->CreateNtupleSColumn("parameter");
        AnalysisManager->CreateNtupleSColumn("value");
        AnalysisManager->FinishNtuple(index);
        mSavedParameters = false;
#endif
        if(SaveParticleMaps())
        {
            G4int index = GetIndex("ParticleMaps");
            AnalysisManager->CreateNtuple("ParticleMaps", "ParticleMaps");
            AnalysisManager->CreateNtupleIColumn("event");
            AnalysisManager->CreateNtupleIColumn("track_id");
            AnalysisManager->CreateNtupleSColumn("particle");
            AnalysisManager->CreateNtupleIColumn("pdg");
            AnalysisManager->CreateNtupleIColumn("parent_track_id");
            AnalysisManager->CreateNtupleIColumn("ancestor_track_id");
            AnalysisManager->FinishNtuple(index);
        }
        if(SavePrimaryInfo())
        {
            G4int index = GetIndex("Primaries");
            AnalysisManager->CreateNtuple("Primaries", "Primaries");
            AnalysisManager->CreateNtupleIColumn("event");
            AnalysisManager->CreateNtupleIColumn("track_id");
            AnalysisManager->CreateNtupleSColumn("particle");
            AnalysisManager->CreateNtupleIColumn("pdg");
            AnalysisManager->CreateNtupleDColumn("init_t");
            AnalysisManager->CreateNtupleDColumn("init_x");
            AnalysisManager->CreateNtupleDColumn("init_y");
            AnalysisManager->CreateNtupleDColumn("init_z");
            AnalysisManager->CreateNtupleDColumn("init_energy");
            AnalysisManager->CreateNtupleDColumn("final_t");
            AnalysisManager->CreateNtupleDColumn("final_x");
            AnalysisManager->CreateNtupleDColumn("final_y");
            AnalysisManager->CreateNtupleDColumn("final_z");
            AnalysisManager->CreateNtupleDColumn("final_energy");
            AnalysisManager->CreateNtupleIColumn("number_of_daughters");
            AnalysisManager->CreateNtupleDColumn("total_daughter_init_energy");
            AnalysisManager->CreateNtupleDColumn("total_daughter_final_energy");
            AnalysisManager->CreateNtupleDColumn("total_daughter_edep");
            AnalysisManager->CreateNtupleIColumn("number_of_edeps");
            AnalysisManager->CreateNtupleDColumn("total_edep");
            AnalysisManager->CreateNtupleIColumn("number_of_hits");
            AnalysisManager->FinishNtuple(index);
        }
        if(SaveParticleInfo())
        {
            G4int index = GetIndex("Particles");
            AnalysisManager->CreateNtuple("Particles", "Particles");
            AnalysisManager->CreateNtupleIColumn("event");
            AnalysisManager->CreateNtupleIColumn("track_id");
            AnalysisManager->CreateNtupleIColumn("pdg");
            AnalysisManager->CreateNtupleIColumn("parent_track_id");
            AnalysisManager->CreateNtupleSColumn("creation_process");
            AnalysisManager->CreateNtupleDColumn("global_init_time");
            AnalysisManager->CreateNtupleDColumn("local_init_time");
            AnalysisManager->CreateNtupleDColumn("init_x");
            AnalysisManager->CreateNtupleDColumn("init_y");
            AnalysisManager->CreateNtupleDColumn("init_z");
            AnalysisManager->CreateNtupleDColumn("init_energy");
            AnalysisManager->CreateNtupleDColumn("init_px");
            AnalysisManager->CreateNtupleDColumn("init_py");
            AnalysisManager->CreateNtupleDColumn("init_pz");
            AnalysisManager->FinishNtuple(index);
        }
        if(SaveEnergyDeposits())
        {
            G4int index = GetIndex("EnergyDeposits");
            AnalysisManager->CreateNtuple("EnergyDeposits", "EnergyDeposits");
            AnalysisManager->CreateNtupleIColumn("event");
            AnalysisManager->CreateNtupleIColumn("track_id");
            AnalysisManager->CreateNtupleSColumn("volume");
            AnalysisManager->CreateNtupleSColumn("particle");
            AnalysisManager->CreateNtupleIColumn("pdg");
            AnalysisManager->CreateNtupleDColumn("local_time");
            AnalysisManager->CreateNtupleDColumn("global_time");
            AnalysisManager->CreateNtupleDColumn("init_x");
            AnalysisManager->CreateNtupleDColumn("init_y");
            AnalysisManager->CreateNtupleDColumn("init_z");
            AnalysisManager->CreateNtupleDColumn("final_x");
            AnalysisManager->CreateNtupleDColumn("final_y");
            AnalysisManager->CreateNtupleDColumn("final_z");
            AnalysisManager->CreateNtupleDColumn("edep");
            AnalysisManager->FinishNtuple(index);
        }
        if(SaveProfileData())
        {
            G4int index = GetIndex("ProfileEventData");
            AnalysisManager->CreateNtuple("ProfileEventData", "ProfileEventData");
            AnalysisManager->CreateNtupleSColumn("name");
            AnalysisManager->CreateNtupleSColumn("volume");
            AnalysisManager->CreateNtupleIColumn("exit");
            AnalysisManager->CreateNtupleDColumn("x");
            AnalysisManager->CreateNtupleDColumn("y");
            AnalysisManager->FinishNtuple(index);
        }
        EndFunctionProfile("CreateTuples");
    }

    void EventManager::FillParticleMaps(G4int EventID)
    {
        if(!SaveParticleMaps()) {
            return;
        }
        StartFunctionProfile();

        auto AnalysisManager = G4AnalysisManager::Instance();
        G4int index = GetIndex("ParticleMaps");
        for(auto const& [key, val] : mParticleName)
        {
            AnalysisManager->FillNtupleIColumn(index, 0, EventID);
            AnalysisManager->FillNtupleIColumn(index, 1, key);
            AnalysisManager->FillNtupleSColumn(index, 2, GetParticleName(key));
            AnalysisManager->FillNtupleIColumn(index, 3, GetParticlePDG(key));
            AnalysisManager->FillNtupleIColumn(index, 4, GetParticleParentTrackID(key));
            AnalysisManager->FillNtupleIColumn(index, 5, GetParticleAncestorTrackID(key));
            AnalysisManager->AddNtupleRow(index);
        }
        EndFunctionProfile("FillParticleMaps");
    }
    void EventManager::FillPrimaryInfo(G4int EventID)
    {
        if(!SavePrimaryInfo()) {
            return;
        }
        StartFunctionProfile();

        auto AnalysisManager = G4AnalysisManager::Instance();
        G4int index = GetIndex("Primaries");
        for(size_t ii = 0; ii < mPrimaryData.size(); ii++)
        {
            AnalysisManager->FillNtupleIColumn(index, 0, EventID);
            AnalysisManager->FillNtupleIColumn(index, 1, mPrimaryData[ii].track_id);
            AnalysisManager->FillNtupleSColumn(index, 2, mPrimaryData[ii].name);
            AnalysisManager->FillNtupleIColumn(index, 3, mPrimaryData[ii].pdg);
            AnalysisManager->FillNtupleDColumn(index, 4, mPrimaryData[ii].init_t);
            AnalysisManager->FillNtupleDColumn(index, 5, mPrimaryData[ii].init_x);
            AnalysisManager->FillNtupleDColumn(index, 6, mPrimaryData[ii].init_y);
            AnalysisManager->FillNtupleDColumn(index, 7, mPrimaryData[ii].init_z);
            AnalysisManager->FillNtupleDColumn(index, 8, mPrimaryData[ii].init_energy);
            AnalysisManager->FillNtupleDColumn(index, 9, mPrimaryData[ii].final_t);
            AnalysisManager->FillNtupleDColumn(index, 10, mPrimaryData[ii].final_x);
            AnalysisManager->FillNtupleDColumn(index, 11, mPrimaryData[ii].final_y);
            AnalysisManager->FillNtupleDColumn(index, 12, mPrimaryData[ii].final_z);
            AnalysisManager->FillNtupleDColumn(index, 13, mPrimaryData[ii].final_energy);
            AnalysisManager->FillNtupleIColumn(index, 14, mPrimaryData[ii].number_of_daughters);
            AnalysisManager->FillNtupleDColumn(index, 15, mPrimaryData[ii].total_daughter_init_energy);
            AnalysisManager->FillNtupleDColumn(index, 16, mPrimaryData[ii].total_daughter_final_energy);
            AnalysisManager->FillNtupleDColumn(index, 17, mPrimaryData[ii].total_daughter_edep);
            AnalysisManager->FillNtupleIColumn(index, 18, mPrimaryData[ii].number_of_edeps);
            AnalysisManager->FillNtupleDColumn(index, 19, mPrimaryData[ii].total_edep);
            AnalysisManager->FillNtupleIColumn(index, 20, mPrimaryData[ii].number_of_hits);
            AnalysisManager->AddNtupleRow(index);
        }
        EndFunctionProfile("FillPrimaryInfo");
    }
    void EventManager::FillParticleInfo(G4int EventID)
    {
        StartFunctionProfile();
        if(!SaveParticleInfo()) {
            return;
        }
        auto AnalysisManager = G4AnalysisManager::Instance();
        G4int index = GetIndex("Particles");
        for(size_t ii = 0; ii < mParticles.size(); ii++)
        {
            AnalysisManager->FillNtupleIColumn(index, 0, EventID);
            AnalysisManager->FillNtupleIColumn(index, 1, mParticles[ii].TrackID());
            AnalysisManager->FillNtupleIColumn(index, 2, mParticles[ii].PDG());
            AnalysisManager->FillNtupleIColumn(index, 3, mParticles[ii].ParentTrackID());
            AnalysisManager->FillNtupleSColumn(index, 4, mParticles[ii].CreationProcess());
            AnalysisManager->FillNtupleDColumn(index, 5, mParticles[ii].GetGlobalCreationTime());
            AnalysisManager->FillNtupleDColumn(index, 6, mParticles[ii].GetT());
            AnalysisManager->FillNtupleDColumn(index, 7, mParticles[ii].GetX());
            AnalysisManager->FillNtupleDColumn(index, 8, mParticles[ii].GetY());
            AnalysisManager->FillNtupleDColumn(index, 9, mParticles[ii].GetZ());
            AnalysisManager->FillNtupleDColumn(index, 10, mParticles[ii].GetE());
            AnalysisManager->FillNtupleDColumn(index, 11, mParticles[ii].GetPx());
            AnalysisManager->FillNtupleDColumn(index, 12, mParticles[ii].GetPy());
            AnalysisManager->FillNtupleDColumn(index, 13, mParticles[ii].GetPz());
            AnalysisManager->AddNtupleRow(index);
        }
        EndFunctionProfile("FillParticleInfo");
    }
    void EventManager::FillEnergyDeposits(G4int EventID)
    {
        StartFunctionProfile();
        if(!SaveEnergyDeposits()) {
            return;
        } 
        auto AnalysisManager = G4AnalysisManager::Instance();
        G4int index = GetIndex("EnergyDeposits");
        for(size_t ii = 0; ii < mEnergyDeposits.size(); ii++)
        {
            AnalysisManager->FillNtupleIColumn(index, 0, EventID);
            AnalysisManager->FillNtupleIColumn(index, 1, mEnergyDeposits[ii].track_id);
            AnalysisManager->FillNtupleSColumn(index, 2, mEnergyDeposits[ii].volume);
            AnalysisManager->FillNtupleSColumn(index, 3, mEnergyDeposits[ii].name);
            AnalysisManager->FillNtupleIColumn(index, 4, mEnergyDeposits[ii].pdg);
            AnalysisManager->FillNtupleDColumn(index, 5, mEnergyDeposits[ii].local_time);
            AnalysisManager->FillNtupleDColumn(index, 6, mEnergyDeposits[ii].global_time);
            AnalysisManager->FillNtupleDColumn(index, 7, mEnergyDeposits[ii].pre_step_position[0]);
            AnalysisManager->FillNtupleDColumn(index, 8, mEnergyDeposits[ii].pre_step_position[1]);
            AnalysisManager->FillNtupleDColumn(index, 9, mEnergyDeposits[ii].pre_step_position[2]);
            AnalysisManager->FillNtupleDColumn(index, 10, mEnergyDeposits[ii].post_step_position[0]);
            AnalysisManager->FillNtupleDColumn(index, 11, mEnergyDeposits[ii].post_step_position[1]);
            AnalysisManager->FillNtupleDColumn(index, 12, mEnergyDeposits[ii].post_step_position[2]);
            AnalysisManager->FillNtupleDColumn(index, 13, mEnergyDeposits[ii].energy);
            AnalysisManager->AddNtupleRow(index);
        }
        EndFunctionProfile("FillEnergyDeposits");
    }
    void EventManager::AddParticleMapsFromTrack(const G4Track* track)
    {
        StartFunctionProfile();
        G4int track_id = track->GetTrackID();
        G4String particle_name = track->GetParticleDefinition()->GetParticleName();
        G4int particle_pdg = track->GetParticleDefinition()->GetPDGEncoding();
        G4int parent_track_id = track->GetParentID();
        
        const G4VProcess* creator_process = track->GetCreatorProcess();
        G4String process = "none";
        if(creator_process) {
            process = creator_process->GetProcessName();
        }

        // Save Particle Info    
        AddParticleName(track_id, particle_name);
        AddParticlePDG(track_id, particle_pdg);
        AddParticleParentTrackID(track_id, parent_track_id);

        // Add ancestor info
        if (parent_track_id == 0) {
            AddParticleAncestorTrackID(track_id, 0);
        }
        else if (GetParticleParentTrackID(parent_track_id) == 0) {
            AddParticleAncestorTrackID(track_id, parent_track_id);
        }
        else {
            AddParticleAncestorTrackID(track_id, GetParticleAncestorTrackID(parent_track_id));
        }

        EndFunctionProfile("AddParticleMapsFromTrack");
    }

    void EventManager::AddPrimaryInfoFromTrackBegin(const G4Track* track)
    {
        StartFunctionProfile();
        G4String particle_name = track->GetParticleDefinition()->GetParticleName();
        G4int particle_pdg = track->GetParticleDefinition()->GetPDGEncoding();
        G4ThreeVector particle_position = track->GetVertexPosition();

        G4double global_time = track->GetGlobalTime();
        G4double kinetic_energy = track->GetKineticEnergy();
        G4int track_id = track->GetTrackID();
        G4int parent_track_id = track->GetParentID();

        const G4VProcess* creator_process = track->GetCreatorProcess();
        G4String process = "none";
        if(creator_process) {
            process = creator_process->GetProcessName();
        }

        if(parent_track_id == 0)
        {
            AddPrimaryData(
                PrimaryData(
                    particle_name, track_id, particle_pdg,
                    global_time, particle_position[0], particle_position[1],
                    particle_position[2], kinetic_energy 
                )
            );
        }
        else
        {
            GetPrimaryData(GetParticleAncestorTrackID(track_id)).total_daughter_init_energy += kinetic_energy;
            GetPrimaryData(GetParticleAncestorTrackID(track_id)).number_of_daughters += 1;
        }
        EndFunctionProfile("AddPrimaryInfoFromTrackBegin");
    }

    void EventManager::AddPrimaryInfoFromTrackEnd(const G4Track* track)
    {
        StartFunctionProfile();
        G4String particle_name = track->GetParticleDefinition()->GetParticleName();
        G4ThreeVector particle_position = track->GetVertexPosition();
        G4double kinetic_energy = track->GetKineticEnergy();
        G4int track_id = track->GetTrackID();
        G4double global_time = track->GetGlobalTime();
        
        if(GetParticleParentTrackID(track_id) == 0)
        {
            GetPrimaryData(track_id).AddFinalTrackData(
                global_time, particle_position[0], particle_position[1],
                particle_position[2], kinetic_energy
            );
        }
        else {
            GetPrimaryData(GetParticleAncestorTrackID(track_id)).total_daughter_final_energy += kinetic_energy;
        }
        EndFunctionProfile("AddPrimaryInfoFromTrackEnd");
    }
    void EventManager::AddParticleInfoFromTrackBegin(const G4Track* track)
    {
        StartFunctionProfile();
        G4int track_id = track->GetTrackID();
        G4int particle_pdg = track->GetParticleDefinition()->GetPDGEncoding();
        G4int parent_track_id = track->GetParentID();
        const G4VProcess* creator_process = track->GetCreatorProcess();
        G4String process = "none";
        if(creator_process) {
            process = creator_process->GetProcessName();
        }
        G4double global_time = track->GetGlobalTime();
        G4double local_time = track->GetLocalTime();
        G4ThreeVector particle_position = track->GetVertexPosition();
        G4double kinetic_energy = track->GetKineticEnergy();
        G4ThreeVector particle_momentum = track->GetMomentum();

        mParticles.emplace_back(
            Particle(
                track_id, 
                particle_pdg,
                parent_track_id,
                process,
                global_time,
                local_time,
                particle_position,
                kinetic_energy,
                particle_momentum
            )
        );
        AddParticleTrackID(track_id, mParticles.size()-1);
        EndFunctionProfile("AddParticleInfoFromTrackBegin");
    }
    void EventManager::AddParticleInfoFromTrackEnd(const G4Track* track)
    {
        StartFunctionProfile();
        // G4int           track_id = track->GetTrackID();
        // G4double        local_time = track->GetLocalTime();
        // G4ThreeVector   end_position = track->GetStep()->GetPostStepPoint()->GetPosition(); //track->GetPosition();
        // G4ThreeVector   end_momentum = track->GetMomentum();
        // G4double        kinetic_energy = track->GetKineticEnergy();

        // mParticles[GetParticleTrackID(track_id)].AddTrajectoryPoint(
        //     local_time, end_position, kinetic_energy, end_momentum
        // );

        EndFunctionProfile("AddParticleInfoFromTrackEnd");
    }
    void EventManager::AddParticleInfoFromStep(const G4Step *step)
    {
        StartFunctionProfile();
        G4VPhysicalVolume* physicalVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
        G4LogicalVolume *logicalVolume = physicalVolume->GetLogicalVolume();
        G4Track* track = step->GetTrack();
        G4StepPoint* pre_step_point = step->GetPreStepPoint();
        G4StepPoint* post_step_point = step->GetPostStepPoint();
        const G4VProcess* pre_process = pre_step_point->GetProcessDefinedStep();
        const G4VProcess* post_process = post_step_point->GetProcessDefinedStep();

        G4int           event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        G4String        particle_name = step->GetTrack()->GetParticleDefinition()->GetParticleName();
        G4int           particle_pdg = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
        G4ThreeVector   particleStartPosition = step->GetPreStepPoint()->GetPosition();
        G4double        local_time = track->GetLocalTime();
        G4double        global_time = track->GetGlobalTime();
        G4double        pre_step_energy = pre_step_point->GetKineticEnergy();
        G4double        post_step_energy = post_step_point->GetKineticEnergy();
        G4int           track_status = (int)track->GetTrackStatus();
        G4int           track_id = track->GetTrackID();
        G4int           parent_track_id = track->GetParentID();
        G4ThreeVector   pre_step_position = pre_step_point->GetPosition();
        G4ThreeVector   post_step_position = post_step_point->GetPosition();
        G4ThreeVector   pre_step_momentum = pre_step_point->GetMomentum();
        G4ThreeVector   post_step_momentum = post_step_point->GetMomentum();
        G4double        edep = step->GetTotalEnergyDeposit();
        G4double        kinetic_energy = track->GetKineticEnergy();
        if(GetParticleParentTrackID(track_id) == 0)
        {
            if(parent_track_id == 0) {
                GetPrimaryData(track_id).number_of_edeps += 1;
                GetPrimaryData(track_id).total_edep += edep;
            }
            else {
                GetPrimaryData(GetParticleAncestorTrackID(track_id)).total_daughter_edep += edep;
            }
        }    

        mParticles[GetParticleTrackID(track_id)].AddTrajectoryPoint(
            local_time, pre_step_position, kinetic_energy, pre_step_momentum
        );

        G4StepStatus pre_step_status = pre_step_point->GetStepStatus();
        G4StepStatus post_step_status = post_step_point->GetStepStatus();
        if (post_step_status == fWorldBoundary)
        {
            mParticles[GetParticleTrackID(track_id)].AddTrajectoryPoint(
                local_time, post_step_position, kinetic_energy, post_step_momentum
            );
        }
        
        EndFunctionProfile("AddParticleInfoFromStep");
    }

    void EventManager::AddEnergyDepositInfoFromStep(const G4Step* step)
    {
        StartFunctionProfile();
        G4VPhysicalVolume* physicalVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
        G4LogicalVolume *logicalVolume = physicalVolume->GetLogicalVolume();
        G4String volume = logicalVolume->GetName();

        G4Track* track = step->GetTrack();
        G4StepPoint* pre_step_point = step->GetPreStepPoint();
        G4StepPoint* post_step_point = step->GetPostStepPoint();
        const G4VProcess* pre_process = pre_step_point->GetProcessDefinedStep();
        const G4VProcess* post_process = post_step_point->GetProcessDefinedStep();

        G4int           event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        G4String        particle_name = step->GetTrack()->GetParticleDefinition()->GetParticleName();
        G4int           particle_pdg = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
        G4ThreeVector   particleStartPosition = step->GetPreStepPoint()->GetPosition();
        G4double        local_time = track->GetLocalTime();
        G4double        global_time = track->GetGlobalTime();
        G4double        pre_step_energy = pre_step_point->GetKineticEnergy();
        G4double        post_step_energy = post_step_point->GetKineticEnergy();
        G4int           track_status = (int)track->GetTrackStatus();
        G4int           track_id = track->GetTrackID();
        G4int           parent_track_id = track->GetParentID();
        G4ThreeVector   pre_step_position = pre_step_point->GetPosition();
        G4ThreeVector   post_step_position = post_step_point->GetPosition();
        G4double        edep = step->GetTotalEnergyDeposit();
        if (edep == 0.0) {
            return;
        }
        mEnergyDeposits.emplace_back(
            EnergyDeposit(
                track_id, volume, particle_name,
                particle_pdg, local_time, global_time,
                pre_step_position, post_step_position,
                edep
            )
        );
        EndFunctionProfile("AddEnergyDepositInfoFromStep");
    }
    void EventManager::EvaluateRunBegin()
    {
        for(size_t ii = 0; ii < mAnalysisRunBeginFunctions.size(); ii++)
        {
            mAnalysisRunBeginFunctions[ii]();
        }
    }
    void EventManager::EvaluateRunEnd()
    {
        for(size_t ii = 0; ii < mAnalysisRunEndFunctions.size(); ii++)
        {
            mAnalysisRunEndFunctions[ii]();
        }
    }
    void EventManager::EvaluateEventBegin()
    {
        for(size_t ii = 0; ii < mAnalysisEventBeginFunctions.size(); ii++)
        {
            mAnalysisEventBeginFunctions[ii]();
        }
    }
    void EventManager::EvaluateEventEnd()
    {
        for(size_t ii = 0; ii < mAnalysisEventEndFunctions.size(); ii++)
        {
            mAnalysisEventEndFunctions[ii]();
        }
    }
}