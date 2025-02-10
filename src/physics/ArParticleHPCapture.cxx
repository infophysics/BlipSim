#include "ArParticleHPCapture.hh"
#include "ArParticleHPCaptureFS.hh"

#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleHPDeExGammas.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Threading.hh"

namespace Blip
{
  ArParticleHPCapture::ArParticleHPCapture()
   :G4HadronicInteraction("NeutronHPCapture")
  ,theCapture(NULL)
  ,numEle(0)
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*MeV );
/*
//    G4cout << "Capture : start of construction!!!!!!!!"<<G4endl;
    if(!std::getenv("G4NEUTRONHPDATA"))
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
    dirName = std::getenv("G4NEUTRONHPDATA");
    G4String tString = "/Capture";
    dirName = dirName + tString;
    numEle = G4Element::GetNumberOfElements();
//    G4cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
//    G4cout <<"Disname="<<dirName<<" numEle="<<numEle<<G4endl;
    //theCapture = new G4ParticleHPChannel[numEle];
//    G4cout <<"G4ParticleHPChannel constructed"<<G4endl;
    ArParticleHPCaptureFS * theFS = new ArParticleHPCaptureFS;
    //for (G4int i=0; i<numEle; i++)
    //{
//  //    G4cout << "initializing theCapture "<<i<<" "<< numEle<<G4endl;
    //  theCapture[i].Init((*(G4Element::GetElementTable()))[i], dirName);
    //  theCapture[i].Register(theFS);
    //}
    for ( G4int i = 0 ; i < numEle ; i++ )
    {
       theCapture.push_back( new G4ParticleHPChannel );
       (*theCapture[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
       (*theCapture[i]).Register(theFS);
    }
    delete theFS;
//    G4cout << "-------------------------------------------------"<<G4endl;
//    G4cout << "Leaving ArParticleHPCapture::ArParticleHPCapture"<<G4endl;
*/
  }

  ArParticleHPCapture::~ArParticleHPCapture()
  {
    //delete [] theCapture;
    //vector is shared, only master deletes
    if ( ! G4Threading::IsWorkerThread() ) {
        if ( theCapture != NULL ) {
            for ( std::vector<G4ParticleHPChannel*>::iterator
                ite = theCapture->begin() ; ite != theCapture->end() ; ite++ ) {
                delete *ite;
            }
            theCapture->clear();
        }
    }

      //<--for(int i=0;i<100;i++){
      //<--    std::cout<<"I am ArParticleHPCapture"<<std::endl;
      //<--}
  }

  #include "G4ParticleHPThermalBoost.hh"
  G4HadFinalState * ArParticleHPCapture::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus )
  {

    //if ( numEle < (G4int)G4Element::GetNumberOfElements() ) addChannelForNewElement();

    G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();
    if(std::getenv("NeutronHPCapture")) G4cout <<" ####### ArParticleHPCapture called"<<G4endl;
    const G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    G4int index = theMaterial->GetElement(0)->GetIndex();
    if(n!=1)
    {
      G4double* xSec = new G4double[n];
      G4double sum=0;
      G4int i;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;
      G4ParticleHPThermalBoost aThermalE;
      for (i=0; i<n; i++)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        //xSec[i] = theCapture[index].GetXsec(aThermalE.GetThermalEnergy(aTrack,
        xSec[i] = ((*theCapture)[index])->GetXsec(aThermalE.GetThermalEnergy(aTrack,
  		                                                     theMaterial->GetElement(i),
  								     theMaterial->GetTemperature()));
        xSec[i] *= rWeight;
        sum+=xSec[i];
      }
      G4double random = G4UniformRand();
      G4double running = 0;
      for (i=0; i<n; i++)
      {
        running += xSec[i];
        index = theMaterial->GetElement(i)->GetIndex();
        //if(random<=running/sum) break;
        if( sum == 0 || random <= running/sum ) break;
      }
      if(i==n) i=std::max(0, n-1);
      delete [] xSec;
    }

    //return theCapture[index].ApplyYourself(aTrack);
    //G4HadFinalState* result = theCapture[index].ApplyYourself(aTrack);
    G4HadFinalState* result = ((*theCapture)[index])->ApplyYourself(aTrack);

    //Overwrite target parameters
    aNucleus.SetParameters(G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
    const G4Element* target_element = (*G4Element::GetElementTable())[index];
    const G4Isotope* target_isotope=NULL;
    G4int iele = target_element->GetNumberOfIsotopes();
    for ( G4int j = 0 ; j != iele ; j++ ) {
       target_isotope=target_element->GetIsotope( j );
       if ( target_isotope->GetN() == G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA() ) break;
    }
    //G4cout << "Target Material of this reaction is " << theMaterial->GetName() << G4endl;
    //G4cout << "Target Element of this reaction is " << target_element->GetName() << G4endl;
    //G4cout << "Target Isotope of this reaction is " << target_isotope->GetName() << G4endl;
    aNucleus.SetIsotope( target_isotope );

    G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();
    return result;
  }

const std::pair<G4double, G4double> ArParticleHPCapture::GetFatalEnergyCheckLevels() const
{
   //return std::pair<G4double, G4double>(10*perCent,10*GeV);
   return std::pair<G4double, G4double>(10*perCent,DBL_MAX);
}

/*
void ArParticleHPCapture::addChannelForNewElement()
{
   ArParticleHPCaptureFS* theFS = new ArParticleHPCaptureFS;
   for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ )
   {
      G4cout << "ArParticleHPCapture Prepairing Data for the new element of " << (*(G4Element::GetElementTable()))[i]->GetName() << G4endl;
      theCapture.push_back( new G4ParticleHPChannel );
      (*theCapture[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
      (*theCapture[i]).Register(theFS);
   }
   delete theFS;
   numEle = (G4int)G4Element::GetNumberOfElements();
}
*/

G4int ArParticleHPCapture::GetVerboseLevel() const
{
   return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}
void ArParticleHPCapture::SetVerboseLevel( G4int newValue )
{
   G4ParticleHPManager::GetInstance()->SetVerboseLevel(newValue);
}

void ArParticleHPCapture::BuildPhysicsTable(const G4ParticleDefinition&)
{

   G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();

   theCapture = hpmanager->GetCaptureFinalStates();

   if ( G4Threading::IsMasterThread() ) {

      if ( theCapture == NULL ) theCapture = new std::vector<G4ParticleHPChannel*>;

      if ( numEle == (G4int)G4Element::GetNumberOfElements() ) return;

      if ( theCapture->size() == G4Element::GetNumberOfElements() ) {
         numEle = G4Element::GetNumberOfElements();
         return;
      }

      if ( !std::getenv("G4NEUTRONHPDATA") )
          throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
      dirName = std::getenv("G4NEUTRONHPDATA");
      G4String tString = "/Capture";
      dirName = dirName + tString;

      G4ParticleHPCaptureFS * theFS = new G4ParticleHPCaptureFS;
      ArParticleHPCaptureFS * theArFS = new ArParticleHPCaptureFS;
      for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ )
      {
	      theCapture->push_back( new G4ParticleHPChannel );
        ((*theCapture)[i])->Init((*(G4Element::GetElementTable()))[i], dirName);
        if((*(G4Element::GetElementTable()))[i]->GetZ() == 18) {
         	((*theCapture)[i])->Register(theArFS);
         	std::cout<<"======= use new Argon Capture ======="<<std::endl;
        }
        else ((*theCapture)[i])->Register(theFS);
      }
      delete theFS;
      delete theArFS;
      hpmanager->RegisterCaptureFinalStates( theCapture );
   }
   numEle = G4Element::GetNumberOfElements();
}

void ArParticleHPCapture::ModelDescription(std::ostream& outFile) const
{
   outFile << "High Precision model based on Evaluated Nuclear Data Files (ENDF) for radiative capture reaction of neutrons below 20MeV\n";
}
}