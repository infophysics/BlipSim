#include <iomanip>

#include <CLHEP/Units/SystemOfUnits.h>

// -- artg4tk includes
#include "MyG4HadronPhysicsQGSP_BERT_ArHP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"

//#include "G4HadronPhysicsQGSP_BERT_HP.hh"

/////////////////////////////////////////////////////////////////////////////
// The following change is the _only_ required changed to move from
// the non-extensible factory to the exensible factory.  All other changes
// relative to the "factory" example are there to demonstrate new features.
/////////////////////////////////////////////////////////////////////////////
//non-extensible:  #include "G4PhysListFactory.hh"
#include "G4PhysListFactoryAlt.hh"
/////////////////////////////////////////////////////////////////////////////
// headers needed to demonstrate new features
/////////////////////////////////////////////////////////////////////////////

// allow ourselves to extend the short names for physics ctor addition/replace
// along the same lines as EMX, EMY, etc
#include "G4PhysListRegistry.hh"

// allow ourselves to give the user extra info about available physics ctors
#include "G4PhysicsConstructorFactory.hh"

/////////////////////////////////////////////////////////////////////////////
// pull in a user defined physics list definition into the main program
// and register it with the factory (doesn't have to be the main program
// but the .o containing the declaration _must_ get linked/loaded)

#include "G4VModularPhysicsList.hh"

#include "G4PhysListStamper.hh"  // defines macro for factory registration
#include "MyQGSP_BERT_ArHP.hh"

// -- Register the physics list
namespace Blip
{
G4_DECLARE_PHYSLIST_FACTORY(MyQGSP_BERT_ArHP);

MyQGSP_BERT_ArHP::MyQGSP_BERT_ArHP(G4int ver)
{

  G4cout << "<<< Geant4 Physics List simulation engine: MyQGSP_BERT_ArHP"<<G4endl;
  G4cout <<G4endl<<G4endl;

  defaultCutValue = 0.7*CLHEP::mm;
  SetVerboseLevel(ver);

  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(ver) );

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  RegisterPhysics( new G4DecayPhysics(ver) );
  RegisterPhysics( new G4RadioactiveDecayPhysics(ver) );

  // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );

  // Hadron Physics
  RegisterPhysics( new MyG4HadronPhysicsQGSP_BERT_ArHP(ver));

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver));

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));

}

void MyQGSP_BERT_ArHP::SetCuts()
{
  if (verboseLevel >1){
    G4cout << "MyQGSP_BERT_ArHP::SetCuts:";
  }
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types

  SetCutsWithDefault();

  //Set proton cut value to 0 for producing low energy recoil nucleus
  SetCutValue(0, "proton");

}
}