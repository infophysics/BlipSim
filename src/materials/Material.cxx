/**
 * @file Material.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-30
 */
#include "Material.hh"

namespace Blip
{
    G4Material* CreateMaterial(G4String material_name, G4String volume_name)
    {
        G4Material* material;
        if(material_name == "gaseous_argon") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar");
        }
        else if(material_name == "gaseous_argon_300bar") {
            G4Isotope* Ar36 = new G4Isotope("Ar36", 18., 36, 35.9675*g/mole);
            G4Isotope* Ar38 = new G4Isotope("Ar38", 18., 38, 37.9627*g/mole);
            G4Isotope* Ar40 = new G4Isotope("Ar40", 18., 40, 39.9624*g/mole);

            G4Element* elAr  = new G4Element("Natural Argon", "Ar", 3);
            elAr->AddIsotope(Ar36,0.334*perCent);
            elAr->AddIsotope(Ar38,0.063*perCent);
            elAr->AddIsotope(Ar40,99.603*perCent);

            material = new G4Material("Gaseous Argon 300bar", 0.4837 * g/cm3, 1);
            material->AddElement(elAr, 1);
        }
        else if(material_name == "liquid_argon") {
            if(EventManager::GetEventManager()->UseG4Definition()) {
                material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
            }
            else {
                auto argon = new Argon("LAr");
                material = argon->GetMaterial();
            }
        }
        else if(material_name == "g4_liquid_argon") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
        }
        else if(material_name == "lithium_glass") {
            G4Element* elSi = new G4Element("Silicon", "Si", 14., 28.09*g/mole);
            G4Element* elO = new G4Element("Oxygen", "O", 8., 16.00*g/mole);
            G4Element* elMg = new G4Element("Magnesium", "Mg", 12., 24.31*g/mole);
            G4Element* elAl = new G4Element("Aluminum", "Al", 13., 26.98*g/mole);

            G4Material* SiO2 = new G4Material("SiO2", 2.2*g/cm3, 2);
            SiO2->AddElement(elSi, 1);
            SiO2->AddElement(elO, 2);

            G4Material* MgO = new G4Material("MgO", 3.58*g/cm3, 2);
            MgO->AddElement(elMg, 1);
            MgO->AddElement(elO, 1);

            G4Material* Al2O3 = new G4Material("Al2O3", 3.97*g/cm3, 2);
            Al2O3->AddElement(elAl, 2);
            Al2O3->AddElement(elO, 3);

            G4Isotope* lig_6 = new G4Isotope("Li6", 3., 6.);
            G4Isotope* lig_7 = new G4Isotope("Li7", 3., 7.);
            G4Element* li = new G4Element("Li", "Li", 2);
            li->AddIsotope(lig_6, 95.0);
            li->AddIsotope(lig_7, 5.0);

            G4Material* Li2O = new G4Material("Li2O", 2.01 * g/cm3, 2);
            Li2O->AddElement(li, 2);
            Li2O->AddElement(elO, 1);

            G4Isotope* Ce3Ion = new G4Isotope("Ce3Ion", 58, 140, 0.0*g/mole, 15.8708*keV);
            G4Element* ce = new G4Element("Ce", "Ce", 1);
            ce->AddIsotope(Ce3Ion, 100.0);

            G4Material* Ce2O3 = new G4Material("Ce2O3", 7.13*g/cm3, 2);
            Ce2O3->AddElement(ce, 2);
            Ce2O3->AddElement(elO, 3);

            material = new G4Material(
                "LiG",
                2.5 * g / cm3,
                5
            );
            material->AddMaterial(SiO2, 0.56);
            material->AddMaterial(MgO, 0.04);
            material->AddMaterial(Al2O3, 0.18);
            material->AddMaterial(Ce2O3, 0.04);
            material->AddMaterial(Li2O, 0.18);
        }
        else if(material_name == "borated_polyethylene") {
            G4NistManager* man = G4NistManager::Instance();
            material = new G4Material("Borated Polyethylene", 1.04*g/cm3, 2);
            G4Material* polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
            G4Material* boron = man->FindOrBuildMaterial("G4_B");
            material->AddMaterial(polyethylene, 0.95);
            material->AddMaterial(boron, 0.05);
        }
        else if(material_name == "graphite") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_GRAPHITE");
        }
        else if(material_name == "concrete") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");
        }
        else if(material_name == "kapton") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
        }
        else if(material_name == "polyurethane") {
            G4NistManager* man = G4NistManager::Instance();
            material = new G4Material("Polyurethane", 1.005*g/cm3, 4);
            G4Element* C = man->FindOrBuildElement("C");
            G4Element* H = man->FindOrBuildElement("H");
            G4Element* N = man->FindOrBuildElement("N");
            G4Element* O = man->FindOrBuildElement("O");
            material->AddElement(C, 3);
            material->AddElement(H, 8);
            material->AddElement(N, 2);
            material->AddElement(O, 1);
        }
        else if(material_name == "stainless_steel") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
        }
        else if(material_name == "rough_vacuum") 
        {
            auto air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
            material = new G4Material(
                "VacuumRough", air->GetDensity() * 0.1, 1
            );
            material->AddMaterial(air, 1.0);
        }
        else if(material_name == "high_vacuum") 
        {
            auto air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
            material = new G4Material(
                "VacuumHigh", air->GetDensity() * 1.0E-11, 1
            );
            material->AddMaterial(air, 1.0);
        }
        else if(material_name == "water") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
        }
        else if(material_name == "air") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
        }
        else if(material_name == "aluminum") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
        }
        else if(material_name == "sulfur") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_S");
        }
        else if(material_name == "carbon") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_C");
        }
        else if(material_name == "tungsten") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
        }
        else if(material_name == "brass") {
            material = G4NistManager::Instance()->FindOrBuildMaterial("G4_BRASS");
        }
        else
        {
            std::cerr << "ERROR! Material (" + material_name + ") not recognized!" << G4endl;
            exit(0);
        }
        material->SetName(volume_name + "_" + material_name);
        return material;
    }
}