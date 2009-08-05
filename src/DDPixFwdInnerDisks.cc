/* 
== CMS Forward Pixels Geometry ==

 
Author:  Pratima Jindal, Purdue University Calumet
         July 2009


  Algorithm for placing blades on inner disks for Phase 1 Pixel Forward

*/

#include <cmath>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DetectorDescription/Base/interface/DDRotationMatrix.h"
#include "DetectorDescription/Base/interface/DDutils.h"
#include "DetectorDescription/Core/interface/DDPosPart.h"
#include "DetectorDescription/Core/interface/DDLogicalPart.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "DetectorDescription/Core/interface/DDMaterial.h"
#include "DetectorDescription/Core/interface/DDCurrentNamespace.h"
#include "DetectorDescription/Core/interface/DDSplit.h"
#include "DetectorDescription/Core/interface/DDConstant.h"
#include "SLHCUpgradeSimulations/Geometry/interface/DDPixFwdInnerDisks.h"
#include "CLHEP/Vector/RotationInterfaces.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

  // -- Input geometry parameters :  -----------------------------------------------------

const int DDPixFwdInnerDisks::nBlades = 24;            // Number of blades
const double DDPixFwdInnerDisks::bladeAngle = 20.*deg;    // Angle of blade rotation around its axis
const double DDPixFwdInnerDisks::zPlane = 4.5*mm;             // Common shift in Z for all blades (with respect to disk center plane)
const double DDPixFwdInnerDisks::bladeZShift = 0.*mm;     // Shift in Z between the axes of two adjacent blades
  
const double DDPixFwdInnerDisks::ancorRadius =  80.6*mm; // Distance from beam line to ancor point defining center of "blade frame"
  
  

std::map<std::string, int> DDPixFwdInnerDisks::copyNumbers;

// -- Constructors & Destructor :  -------------------------------------------------------

DDPixFwdInnerDisks::DDPixFwdInnerDisks() {}
DDPixFwdInnerDisks::~DDPixFwdInnerDisks() {}

// Initialization :  ---------------------------------------------------------------------

void DDPixFwdInnerDisks::initialize(const DDNumericArguments & nArgs,
				const DDVectorArguments & vArgs,
				const DDMapArguments & ,
				const DDStringArguments & sArgs,
				const DDStringVectorArguments & ) {
				  	
  if ( nArgs.find("Endcap") != nArgs.end() ) {
    endcap = nArgs["Endcap"];
  } else {
    endcap = 1.;
  }

  if ( sArgs.find("FlagString") != sArgs.end() ) {
    flagString = sArgs["FlagString"];
    flagSelector = sArgs["FlagSelector"];
  } else {
    flagString = "YYYYYYYYYYYYYYYYYYYYYYYY";
    flagSelector = "Y";
  }

  if ( sArgs.find("Child") != sArgs.end() ) {
    childName   = sArgs["Child"];
  } else {
    childName   = "";
  }

  if ( vArgs.find("ChildTranslation") != vArgs.end() ) {
    childTranslationVector = vArgs["ChildTranslation"];
  } else {
    childTranslationVector = std::vector<double>(3, 0.);
  }

  if ( sArgs.find("ChildRotation") != sArgs.end() ) {
    childRotationName = sArgs["ChildRotation"];
  } else {
    childRotationName = "";
  }

  idNameSpace = DDCurrentNamespace::ns();
}
  
// Execution :  --------------------------------------------------------------------------

void DDPixFwdInnerDisks::execute() {

  // -- Compute Nipple parameters if not already computed :

	if (childName == "") return;
  
  // -- Signed versions of blade angle and z-shift :
  
  double effBladeAngle = - endcap * bladeAngle;
  double effBladeZShift = endcap * bladeZShift;
  
  // -- Names of mother and child volumes :

  DDName mother = parent().name();
  DDName child(DDSplit(childName).first, DDSplit(childName).second);
  
  // -- Get translation and rotation from "blade frame" to "child frame", if any :
  
  HepRotation childRotMatrix = HepRotation();
  if (childRotationName != "") {
    DDRotation childRotation = DDRotation(DDName(DDSplit(childRotationName).first, DDSplit(childRotationName).second));
    // due to conversion to ROOT::Math::Rotation3D -- Michael Case
    DD3Vector x, y, z;
    childRotation.rotation()->GetComponents(x, y, z); // these are the orthonormal columns.
    HepRep3x3 tr(x.X(), y.X(), z.X(), x.Y(), y.Y(), z.Y(), x.Z(), y.Z(), z.Z());
    childRotMatrix = HepRotation(tr);
  }

  
  Hep3Vector childTranslation;

  	childTranslation = Hep3Vector(childTranslationVector[0],childTranslationVector[1],childTranslationVector[2]);
 
  // Create a matrix for rotation around blade axis (to "blade frame") :
  
  HepRotation bladeRotMatrix(Hep3Vector(0.,1.,0.), effBladeAngle);
  
  // Cycle over Phi positions, placing copies of the child volume :

  double deltaPhi = (360./nBlades)*deg;
  int nQuarter = nBlades/4;
  double zShiftMax = effBladeZShift*((nQuarter-1)/2.);

  for (int iBlade=0; iBlade < nBlades; iBlade++) {
    
    // check if this blade position should be skipped :
  	
  	if (flagString[iBlade] != flagSelector[0]) continue;
    int copy = issueCopyNumber();
    
    // calculate Phi and Z shift for this blade :

    double phi = (iBlade + 0.5) * deltaPhi - 90.*deg;
    int iQuarter = iBlade % nQuarter;
    double zShift = - zShiftMax + iQuarter * effBladeZShift;
    
    // compute rotation matrix from mother to blade frame :
    
    HepRotation* rotMatrix = new HepRotation(Hep3Vector(0.,0.,1.), phi);
    (*rotMatrix) *= bladeRotMatrix;
    
    // convert translation vector from blade frame to mother frame, and add Z shift :
    
    Hep3Vector translation = (*rotMatrix)(childTranslation + Hep3Vector(0., ancorRadius, 0.));
    translation += Hep3Vector(0., 0., zShift + zPlane);
    
    // create DDRotation for placing the child if not already existent :

    DDRotation rotation;   
    std::string rotstr = DDSplit(childName).first + int_to_string(copy);
    rotation = DDRotation(DDName(rotstr, idNameSpace));

    if (!rotation) {
      (*rotMatrix) *= childRotMatrix;
      DDRotationMatrix* temp = new DDRotationMatrix(rotMatrix->xx(), rotMatrix->xy(), rotMatrix->xz(),
						    rotMatrix->yx(), rotMatrix->yy(), rotMatrix->yz(),
						    rotMatrix->zx(), rotMatrix->zy(), rotMatrix->zz() );
      rotation = DDrot(DDName(rotstr, idNameSpace), temp);
    } else {
      delete rotMatrix;
    }

    // position the child :

    DDTranslation ddtran(translation.x(), translation.y(), translation.z());
    DDpos(child, mother, copy, ddtran, rotation);
    // LogDebug("PixelGeom") << "DDPixFwdInnerDisks: " << child << " Copy " << copy << " positioned in " << mother << " at " << translation << " with rotation " << rotation;
  }

  // End of cycle over Phi positions

}

// -- Helpers :  -------------------------------------------------------------------------

int DDPixFwdInnerDisks::issueCopyNumber() {
  if (copyNumbers.count(childName) == 0) copyNumbers[childName] = 0;
  return ++copyNumbers[childName];
}


